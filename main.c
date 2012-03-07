/*
**
** A model of growth where each axon pathway is a simulated as a 
** series of segments between cells. Each cell tries to follow the path
** of its nearest neighbour, but if out of room searches with +-90 degrees
** of desired trajectory for a new node/cell in its path.
**
** Uses a 1 micron square grid: lat/y and lon/x, with fovea (0,0).
**
** Andrew Turpin
** Mon Nov  7 11:11:28 EST 2011
**
*/

#include <gsl/gsl_rng.h>
#include <gdk/gdk.h>
#include <glib.h>
#include <math.h>
#include <values.h>
#include <assert.h>
#include "setup.h"
#include "queue.h"
#include "density.h"
#include "types.h"

int debug = 0; 

#define PRINT_ENDPOINTS
#define PRINT_PATHS 100
//#define PRINT_OCT_PROFILE

   // macros for handling byte for count and thickness
//#define IS_ROOM(_c) (((_c)->thickness != UCHAR_MAX) && ( (_c)->count < (_c)->thickness))
//#define INC_COUNT(_c) do { (_c)->count += (((_c)->count) < UCHAR_MAX) ? 1 : 0; } while (0);
#define IS_ROOM(_c) ((_c)->count < (_c)->thickness)
#define INC_COUNT(_c) do { (_c)->count += 1; } while (0);

#define SCAN_POINTS_RADIUS (SIZE/2)   // size of search lookup table
PointD *scanPoints; // list of deltaX, deltaY, theta to use for searching grid
int scanPointLen;   // scanPoints[0..scanPointLen-1] are valid

extern Cell *cellBlock; // from setup.c
extern int numCells;    // from setup.c

/*
** Initialise scanPoints array
**    p.x = delta x from 0
**    p.y = delta y from 0
**    dist = angle of (x,y) from centre (radians)
** Elements are sorted by distance from (0,0)
**
** Sets scanPointLen.
*/
void
init_scanPoints() {
   scanPoints = (PointD *) malloc(sizeof(PointD)*(2*SCAN_POINTS_RADIUS+1)*(2*SCAN_POINTS_RADIUS+1));
   assert(scanPoints != NULL);

      // use dist temporarily for sorting
   int index = 0;
   for(int i = -NEW_PATH_RADIUS_LIMIT ; i <= +NEW_PATH_RADIUS_LIMIT ; i++)
      for(int j = -NEW_PATH_RADIUS_LIMIT ; j <= +NEW_PATH_RADIUS_LIMIT ; j++) {
         if (i == 0 && j == 0) continue;  // exclude (0,0)
         Point p = {i,j};
         Point po = {0,0};
         float dist = DIST(p,po);
         if (dist > NEW_PATH_RADIUS_LIMIT) continue;    // exclude corner points
         scanPoints[index].p.x = i;
         scanPoints[index].p.y = j;
         scanPoints[index].dist = dist;
         index++;
      }

   scanPointLen = index;
   qsort(scanPoints, scanPointLen, sizeof(PointD), cmp_PointD);

      // now replace dist with theta
   for(int i = 0 ; i < scanPointLen ; i++)
      scanPoints[i].dist = atan2(scanPoints[i].p.y, scanPoints[i].p.x);
}//init_scanPoints()

/*
** Return 1 if (x,y) is inside FOVEA_RADIUS, 0 otherwise
*/
int
in_fovea(int x, int y) {
   float distSqr = (x-SIZE/2)*(x-SIZE/2) + (y-SIZE/2)*(y-SIZE/2);
   return (distSqr <= FOVEA_RADIUS*FOVEA_RADIUS);
}//in_fovea()

/*
** Return 1 if a and b are on different sides of the raphe, 0 otherwise
*/
int
cross_raphe(Point a, Point b) {
    if (a.x > SIZE/2 || b.x > SIZE/2) return 0;  // there is no spoon (raphe)

    if (a.y < SIZE/2) 
        return (b.y > SIZE/2);

    return (b.y < SIZE/2);
}//cross_raphe()

/*
** Just print cell position and last point of path,
** or fll path if full=TRUE
*/
void print_path(Cell *c, char full) {
   if (full) {
      Node *p = c->path; 
      while (p != NULL) {
         printf("%6d %6d\n",p->c->p.x, p->c->p.y);
         p = p->next;
      }
      printf("-1 -1\n");
   } else {
      if (c->path != NULL) {
            // just print first and last
         Node *p = c->path;
         printf("%6d %6d ",p->c->p.x, p->c->p.y);
         while(p->next != NULL)
            p = p->next;

         printf("%6d %6d\n",p->c->p.x, p->c->p.y);
      }
   }
   fflush(stdout);
}//print_path()

/*
** Find closest cell in the direction of target from current (+-THETA_LIMIT)
** that is not target, and that has room to be part of a path.
**
** Use a spiral search (well, they're squares of increasing dist) outwards from
** current ignoring cells with flag==1, points outside theta range, and 
** count >= thickness.
*/
Cell *
findNewPath(Cell *current, Cell *target, Grid **grid) {
if (debug)printf("# current = %5d %5d ",current->p.x, current->p.y);
if (debug)printf(" wants %5d %5d ",target->p.x, target->p.y);

   double theta = atan2(target->p.y - current->p.y, target->p.x - current->p.x);

   target->flag = 1; // rule out current
   Cell *c, *minC = NULL;
   for(PointD *s = scanPoints ; (s < scanPoints + scanPointLen) && (minC == NULL) ; s++) {
      if (s->dist < theta - THETA_LIMIT) continue;    // outside theta range
      if (s->dist > theta + THETA_LIMIT) continue;
      if ((target->p.x < SIZE/2) && (target->p.y < SIZE/2) && (s->p.y > 0)) continue; // away from raphe, x < fov only
      if ((target->p.x < SIZE/2) && (target->p.y > SIZE/2) && (s->p.y < 0)) continue; // away from raphe, x < fov only
      int x = target->p.x + s->p.x;
      int y = target->p.y + s->p.y;
      if ((x < 0) || (x >= SIZE)) continue;        // off grid
      if ((y < 0) || (y >= SIZE)) continue;
      //if ((target->p.y < SIZE/2) && (y > SIZE/2)) continue; // no cross raphe
      //if ((target->p.y > SIZE/2) && (y < SIZE/2)) continue; // no cross raphe

      if ((c = grid[x][y].soma) == NULL)  {   // blanko - make a fake cell 
         if (in_fovea(x,y)) continue;         // but not if in fovea

         c = (Cell *)malloc(sizeof(Cell));
         c->p.x       = x;
         c->p.y       = y;
         //Point po = {ONH_X, ONH_Y};
         //double theta = atan2((double) y - (double)ONH_Y, (double) x - (double)ONH_X);
         //c->distToOnh = DIST(c->p,po) - ONH_EDGE(theta);
         c->count     = 0;
         int distFromFovea = MACULAR_DIST_SQ(c->p);
         c->thickness = MAX_AXON_COUNT(sqrt(distFromFovea));
         c->flag      = 0;
         c->alternate = NULL;
if (debug) {
if (target->path->next == NULL)
printf("Fake cell (%5d,%5d) -> (%5d,%5d)\n",x,y,-1, -1);
else
printf("Fake cell (%5d,%5d) -> (%5d,%5d) th=%u\n",x,y,target->path->next->c->p.x, target->path->next->c->p.y,c->thickness);
}

         Node *p = malloc(sizeof(Node));  // self then target.path->next
         p->c = c;
         p->next = target->path->next;
         c->path = p;

         grid[x][y].soma = c;
      } else {
         if (c->flag)         continue;    // on current path
         if (c->path == NULL) continue;    // has no path
         if (!IS_ROOM(c))     continue;    // no room
      }
      minC = c;
   }
//if (minC == NULL)
//printf(" gets NULL  NULL\n");
//else
//printf(" gets %5d %5d\n",minC->p.x, minC->p.y);

   target->flag = 0; // reset current flag

   return minC;
}//findNewPath()

/*
** For cell cc, whose nearest (completed) neighbour is c, begin a path at c->path.
** Follow c's path as far as possible, and jump to a new c if needed.
** Return 1 if success, 0 if fail to find path.
*/
int 
makeOnePath(int icc, Cell *target, Grid **grid) {
   Cell *current = cellBlock + icc;
/*
if (grid[9997][10445].soma != NULL) {
    Node *n = grid[9997][10445].soma->path;
    if (grid[9997][10445].soma->alternate != NULL)
        n = grid[9997][10445].soma->alternate->path;
    if (n != NULL) {
        printf("(%5d,%5d)", n->c->p.x, n->c->p.y);
        if (n->next != NULL) printf("->(%5d,%5d)", n->next->c->p.x, n->next->c->p.y);
        while (n->next != NULL)
            n = n->next;
        printf("...(%5d,%5d)\n", n->c->p.x, n->c->p.y);
    } else
        printf("( 9997,10445)...NULL?)\n");
}
*/

if (debug) printf("Current = (%5d,%5d)\n", current->p.x, current->p.y);
      // First, put in a start of path node at self
      // Note no space checking (assuming axon can start here)
   current->path = (Node *)malloc(sizeof(Node));
   current->path->c         = current;
   current->path->next      = NULL;
   INC_COUNT(current);

   Node *tail = current->path; // last entry in current path (new)
   current->flag = 1;
   int result = 0;

   while (target != NULL) {
if (debug)printf("\ttarget = (%5d,%5d)\n", target->p.x, target->p.y);
      result = 0; // assume fail
             // check if we can just follow same path (to save memory)
      Node *n = target->path;
if (debug)printf("\t\tcheck (%5d, %5d) count=%d < th=%u\n", n->c->p.x, n->c->p.y,n->c->count, n->c->thickness);
      while (n != NULL && IS_ROOM(n->c))
{
         n = n->next;
if (debug)if (n!=NULL) printf("\t\tcheck (%5d, %5d) count=%d < th=%u\n", n->c->p.x, n->c->p.y,n->c->count, n->c->thickness);
}
      if (n == NULL) { // hooray! just incrememnt count in each cell on path 
if (debug)printf("\tNo Copy required\n");
         tail->next = target->path;
         n = target->path;
         while ((n != NULL) && (n->next != NULL)) {
            INC_COUNT(n->c);
            n = n->next;
         }
         target = NULL;
         result = 1; // yay, succeed
      } else {
            // Path has to branch, so 
            // we need to take a copy of path up to the point where there's
            // no room, make a new node there and follow on
            // Mark all path nodes with flag = 1 to assist findNewPath
         Node *follow = target->path;
         while (IS_ROOM(follow->c)) {
            Node *n = (Node *)malloc(sizeof(Node)); assert(n != NULL);
            n->c    = follow->c;
            n->next = NULL;
            tail->next = n;
            tail = n;
            INC_COUNT(n->c);
            n->c->flag  = 1;
            follow = follow->next;
         }

           // Note alternate could end up being NULL
         if ((follow->c->alternate == NULL) || !IS_ROOM(follow->c->alternate))
            follow->c->alternate = findNewPath(tail->c, follow->c, grid);

         target = follow->c->alternate;
      }
   }

      // reset all the flags along path
   Node *n = current->path;
   while (n != NULL) {
      n->c->flag = 0;
      n = n->next;
   }
   return result;
}//makeOnePath()

/*
** This version uses precomputed scanPoints array, but it might not be 
** big enough when DENSE_SCALE is small
*/
Cell *
findClosestCompleted_restrictedArea(int i, Grid **grid) {
   Cell *target = cellBlock + i;
   Cell *minC = NULL;
   for(PointD *s = scanPoints ; (s < scanPoints + scanPointLen) && (minC == NULL) ; s++) {
      int x = target->p.x + s->p.x;
      int y = target->p.y + s->p.y;
      if ((x < 0) || (x >= SIZE)) continue;        // off grid
      if ((y < 0) || (y >= SIZE)) continue;
      Point p = {x,y};
      if (cross_raphe(target->p, p)) continue;    

      minC = grid[x][y].soma; // note could be NULL

      if (minC != NULL && minC->path == NULL)                  // no path yet
         minC = NULL;
      if (minC != NULL && (minC->count >= minC->thickness))    // no room
         minC = NULL;
   }
   return minC;
}//findClosestCompleted_restrictedArea()

/*
** Find the closest cell to cellBlock[i] in cellBlock[0..i-1]
** towards the ONH
** ASSUMES cellBlock is sorted in increasing distToOnh
*/
Cell *
findClosestCompleted(int i, Grid **grid) {
   Cell *res = findClosestCompleted_restrictedArea(i, grid);
   if (res != NULL)
      return res;

   float minD = DIST(cellBlock[i-1].p, cellBlock[i].p);
   int   minJ = i-1;
   for(int j = i-2 ; j >= 0 ; j--) {
      if (cellBlock[j].path  == NULL)                   continue;  // no path yet
      if (cellBlock[j].count >= cellBlock[j].thickness) continue;  // no room
      if (cross_raphe(cellBlock[i].p, cellBlock[j].p))  continue;  // crosses raphe
      float dist = DIST(cellBlock[j].p, cellBlock[i].p);
      if ((dist < minD)) {
         minD = dist;
         minJ = j;
      }
   }
   return cellBlock + minJ;
}//findClosestCompleted()

/*
** Print thickness values at each degree on a circle of radius radius mm

needs work as count now in cells, so have to find closest, etc.. Horrible!

void
print_oct(Grid **grid, float radius) {
   for(float theta = 0 ; theta < 360 ; theta++) {
      int x = (int)round(radius * cos(theta)*PIXELS_PER_MM + ONH_X);
      int y = (int)round(radius * sin(theta)*PIXELS_PER_MM + ONH_X);
      printf("%3.0f %d\n",theta, grid[x][y]);
   }
}//print_oct()
*/

/*
** For each cell in cellBlock (in order of increasing dist from ONH)
**   If within START_DIST, make a one node path
**   else find the closest completed and join paths with it
*/
void 
process(int size, Grid **grid) {
   int i = 0;
   float distToOnh = 0;
   for( ; i < numCells && distToOnh < START_DIST ; i++) {
      Point po = {ONH_X, ONH_Y};
      double theta = atan2((double) cellBlock[i].p.y - (double)ONH_Y, (double) cellBlock[i].p.x - (double)ONH_X);
      distToOnh = DIST(cellBlock[i].p,po) - ONH_EDGE(theta);

      if (distToOnh < START_DIST) {
            // just put self in path
         cellBlock[i].path             = (Node *)malloc(sizeof(Node));
         cellBlock[i].path->c          = cellBlock + i;
         cellBlock[i].path->next       = NULL;
         cellBlock[i].thickness        = 10000000; // UINT_MAX - 1;
         #ifdef PRINT_ENDPOINTS
         printf("# S ");
         print_path(cellBlock + i, FALSE);
         #endif
      }
   }

   for( ; i < numCells ; i++) {
      if (cellBlock[i].p.x > 14000) continue;
      //if (cellBlock[i].p.x > ONH_X) continue;
      //if (cellBlock[i].p.x < 5000) continue;
      //if (cellBlock[i].p.y < 5000) continue;
      //if (cellBlock[i].p.y > 15000) continue;
      //Cell *closest = cellBlock + findClosestCompleted(i);
      Cell *closest = findClosestCompleted(i, grid);
      if (closest == NULL) {
         printf("# Kn %d %d\n",cellBlock[i].p.x,cellBlock[i].p.y);
         continue;
      }
      if (makeOnePath(i, closest, grid)) { 
         #ifdef PRINT_ENDPOINTS
         print_path(cellBlock + i, FALSE);
         #endif
         #ifdef PRINT_PATHS
         if (i % PRINT_PATHS == 0)
            print_path(cellBlock + i, TRUE);
         #endif
      } else {
         printf("# K %d %d\n",cellBlock[i].p.x, cellBlock[i].p.y);
         grid[cellBlock[i].p.x][cellBlock[i].p.y].soma = NULL;
      }
   }
   #ifdef PRINT_OCT_PROFILE
   print_oct(grid, 3.4/2.0);
   #endif

}//process()

/*
** 
*/
int
main() {
#ifdef G_THREADS_ENABLED
   fprintf(stderr,"# threads enabled\n");
#else
   fprintf(stderr,"# threads not enabled\n");
#endif

#ifdef G_THREADS_IMPL_POSIX
   fprintf(stderr,"# threads posix\n");
#else
   fprintf(stderr,"# threads not posix\n");
#endif

   fprintf(stdout,"# SCAN_POINTS_RADIUS %10d\n",SCAN_POINTS_RADIUS);
   fprintf(stdout,"# DENSE_SCALE        %10.4f\n",DENSE_SCALE);
   fprintf(stdout,"# MAX_THICK          %10d\n",MAX_THICK);
   fprintf(stdout,"# MACULAR_RADIUS     %10.4f mm\n",(float)MACULAR_RADIUS/(float)PIXELS_PER_MM);
   fprintf(stdout,"# THETA_LIMIT        %10.4f degrees\n",THETA_LIMIT*180.0/M_PI);
   fprintf(stdout,"# ONH_X              %10d\n",ONH_X);
   fprintf(stdout,"# ONH_Y              %10d\n",ONH_Y);
   fprintf(stdout,"# MAJOR AXIS         %10d\n",ONH_MAJOR);
   fprintf(stdout,"# MINOR AXIS         %10d\n",ONH_MINOR);

   g_thread_init(NULL);
   gdk_threads_init();     /* Secure gtk */
   gdk_threads_enter();    /* Obtain gtk's global lock */

   int size;
   Grid **grid;
   init_grid(&size, &grid);

   init_cells();

   init_scanPoints();

//for(int i = 0 ; i < numCells ; i++)
//if (MACULAR_DIST(cellBlock[i].p) < MACULAR_RADIUS)
//printf("im %d\n",i);
//return 0;
   fprintf(stdout, "# Number of cells: %d\n",numCells);
   process(size, grid);

   /* Release gtk's global lock */
   gdk_threads_leave();

   return 0;
}
