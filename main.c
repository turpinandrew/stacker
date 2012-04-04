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
#include <stdlib.h>
#include "setup.h"
#include "queue.h"
#include "density.h"
#include "types.h"
#include "vector.h"

#define min(_a, _b) ((_a) < (_b) ? (_a) : (_b))
#define max(_a, _b) ((_a) > (_b) ? (_a) : (_b))

int debug = 0; 

#define PRINT_ENDPOINTS
#define PRINT_PATHS 500
//#define PRINT_OCT_PROFILE

#define IS_ROOM(_ge) ((_ge)->count < (_ge)->thickness)

   //An LUT for efficiency.
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
   scanPoints = (PointD *) malloc(sizeof(PointD)*(2*NEW_PATH_RADIUS_LIMIT+1)*(2*NEW_PATH_RADIUS_LIMIT+1));
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
** Convert cell location into grid location
*/
Grid *
toGridCoords(Point p, Grid **grid) {
   int x = (int)floor(p.x / (float)SCALE);
   int y = (int)floor(p.y / (float)SCALE);
   return &grid[x][y];
}//toGridCoords()

   
/*
** Return 1 if a and b are on different sides of the raphe, 0 otherwise
** NOTE a and b are in cell coords (0..SIZE) not grid coords
*/
int
cross_raphe(Point a, Point b) {
    if (a.x > SIZE/2 || b.x > SIZE/2) return 0;  // there is no spoon (raphe)
    if (a.y < SIZE/2) 
        return (b.y > SIZE/2);

    return (b.y < SIZE/2);
}//cross_raphe()


/*
** Find closest grid element in the direction of target from current (+-THETA_LIMIT)
** that is not target, that has room, and is not already on current path (flag == 1).
**
** Use a spiral search (well, they're squares of increasing dist) outwards from
** current ignoring cells with flag==1, points outside theta range, and 
** count >= thickness.
**
** If aboveRaphe is 1, can only find a path with y-coord <= GRID_SIZE/2
** If aboveRaphe is 0, can only find a path with y-coord >= GRID_SIZE/2
*/
Grid *
findNewPath(Grid *current, Grid *target, Grid **grid, int aboveRaphe, Grid *onh) {
   if (debug)printf("\n# current = %5d %5d ",current->p.x, current->p.y);
   if (debug)printf(" wants %5d %5d",target->p.x, target->p.y);
   if (debug)printf(" with aboveRaphe=%d\n",aboveRaphe);

   double theta = atan2(target->p.y - current->p.y, target->p.x - current->p.x);
   double hiTheta = theta + THETA_LIMIT;
   double loTheta = theta - THETA_LIMIT;
   if (debug)printf("\ttheta range (%10.6f, %10.6f)\n",loTheta, hiTheta);

   target->flag = 1; // rule out current
   Grid *g, *minG = NULL;
   for(PointD *s = scanPoints ; (s < scanPoints + scanPointLen) && (minG == NULL) ; s++) {
      if (s->dist < loTheta) continue;    // outside theta range
      if (s->dist > hiTheta) continue;
      int x = target->p.x + s->p.x;
      int y = target->p.y + s->p.y;
      if ((x < 0) || (x >= GRID_SIZE)) continue;        // off grid
      if ((y < 0) || (y >= GRID_SIZE)) continue;
      if (aboveRaphe  && (current->p.x < onh->p.x) && (y > GRID_SIZE/2)) continue;
      if (!aboveRaphe && (current->p.x < onh->p.x) && (y < GRID_SIZE/2)) continue;

      g = &grid[x][y];

      if (debug)printf("   check (%4d,%4d) thick=%d flag=%d\n",g->p.x, g->p.y, g->thickness,g->flag);

      if (g->flag)         continue;    // on current path
      if (!IS_ROOM(g))     continue;    // no room

      minG = g;
   }

   if (debug && minG)printf("\tChosen (%10d, %10d)\n",minG->p.x, minG->p.y);
   if (debug && !minG)printf("\tChosen (NULL, NULL)\n");
   target->flag = 0; // reset current flag

   if (minG == NULL)
        return NULL;
   if (minG->lastPathThrough == NULL)  // No axons have been through here before
      minG->lastPathThrough = target->lastPathThrough;   

   return minG;
}//findNewPath()

void unsetFlag(void *v) { 
    Grid *g = *((Grid **)v);
    g->flag = 0; 
    printf("    Unset (%d,%d) %lx\n",g->p.x,g->p.y,g);
}
void printGrid(void *v) { 
   Grid *g = *((Grid **)v); 
   printf("%6d %6d\n",g->p.x, g->p.y);
}

/*
** Make a path from grid element current towards onh, with first step being
** through grid element target.
** Note: current must not equal target.
**
** Print path if printIt == 1
*/
void 
makeOnePath(Grid *current, Grid *target, Grid **grid, char printIt) {
   assert(current != target);

   int aboveRaphe = rand() % 2;
   if (current->p.y < GRID_SIZE/2) aboveRaphe = 1;
   if (current->p.y > GRID_SIZE/2) aboveRaphe = 0;

   if (debug) printf("Current = (%5d,%5d) aboveRaphe=%d\n", current->p.x, current->p.y, aboveRaphe);
   if (debug) printf("Target  = (%5d,%5d)\n", target->p.x, target->p.y);

   Vector *gridUsed = vector_new(100, sizeof(Grid *));
   current->flag = 1;
   current->count += 1;
   vector_add(gridUsed, &current);

   Point onhP   = {ONH_X, ONH_Y};
   Grid *onh    = toGridCoords(onhP, grid);

   if (target == onh)                   // special case initial target
      vector_add(gridUsed, &target);     // is onh

   while ((target != onh) && (target != NULL)) {
      if (debug)printf("\ttarget = (%5d,%5d)\n", target->p.x, target->p.y);

      if (!IS_ROOM(target)) 
         target = findNewPath(current, target, grid, aboveRaphe, onh);

      if (target != NULL) {
         current->lastPathThrough = target;
         target->count += 1;
         target->flag = 1;
         vector_add(gridUsed, &target);
         current = target;
         target = target->lastPathThrough;
      } else {
         current->lastPathThrough = NULL; // don't come to this target again
      }
//if (vector_length(gridUsed) > 100) target = NULL;   // saftey net against the infinite
   }

   //printf("Grid\n");
   //vector_apply(gridUsed, printGrid); // reset all the flags along path
   //printf("End Grid\n");
   vector_apply(gridUsed, unsetFlag); // reset all the flags along path

   //if (printIt && vector_length(gridUsed) > 3) {
   if (printIt && target == onh) {
      vector_apply(gridUsed, printGrid);
      printf("%6d %6d\n",onh->p.x,onh->p.y);
      printf("%6d %6d\n",-1,-1);
   }

   #ifdef PRINT_ENDPOINTS
   Grid *first = *((Grid **)vector_item_first(gridUsed));
   if (target == onh) {
      Grid *last  = *((Grid **)vector_item_last(gridUsed));
      printf("%6d %6d ",first->p.x, first->p.y);
      printf("%6d %6d ",last->p.x, last->p.y);
      printf("%6d\n",vector_length(gridUsed));
   } else {
      printf("# K %d %d\n",first->p.x, first->p.y);
   }
   #endif
   vector_free(gridUsed);

   return;
}//makeOnePath()

/*
** Find the closest cell to cellBlock[i] in cellBlock[0..i-1]
** towards the ONH. If there is already a path through this grid cell
** then use it.
** ASSUMES cellBlock is sorted in increasing distToOnh
** ASSUMES cells in [0..i-1] have all completed
*/
Cell *
findClosestCompleted(int i, Grid **grid) {
   float minD = DIST(cellBlock[i-1].p, cellBlock[i].p);
   int   minJ = i-1;
   for(int j = i-2 ; j >= 0 ; j--) {
      if (cross_raphe(cellBlock[i].p, cellBlock[j].p))  continue;  // crosses raphe
      Grid *g = toGridCoords(cellBlock[j].p, grid);
      if (!IS_ROOM(g)) continue;  // no room
      if (g->lastPathThrough == NULL) continue;       // just in case some of cellBlock has been skipped
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
process(Grid **grid) {
   int i = 0;
   float distToOnh = 0;
   for( ; i < numCells && distToOnh < START_DIST ; i++) {
      Point po = {ONH_X, ONH_Y};
      double theta = atan2((double) cellBlock[i].p.y - (double)ONH_Y, (double) cellBlock[i].p.x - (double)ONH_X);
      distToOnh = DIST(cellBlock[i].p,po) - ONH_EDGE(theta);

      if (distToOnh < START_DIST) {
            // just a one step path to onh
         Grid *self = toGridCoords(cellBlock[i].p, grid);
         Grid *onh  = toGridCoords(po, grid);
         self->lastPathThrough = onh;
         self->count += 1;
         #ifdef PRINT_ENDPOINTS
         printf("# S ");
         printf("%6d %6d ",self->p.x, self->p.y);
         printf("%6d %6d\n",onh->p.x, onh->p.y);
         #endif
      }
   }

   for( ; i < numCells ; i++) {
      //if (cellBlock[i].p.x > ONH_X) continue;
      //if (cellBlock[i].p.x < 5000) continue;
      //if (cellBlock[i].p.y < 5000) continue;
      //if (cellBlock[i].p.y > 15000) continue;
      //Cell *closest = cellBlock + findClosestCompleted(i);
      Grid *curr = toGridCoords(cellBlock[i].p, grid);
   if (debug)printf("curr cell = %5d %5d ",cellBlock[i].p.x, cellBlock[i].p.y);
      Grid *targ = NULL;
      if (curr->lastPathThrough == NULL) {
         Cell *closest = findClosestCompleted(i, grid);
         if (closest == NULL)
            printf("# Kn %d %d\n",cellBlock[i].p.x,cellBlock[i].p.y);
         else {
            targ = toGridCoords(closest->p, grid);
            while (targ == curr) {  // try and find a cell outside current grid element
               closest = findClosestCompleted(closest - cellBlock, grid);
               targ = toGridCoords(closest->p, grid);
            }
         }
      } else
         targ = curr->lastPathThrough;

      if (targ != NULL)
         makeOnePath(curr, targ, grid, i % PRINT_PATHS == 0);
   }
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

   fprintf(stdout,"# DENSE_SCALE           %10.4f\n",DENSE_SCALE);
   fprintf(stdout,"# FOVEA_RADIUS          %10.4f mm %10.0f pixels\n",(float)FOVEA_RADIUS  /(float)PIXELS_PER_MM, FOVEA_RADIUS);
   fprintf(stdout,"# THETA_LIMIT           %10.4f degrees\n",THETA_LIMIT*180.0/M_PI);
   fprintf(stdout,"# ONH_X                 %10d\n",ONH_X);
   fprintf(stdout,"# ONH_Y                 %10d\n",ONH_Y);
   fprintf(stdout,"# MAJOR AXIS            %10d\n",ONH_MAJOR);
   fprintf(stdout,"# MINOR AXIS            %10d\n",ONH_MINOR);
   fprintf(stdout,"# SCALE                 %10.4f\n",SCALE);
   fprintf(stdout,"# GRID_SIZE             %10d\n",GRID_SIZE);
   fprintf(stdout,"# NEW_PATH_RADIUS_LIMIT %10d\n",NEW_PATH_RADIUS_LIMIT);

   g_thread_init(NULL);
   gdk_threads_init();     /* Secure gtk */
   gdk_threads_enter();    /* Obtain gtk's global lock */

   Grid **grid;
   init_grid(&grid);

   Point onhP   = {ONH_X, ONH_Y};
   Grid *onh    = toGridCoords(onhP, grid);
   onh->thickness = 10000000;

   init_cells();

   init_scanPoints();

//for(int i = 0 ; i < numCells ; i++)
//if (MACULAR_DIST(cellBlock[i].p) < MACULAR_RADIUS)
//printf("im %d\n",i);
//return 0;
   srand(time(NULL));
   fprintf(stdout, "# Number of cells: %d\n",numCells);
   process(grid);

   /* Release gtk's global lock */
   gdk_threads_leave();

   return 0;
}
