/*
** Initialise a few key structures
**    - grid : list of cells with axons going through [x][y]
**    - cellBlock:  array of numCells cells
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <glib.h>
#include "types.h"
#include "setup.h"
#include "density.h"
#include "queue.h"

Grid **grid;      // grid[x][y] of axons passing through 
Cell *cellBlock;  // real cells [0..numCells-1]
int numCells;     // length of cellBlock

typedef struct bb { 
   int tlx,tly,brx,bry; // bounding box top-left and bottom-right
   gsl_rng *rng;        // random numbers
   int *numCells;       // num cells created in this bb
   PointD *loc;         // array of numCells locations of cells
} BB;  

static gpointer init_grid_piece(gpointer data) {
   BB *bb = (BB *)data;
//printf("Init gp: (%d,%d) -> (%d,%d)\n",bb->tlx,bb->tly,bb->brx, bb->bry);
   for(int x = bb->tlx ; x <= bb->brx ; x++) {
      for(int y = bb->tly ; y <= bb->bry ; y++) {
         Point p_xy = {x*SCALE,y*SCALE};
         int distFromFovea = MACULAR_DIST_SQ(p_xy);
         if (distFromFovea <= FOVEA_RADIUS * FOVEA_RADIUS)
            grid[x][y].thickness = 0;
         else
            grid[x][y].thickness = MAX_AXON_COUNT(sqrt(distFromFovea)) * SCALE;
         grid[x][y].count           = 0;
         grid[x][y].flag            = 0;
         grid[x][y].lastPathThrough = NULL;
         grid[x][y].p.x             = x;
         grid[x][y].p.y             = y;
      }
   }
   return NULL;
}//init_grid_piece()

/* 
   Set thickness for square grid GRID_SIZE
*/
void init_grid(Grid ***inGrid) 
{
   grid = (Grid **) malloc(sizeof(Grid *)*GRID_SIZE);
   assert(grid != NULL);
   for(int x = 0 ; x < GRID_SIZE ; x++) {
      grid[x] = (Grid *) malloc(sizeof(Grid)*GRID_SIZE);
      assert(grid[x] != NULL);
   }

   fprintf(stderr,"Initialising grid\n");

        //Create threads into an array threads[0..THREADS-1]
   BB *bb = (BB *)malloc(sizeof(BB) * (THREADS+1));
   for(int i = 0 ; i < THREADS+1 ; i++) {
      bb[i].tlx = (int)round((float)    i   * (float)(GRID_SIZE) / ((float)THREADS+1.0))    ;
      bb[i].brx = (int)round(((float)i+1.0) * (float)(GRID_SIZE) / ((float)THREADS+1.0)) - 1;
      bb[i].tly = 0;
      bb[i].bry = GRID_SIZE - 1;
      bb[i].rng = NULL; 
   }
   GError    *error = NULL;
   GThread **threads = (GThread **)malloc(sizeof(GThread *) * (THREADS));
   for(int i = 0 ; i < THREADS ; i++)
      threads[i] = g_thread_create( init_grid_piece, (gpointer)(bb + i) , TRUE, &error );
   init_grid_piece((gpointer) (bb + THREADS));

   for(int i = 0 ; i < THREADS ; i++)
      g_thread_join(threads[i]);
   free(bb);
   free(threads);

/*
      // block out fovea
   int fRad = round((float)FOVEA_RADIUS / (float)SCALE);
   for(int x = -fRad ; x <= +fRad ; x++)
      for(int y = -fRad ; y <= +fRad ; y++)
         if (x*x + y*y <= fRad * fRad)
            grid[GRID_SIZE/2+x][GRID_SIZE/2+y].thickness = 0;

      // add horizontal raphe
   for(int x = 0 ; x < GRID_SIZE/2 ; x++) {
      int y = GRID_SIZE/2;  // XXX might want to change this to be dependant on Fov->ONH angle
      grid[x][y].thickness = 0;
   }

      // block out ONH (+-2 in loops to catch rounding errors)
   for(int x = ONH_X - ONH_MAJOR -2 ; x <= ONH_X + ONH_MAJOR + 2; x++)
      for(int y = ONH_Y - ONH_MINOR -2 ; y <= ONH_Y + ONH_MINOR +2 ; y++) {
         double theta = atan2((double) y - (double)ONH_Y, (double) x - (double)ONH_X);
         Point p = {x,y};
         Point po = {ONH_X, ONH_Y};
         if (DIST(p,po) <= ONH_EDGE(theta) + 1)  // +1 just to get at least one pixel away
            grid[(int)floor((float)x/(float)SCALE)][(int)floor((float)y/(float)SCALE)].thickness = 0;
      }
*/

   *inGrid = grid;
}//init_grid()

static gpointer make_cell_piece(gpointer data) {
   BB *bb = (BB *)data;
//printf("make Cell gp: (%d,%d) -> (%d,%d)\n",bb->tlx,bb->tly,bb->brx, bb->bry);fflush(stdout);
   Point po = {ONH_X, ONH_Y};
   for(int x = bb->tlx ; x <= bb->brx ; x++) {
      for(int y = bb->tly ; y <= bb->bry ; y++) {
            // check room, not in fovea, not in ONH, not on raphe
         if ((x-SIZE/2)*(x-SIZE/2) + (y-SIZE/2)*(y-SIZE/2) <= FOVEA_RADIUS * FOVEA_RADIUS)
            continue;
         if ((x < SIZE/2) && (fabs(y-SIZE/2) < FOVEA_RADIUS/2))   // assume raphe as thick as fovea/2
            continue;
         Point p = {x,y};
         double theta = atan2((double) y - (double)ONH_Y, (double) x - (double)ONH_X);
         if (DIST(p,po) <= ONH_EDGE(theta) + 1)  // +1 just to get at least one pixel away
            continue;

            // flip coin...
         double prob = find_density(((float)x-(float)SIZE/2.0)/(float)PIXELS_PER_MM, ((float)y-(float)SIZE/2.0)/(float)PIXELS_PER_MM) / (double)PIXELS_PER_MM / (double)PIXELS_PER_MM*DENSE_SCALE;
         if (gsl_rng_uniform(bb->rng) < prob) {
//printf("\t(%d,%d) %d\n",bb->tlx,bb->tly,*(bb->numCells));fflush(stdout);
            double theta = atan2((double) y - (double)ONH_Y, (double) x - (double)ONH_X);
            Point p = {x,y};
            Point po = {ONH_X, ONH_Y};
            bb->loc[*(bb->numCells)].p = p;
            bb->loc[*(bb->numCells)].dist = DIST(p,po) - ONH_EDGE(theta);
            *(bb->numCells) += 1;
         }
      }
   }
   return NULL;
}//make_cell_piece()

// sort by increasing dist
int cmp_PointD(const void *a, const void *b)
{
   PointD *aa = (PointD *)a;
   PointD *bb = (PointD *)b;
   if (aa->dist < bb->dist)
      return -1;
   if (aa->dist > bb->dist)
      return +1;
   return 0;
}

/*
** Allocate all cell memory, initialise cells, link them to grid[][].soma
**  - all cells are in outCellBlock[0..numCells-1]
**    sorted by increasing distToOnh
*/
void init_cells()
{
   fprintf(stderr,"\nMaking cells\n");
        //Create threads into an array threads[0..THREADS-1]
   BB *bb = (BB *)malloc(sizeof(BB) * (THREADS+1));
   for(int i = 0 ; i < THREADS+1 ; i++) {
      bb[i].tlx = (int)round((float)    i   * (float)SIZE / ((float)THREADS+1.0))    ;
      bb[i].brx = (int)round(((float)i+1.0) * (float)SIZE / ((float)THREADS+1.0)) - 1;
      bb[i].tly = 0;
      bb[i].bry = SIZE - 1;

      bb[i].rng = gsl_rng_alloc(gsl_rng_default);
#ifdef RANDOM_SEED_IS_TIME
      gsl_rng_set(bb[i].rng, time(NULL) * (i+1));
#else
      gsl_rng_set(bb[i].rng, 1 * (i+1));
#endif

      bb[i].numCells = (int *)malloc(sizeof(int));
      *(bb[i].numCells) = 0;
      bb[i].loc      = (PointD *)malloc(sizeof(PointD)*(bb[i].brx-bb[i].tlx+1)*(bb[i].bry-bb[i].tly+1));
   }
   GError    *error = NULL;
   GThread **threads = (GThread **)malloc(sizeof(GThread *) * (THREADS));
   for(int i = 0 ; i < THREADS ; i++)
      threads[i] = g_thread_create( make_cell_piece, (gpointer)(bb + i) , TRUE, &error );
   make_cell_piece((gpointer) (bb + THREADS));

      // for each pixel, set with prob = #rgs-in-this-square/PIXELS_PER_MM^2
   PointD *loc = (PointD *)malloc(sizeof(PointD) * SIZE * SIZE);  // a temporary list of locations
   numCells = 0;
   for(int i = 0 ; i < THREADS+1 ; i++) {
      if (i < THREADS)
         g_thread_join(threads[i]);

      for(int j = numCells ; j < numCells + *(bb[i].numCells) ; j++)
         loc[j] = bb[i].loc[j - numCells];

      numCells += *(bb[i].numCells);

      free(bb[i].numCells);
      free(bb[i].loc);
   }

   fprintf(stderr,"# Number of cells = %d\n",numCells);

   qsort(loc, numCells, sizeof(PointD), cmp_PointD);

   cellBlock = (Cell *)malloc(sizeof(Cell) * numCells);
   assert(cellBlock != NULL);
   for (int i = 0 ; i < numCells ; i++) 
      cellBlock[i].p = loc[i].p;

    free(loc);
   return;
}//init_cells()
