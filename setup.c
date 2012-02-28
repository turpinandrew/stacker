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
         grid[x][y].soma = NULL;
      }
   }
   return NULL;
}//init_grid_piece()

/* 
   Set square grid SIZE*SIZE all to NULL, *** except in fovea, raphe and ONH where set to 1 
   Assumes fovea is at (SIZE/2, SIZE/2)
   Assumes raphe is at (0...SIZE/2, SIZE/2)
   Sets *size to SIZE
*/
void init_grid(int *size, Grid ***inGrid) 
{
   *size = SIZE;

   grid = (Grid **) malloc(sizeof(Grid *)*SIZE);
   assert(grid != NULL);
   for(int x = 0 ; x < SIZE ; x++) {
      grid[x] = (Grid *) malloc(sizeof(Grid)*SIZE);
      assert(grid[x] != NULL);
   }

   fprintf(stderr,"Initialising grid\n");

        //Create threads into an array threads[0..THREADS-1]
   BB *bb = (BB *)malloc(sizeof(BB) * (THREADS+1));
   for(int i = 0 ; i < THREADS+1 ; i++) {
      bb[i].tlx = (int)round((float)    i   * (float)SIZE / ((float)THREADS+1.0))    ;
      bb[i].brx = (int)round(((float)i+1.0) * (float)SIZE / ((float)THREADS+1.0)) - 1;
      bb[i].tly = 0;
      bb[i].bry = SIZE - 1;
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

      // block out fovea
   for(int x = -FOVEA_RADIUS ; x <= +FOVEA_RADIUS ; x++)
      for(int y = -FOVEA_RADIUS ; y <= +FOVEA_RADIUS ; y++)
         if (x*x + y*y <= FOVEA_RADIUS * FOVEA_RADIUS)
            grid[SIZE/2+x][SIZE/2+y].soma = (Cell *)1;

      // add horizontal raphe
   for(int x = 0 ; x < SIZE/2 ; x++) {
      int y = SIZE/2;      // XXX might want to change this to be dependant on Fov->ONH angle
      grid[x][y].soma = (Cell *)1;
   }

      // block out ONH (+-2 in loops to catch rounding errors)
   for(int x = ONH_X - ONH_MAJOR -2 ; x <= ONH_X + ONH_MAJOR + 2; x++)
      for(int y = ONH_Y - ONH_MINOR -2 ; y <= ONH_Y + ONH_MINOR +2 ; y++) {
         double theta = atan2((double) y - (double)ONH_Y, (double) x - (double)ONH_X);
         Point p = {x,y};
         Point po = {ONH_X, ONH_Y};
         if (DIST(p,po) <= ONH_EDGE(theta) + 1)  // +1 just to get at least one pixel away
            grid[x][y].soma = (Cell *)1;
      }

   *inGrid = grid;

/*
   for(int x = 0 ; x < SIZE ; x++)
      for(int y = 0 ; y < SIZE ; y++)
         if (grid[x][y] > 0)
            printf("%d %d %d\n",x,y,grid[x][y]);
*/

}//init_grid()

// *** Note reset soma to NULL if soma == 1
static gpointer make_cell_piece(gpointer data) {
   BB *bb = (BB *)data;
//printf("make Cell gp: (%d,%d) -> (%d,%d)\n",bb->tlx,bb->tly,bb->brx, bb->bry);fflush(stdout);
   for(int x = bb->tlx ; x <= bb->brx ; x++) {
      for(int y = bb->tly ; y <= bb->bry ; y++) {
            // check room, not in fovea, not in ONH, not on raphe
         if (grid[x][y].soma) {
            grid[x][y].soma = NULL;      // reset soma
            continue;
         }
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
      gsl_rng_set(bb[i].rng, time(NULL) * (i+1));
      //gsl_rng_set(bb[i].rng, 1 * (i+1));

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

   fprintf(stderr,"Number of cells = %d\n",numCells);

   qsort(loc, numCells, sizeof(PointD), cmp_PointD);

   cellBlock = (Cell *)malloc(sizeof(Cell) * numCells);
   assert(cellBlock != NULL);
   for (int i = 0 ; i < numCells ; i++) {
      cellBlock[i].p         = loc[i].p;
      cellBlock[i].path      = NULL;
      //cellBlock[i].distToOnh = loc[i].dist;
      cellBlock[i].count     = 0;
      cellBlock[i].flag      = 0;
      int distFromFovea = MACULAR_DIST_SQ(loc[i].p);
      if (distFromFovea > MACULAR_RADIUS_SQ)
         cellBlock[i].thickness = UCHAR_MAX;
      else
         cellBlock[i].thickness = MAX_AXON_COUNT(sqrt(distFromFovea));
      cellBlock[i].alternate = NULL;
   
      grid[loc[i].p.x][loc[i].p.y].soma = cellBlock + i;
   }
   return;
}//init_cells()
