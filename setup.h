#include "types.h"
#include "queue.h"

#define PIXELS_PER_MM 1000

#define DENSE_SCALE 1.2 // multiply density of cells by this factor

#define SIZE (20 * PIXELS_PER_MM) //  20 mm square for cells

#define AXIAL_LENGTH 25.0    // mm

#define FOVEA_RADIUS ( 0.2 * (float)PIXELS_PER_MM )

#define ONH_X (int)((SIZE/2 + (15.0 /180.0*M_PI*AXIAL_LENGTH/2.0*(float)PIXELS_PER_MM))) // degrees -> grid
#define ONH_Y (int)((SIZE/2 + ( 2.0 /180.0*M_PI*AXIAL_LENGTH/2.0*(float)PIXELS_PER_MM))) // degrees -> grid
#define ONH_MAJOR (int)(1.66/2.0 * (double)PIXELS_PER_MM) /* 1.66mm 2*major axis (x) of optic nerve */
#define ONH_MINOR (int)(1.94/2.0 * (double)PIXELS_PER_MM) /* 1.94mm 2*minor axis (y) of optic nerve */

#define DIST(_p1, _p2)  (float)sqrt( ((double)(_p1).x-(double)(_p2).x)*((double)(_p1).x-(double)(_p2).x) + ((double)(_p1).y-(double)(_p2).y)*((double)(_p1).y-(double)(_p2).y) )

#define ONH_EDGE(_theta) ((float)ONH_MAJOR*(float)ONH_MINOR/sqrt( \
                          (float)ONH_MINOR*cos(_theta)*(float)ONH_MINOR*cos(_theta) + \
                          (float)ONH_MAJOR*sin(_theta)*(float)ONH_MAJOR*sin(_theta)))

#define START_DIST 0.2*PIXELS_PER_MM    // number of pixels from ONH_EDGE to be included in first pool

   // macros to control thickness
#define MACULAR_RADIUS (3 * PIXELS_PER_MM)
#define MACULAR_RADIUS_SQ (MACULAR_RADIUS * MACULAR_RADIUS)
#define MACULAR_DIST(_p) (int)round(sqrt( ((_p).x-SIZE/2)*((_p).x-SIZE/2)+((_p).y-SIZE/2)*((_p).y-SIZE/2)))
#define MACULAR_DIST_SQ(_p) ( ((_p).x - SIZE/2)*((_p).x - SIZE/2) + ((_p).y - SIZE/2)*((_p).y - SIZE/2))

#define min(_a, _b) ((_a) < (_b) ? (_a) : (_b))
#define max(_a, _b) ((_a) > (_b) ? (_a) : (_b))
//#define MAX_THICK 20
//#define MAX_AXON_COUNT(_dist) (int)round(((float)(_dist)-(float)FOVEA_RADIUS)/(float)MACULAR_RADIUS * (float)MAX_THICK*(float)DENSE_SCALE)
   // linear from (FOVEA_RADIUS,0) to (MACULAR_RADIUS, MAX_THICK)
#define MAX_THICK 30
#define MAX_AXON_COUNT(_dist) (min(MAX_THICK, (int)round(((float)(_dist)-(float)FOVEA_RADIUS)/(float)MACULAR_RADIUS * (float)MAX_THICK)) * SCALE * SCALE)

#define SCALE 200.0   // sqrt of how many cell sqaures per Grid square (so grid is SIZE/SCALE * SIZE/SCALE)

#define GRID_SIZE ((int)ceil((float)SIZE/(float)SCALE))

    // how far to search for a new path during growth?
#define NEW_PATH_RADIUS_LIMIT ((int)round(0.2*(double)SCALE))

    // don't search outside +- this from proposed trajectory during growth
#define THETA_LIMIT  (M_PI/3.0)     

   // all in pixels
void init_grid(Grid ***grid);
void init_cells();
int cmp_PointD(const void *a, const void *b);
