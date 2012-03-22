#ifndef _TYPES_H_ 
#define _TYPES_H_ 

#include "queue.h"

typedef unsigned char uchar;

typedef struct point {
   short x;
   short y;
} Point;

typedef struct pointD {
   Point p;
   float dist;
} PointD;

typedef struct grid Grid;
struct grid {
   Point p;                   // location of self. Could infer it, but will store for now.
   unsigned int count;        // count of paths that go through this cell
   unsigned int thickness;    // number of paths that are allowed to go through this cell

   char flag;        // used for various, initially all 0 

   Grid *lastPathThrough;  // record of where the last path went that came through here.
}; 

typedef struct cell {
   Point p;
} Cell;

#endif
