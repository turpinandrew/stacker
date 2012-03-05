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

typedef struct cell Cell; 
typedef struct node Node;
struct node {
   Cell *c;
   Node *next;
};

struct cell {
   Point p;
   Node *path;       // a linked list of cells that path jumps along
                     // first element of path is self
   //float distToOnh;
   unsigned int count;        // count of paths that go through this cell
   unsigned int thickness;    // number of paths that are allowed to go through this cell

   char flag;        // used for various, initially all 0 

   Cell *alternate;  // if a path tried to come through this, but was full so had to do find,
                     // this records the result of the find.
};

typedef struct grid {
   Cell *soma;  // ptr to cell at this location
} Grid; 

#endif
