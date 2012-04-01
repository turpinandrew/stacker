#ifndef VECTOR_H
#define VECTOR_H

#include <string.h>
#include <assert.h>

typedef struct v
{
   int max_n;   // max number of items before realloc
   int n;       // number of items in vector
   int size;    // size of one element in p[]
   void *p;     // array of objects of size size.
} Vector;

#define vector_item_last(_v)   vector_item_at((_v), (_v)->n-1)       // return pointer to last item in p
#define vector_item_first(_v)  vector_item_at((_v), 0)               // return pointer to first item in p
#define vector_is_empty(_v)    ( ((_v)==NULL) ? 1 : (_v)->n == 0 )
#define vector_length(_v)      ( ((_v)==NULL) ? 0 : (_v)->n )

    // apply function _f(x) to each &x in (Vector *)_v
#define vector_apply(_v, _f) do {         \
   if ((_v) != NULL)                      \
      for(int _i = 0 ; _i < (_v)->n ; _i++)  \
          (_f)(vector_item_at((_v), _i));  \
} while(0)

    // apply function _f(x,_a) to each &x in _v
#define vector_apply2(_v, _f, _a) do {       \
   if ((_v) != NULL)                         \
      for(int _i = 0 ; _i < (_v)->n ; _i++)     \
         (_f)(vector_item_at((_v), _i), _a);  \
} while(0)

extern Vector *vector_new(int nelem, int size);
extern void    vector_empty(Vector *v);            // set n to 0, but do not free
extern void    vector_free(Vector *v);             // chuck all memory
extern void    vector_add(Vector *v, void *p);
extern void   *vector_item_at(Vector *v, int i);
extern void    vector_delete_at(Vector *v, int i);

#endif
