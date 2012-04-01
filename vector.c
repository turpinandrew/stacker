#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

Vector *vector_new(int nelem, int size) {
    if (nelem < 1) nelem = 1;

    Vector *v = (Vector *)malloc(sizeof(Vector));
    assert(v != NULL);
    v->p     = malloc(size * nelem);
    assert(v->p != NULL);
    v->n     = 0;
    v->size  = size;
    v->max_n = nelem;
    return v;
}//vector_new()
    
void vector_empty(Vector *v) 
{ 
   if (v != NULL)
      v->n = 0; 
}

void vector_free(Vector *v) 
{
   if (v != NULL) {
      free(v->p);
       v->max_n = 0;
       v->n = 0;
   }
}

   // add object *p to vector v. 
void vector_add(Vector *v, void *p) 
{
   if (v != NULL) {
      if (v->n == v->max_n) {
          v->p = realloc(v->p, v->size * v->max_n * 2);
          v->max_n *= 2;
      }
      memcpy(vector_item_at(v, v->n), p, v->size);
      v->n += 1;
   }
}

   // return pointer to i'th item in v->p
void *vector_item_at(Vector *v, int i) {
   if (v == NULL)
      return NULL;

   return (void *)(  ((unsigned char *)(v->p)) + ((i) * (v)->size) ) ;
}

   // Remove element at i by simply swapping the final entry to position i, 
   // and reducing n
void vector_delete_at(Vector *v, int i) {
   if (v == NULL || i >= v->n)
      return;

   memcpy(vector_item_at(v, i), vector_item_at(v, v->n-1), v->size);
   v->n -= 1;
}
   
