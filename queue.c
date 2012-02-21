#include <glib.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include "queue.h"

G_LOCK_DEFINE(queue_lock);

static Queue_node *new_node(void);

static Queue_node *new_node(void)
{
   Queue_node *n;
   n = (Queue_node *)malloc(sizeof(Queue_node));

   assert (n != NULL);

   return n;
}

Queue *new_empty_queue(void)
{
   Queue *q;
   q = (Queue *)malloc(sizeof(Queue));

   assert(q != NULL);

   q->head = NULL;
   q->tail = NULL;
   q->count = 0;
   return q;
}

void *remove_first(Queue *q)
{
   Queue_node *tmp;
   void *result;

   assert(q != NULL);

   #ifdef ICC
   #pragma warning (disable:1293) // icc complains about a may_alias attribute, which I think is gcc specific
   #pragma warning (disable:1292)
   #endif
   G_LOCK(queue_lock);

   /* nothing in the queue */
   if (q->head == NULL)
   {
      assert (q->tail == NULL);
      result = NULL;
   }
   /* one thing in the queue */
   else if (q->head == q->tail)
   {
      q->count--;
      tmp = q->head;
      result = tmp->cell;
      q->head = q->tail = NULL;
      free(tmp);
   }
   /* more than on thing in the queue */
   else
   {
      q->count--;
      tmp = q->head;
      result = tmp->cell;
      q->head = q->head->next;
      free(tmp);
   }
   #ifdef ICC
   #pragma warning (disable:1293) // icc complains about a may_alias attribute, which I think is gcc specific
   #pragma warning (disable:1292)
   #endif
   G_UNLOCK(queue_lock);

   return (result);
}

void insert_last(Queue *q, void *cell)
{
   Queue_node *node;

   assert (q != NULL);

   #ifdef ICC
   #pragma warning (disable:1293) // icc complains about a may_alias attribute, which I think is gcc specific
   #pragma warning (disable:1292)
   #endif
   G_LOCK(queue_lock);

   node = new_node();
   node->cell = cell;
   node->next = NULL;

   q->count++;
   /* nothing in the queue */
   if (q->head == NULL)
   {
      assert (q->tail == NULL);
      q->head = q->tail = node;
   }
   /* at least one thing in the queue */
   else
   {
      q->tail->next = node;
      q->tail = node;
   }
   #ifdef ICC
   #pragma warning (disable:1293) // icc complains about a may_alias attribute, which I think is gcc specific
   #pragma warning (disable:1292)
   #endif
   G_UNLOCK(queue_lock);
}
