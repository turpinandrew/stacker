#ifndef QUEUE_H
#define QUEUE_H

#include <glib.h>

G_LOCK_EXTERN(queue_lock);

typedef struct qnode Queue_node;
typedef struct queue Queue;

struct qnode
{
   void *cell;
   Queue_node *next;
};

struct queue
{
   Queue_node *head;
   Queue_node *tail;
   int count;
};

extern void insert_last(Queue *, void *);
extern void *remove_first(Queue *);
extern Queue *new_empty_queue(void);

#endif
