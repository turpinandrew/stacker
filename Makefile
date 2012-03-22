GTK_INCLUDES = `pkg-config --cflags gtk+-3.0`
GTK_LIBS = `pkg-config --libs gtk+-3.0`

DEFS = -DTHREADS=10    # number of EXTRA threads to use

#for gcc
CC = gcc
#CFLAGS = -O3 -Wall -std=c99 $(DEFS)
CFLAGS = -g -Wall -std=c99 $(DEFS)

# for icc
   # remark #981: operands are evaluated in unspecified order
   # remark #869: parameter "widget" was never referenced
#CC = icc
#CFLAGS = -O3 -Wall -std=c99 -DICC -wd981 -wd869 $(DEFS)

#CPPFLAGS =
#LD_FLAGS = $(GTK_LIBS) -lm -lgsl -lgslcblas
#LD_FLAGS = -lm -lgsl -lgslcblas
HDRS = main.h density.h queue.h setup.h types.h vector.h
OBJS = main.o density.o queue.o setup.o vector.o
SRCS = main.c queue.c density.c setup.c vector.c
EXE = stack

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) -lm -lgsl -lgslcblas $(GTK_LIBS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $< $(CPPFLAGS) $(GTK_INCLUDES)

clean:
	/bin/rm -fr $(OBJS)

clobber: clean
	/bin/rm -fr $(EXE)

main.o: main.c queue.h Makefile setup.h types.h vector.h
queue.o: queue.c queue.h Makefile
density.o: density.c density.h Makefile
setup.o: setup.c setup.h Makefile queue.h density.h types.h
vector.o: vector.h Makefile
