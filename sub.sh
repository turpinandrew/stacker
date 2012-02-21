#!/bin/bash

#PBS -l walltime=15:00:00
#PBS -N Z
#PBS -l procs=1
#PBS -l pvmem=15GB

#  qsub -I -X -l pvmem=20gb  # interactive

module load gsl-gcc/1.14 gtk-gcc/2.22.0 

./stack
