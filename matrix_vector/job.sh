#!/bin/bash
#PBS -l nodes=1:ppn=4,walltime=0:10:00 -q dssc

module load openmpi-4.1.1+gnu-9.3.0
cd /u/dssc/tfonda/fast/hpc/matrix_vector
mpirun -np 4 --oversubscribe matrix_vector
