#!/bin/bash
#PBS -l nodes=1:ppn=4,walltime=0:20:00 -q dssc

module load openmpi
cd /u/dssc/tfonda/fast/hpc/git/ring
mpirun -np 4 --oversubscribe ring
