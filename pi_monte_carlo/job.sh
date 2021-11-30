#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=0:20:00 -q dssc

module load openmpi-4.1.1+gnu-9.3.0
cd /u/dssc/tfonda/fast/hpc/git/pi_monte_carlo
mpirun -np 8 --oversubscribe pi_monte_carlo ${n}
