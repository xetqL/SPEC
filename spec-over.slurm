#!/bin/sh 
#SBATCH --cpus-per-task=1
#SBATCH --output=%J.out
#SBATCH --exclusive
let sqr=$(($3**2))
rm log/*.log
rm out.npz
ml purge; ml GCC/6.4.0-2.28 OpenMPI/2.1.1
mpirun ./bin/SPEC $1 $2 $sqr $4 $5

