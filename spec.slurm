#!/bin/sh 
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
let sqr=$(($4**2))
rm log/*.log
rm out.npz
ml purge; ml GCC/6.4.0-2.28 OpenMPI/2.1.1
mpirun --oversubscribe -np $2  $1 --yprocs $2 --xprocs $3 -c $sqr -t $5 -s $6 -N $7 -a $8
#mpirun $1 $2 $3 $sqr $5 $6

