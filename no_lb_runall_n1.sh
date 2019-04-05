#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do
sbatch --constraint=E5-2630V4 -p $partition --output N1_1_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_1 $P 1 1000 100 12 1 0.18
sbatch --constraint=E5-2630V4 -p $partition --output N1_0_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_NO_LB $P 1 1000 100 12 1 0.18
done 
