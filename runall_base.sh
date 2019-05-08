#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do
sbatch --constraint=E5-2630V4 -p $partition --output N?_1_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_1 $P 1 1000 450 12 ? 0.2
#sbatch --constraint=E5-2630V4 -p $partition --output N?_3_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_3 $P 1 1000 450 12 ? 0.45
#sbatch --constraint=E5-2630V4 -p $partition --output N?_0_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_NO_LB $P 1 1000 100 12 ? 0.45
done 
