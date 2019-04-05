#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do
sbatch --constraint=E5-2630V4 -p $partition --output N5_5_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 $P 1 1000 300 12 5 0.23
sbatch --constraint=E5-2630V4 -p $partition --output N5_3_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_3 $P 1 1000 300 12 5 0.23
#sbatch --constraint=E5-2630V4 -p $partition --output N5_2_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_2 $P 1 1000 100 12 5 0.18
done 
