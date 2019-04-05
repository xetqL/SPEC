#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do
sbatch -p $partition --output N10_5_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 $P 1 1000 100 12 10 0.18
sbatch -p $partition --output N10_3_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_3 $P 1 1000 100 12 10 0.18
sbatch -p $partition --output N10_2_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_2 $P 1 1000 100 12 10 0.18
done 
