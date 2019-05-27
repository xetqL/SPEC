#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do

sbatch --constraint=E5-2630V4 -p $partition --output N1_ULBA_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_ULBA $P 1 2000 450 12 1 0.3

sbatch --constraint=E5-2630V4 -p $partition --output N1_STD_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_STD $P 1 2000 450 12 1 0.3

done 
