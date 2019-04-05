#!/bin/bash
set -o xtrace
partition=$1
time="0-00:13:00"
shift
for P in "$@"
do
#sbatch  -p $partition --output N1_0_"$P"_%j.out -N 2 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_NO_LB $P 1 1000 500 12 1 0.23
sbatch --constraint=E5-2630V4 -p $partition --output N1_5_"$P"_%j.out -N 2 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 $P 1 1000 500 12 1 0.18
sbatch --constraint=E5-2630V4 -p $partition --output N1_3_"$P"_%j.out -N 2 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_3 $P 1 1000 500 12 1 0.18

sbatch -p debug --output N1_5_"$P"_%j.out -N 2 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 $P 1 1000 410 12 1 0.4

sbatch -p debug --output N1_3_"$P"_%j.out -N 2 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_3 $P 1 1000 410 12 1 0.4

#sbatch --constraint=E5-2630V4 -p $partition --output N1_2_"$P"_%j.out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_2 $P 1 1000 100 12 1 0.18
done 
