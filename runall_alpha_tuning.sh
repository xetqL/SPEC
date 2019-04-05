#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for A in "$@"
do
sbatch --constraint=E5-2630V4 -p $partition --output alpha_N1_5_32_"$A"_%j.out -n 32 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 32 1 1000 450 12 1 $A
sbatch --constraint=E5-2630V4 -p $partition --output alpha_N1_5_64_"$A"_%j.out -n 64 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 64 1 1000 450 12 1 $A
sbatch --constraint=E5-2630V4 -p $partition --output alpha_N1_5_128_"$A"_%j.out -n 128 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 128 1 1000 450 12 1 $A
sbatch --constraint=E5-2630V4 -p $partition --output alpha_N1_5_256_"$A"_%j.out -n 256 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 256 1 1000 450 12 1 $A
sbatch --constraint=E5-2630V4 -p $partition --output alpha_N1_5_512_"$A"_%j.out -n 512 -c 1 -J SPEC -t $time spec.slurm bin/SPEC_5 512 1 1000 450 12 1 $A
done 
