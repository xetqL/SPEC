#!/bin/bash
set -o xtrace
partition=$1
shift
for P in "$@"
do
sbatch  --mincpus 16 -p $partition --output debug-3_"$P"_%j.out -n $P -c 1 -J SPEC -t 0-00:08:00 spec-new.slurm bin/SPEC_Debug3 $P 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 5_64_%j.out -n 64 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_5 64 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 5_128_%j.out -n 128 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_5 128 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 5_256_%j.out -n 256 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_5 256 1 1000 100 12
#sbatch  --mincpus 16 -p $partition --output debug-3_"$P"_%j.out -n $P -c 1 -J SPEC -t 0-00:08:00 spec.slurm bin/SPEC_3 $P 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 3_64_%j.out -n 64 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_3 64 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 3_128_%j.out -n 128 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_3 128 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 3_256_%j.out -n 256 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_3 256 1 1000 100 12
#sbatch  --mincpus 16 -p $partition --output debug-2_"$P"_%j.out -n $P -c 1 -J SPEC -t 0-00:08:00 spec.slurm bin/SPEC_2 $P 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 2_128_%j.out -n 128 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_2 128 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 2_64_%j.out -n 64 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_2 64 1 1000 100 12
#sbatch --mincpus 16 -p shared --output 2_32_%j.out -n 32 -c 1 -J SPEC -t 0-00:30:00 spec.slurm bin/SPEC_2 32 1 1000 100 12
done 
