#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
seed=$RANDOM
for P in "$@"
do

sbatch --constraint=E5-2630V4 -p $partition --output zoltan_paper_%j_N?_ULBA_"$P".out -n $P -c 1 -J zSPEC_ULBA -t $time spec.slurm bin/SPEC_ULBA $P 1 1000 450 $seed ? 0.2

sbatch --constraint=E5-2630V4 -p $partition --output zoltan_paper_%j_N?_STD_"$P".out -n $P -c 1 -J  zSPEC_STD -t $time spec.slurm bin/SPEC_STD -t $time spec.slurm bin/SPEC_STD $P 1 1000 450 $seed ? 0.2

done 
