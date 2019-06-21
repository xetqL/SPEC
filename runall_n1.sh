#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do

sbatch -p $partition --output %j_N1_ULBA_"$P".out -n $P -c 1 -J SPEC_ULBA -t $time spec.slurm bin/SPEC_ULBA $P 1 3000 200 12 1 0.3

sbatch -p $partition --output %j_N1_STD_"$P".out -n $P -c 1 -J SPEC_STD -t $time  spec.slurm bin/SPEC_STD $P 1 3000 200 12 1 0.3

done 
