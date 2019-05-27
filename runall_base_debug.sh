#!/bin/bash
set -o xtrace
partition=$1
time="0-00:15:00"
shift
for P in "$@"
do

sbatch -p $partition --output %j_N?_ULBA_"$P".out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_ULBA $P 1 2000 450 12 ? 0.3

sbatch -p $partition --output %j_N?_STD_"$P".out -n $P -c 1 -J SPEC -t $time spec.slurm bin/SPEC_STD $P 1 2000 450 12 ? 0.3

done 
