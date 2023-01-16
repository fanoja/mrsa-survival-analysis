#!/bin/bash

#SBATCH --time=1:30:00
#SBATCH --mem=5G
#SBATCH -c 4
#SBATCH --array=1-8
#SBATCH --output=outputs/slurm_raneff_parallel_%A_%a.out

cd scripts/cluster/
source r4env.sh

n=$SLURM_ARRAY_TASK_ID

raneff=`sed -n "${n} p" raneff.txt`
adapt_delta=`sed -n "${n} p" adapt_delta.txt`

Rscript mod2_parallel.R $raneff $adapt_delta

