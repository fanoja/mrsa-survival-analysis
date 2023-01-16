#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --mem=15G
#SBATCH -c 4
#SBATCH --array=1-8
#SBATCH --output=outputs/slurm_raneff_parallel_%A_%a.out

cd scripts/cluster/
source r4env.sh

n=$SLURM_ARRAY_TASK_ID

raneff=`sed -n "${n} p" raneff.txt`

Rscript mod2_parallel_sites.R $raneff

