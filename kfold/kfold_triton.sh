#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH --mem=15G
#SBATCH -c 4
#SBATCH --array=1-8
#SBATCH --output=outputs/slurm_kfold_parallel_%A_%a.out

cd $WRKDIR/MRSA/stan_surv
source r4env.sh

cd $WRKDIR/MRSA/stan_surv/kfold

n=$SLURM_ARRAY_TASK_ID

raneff_formula=`sed -n "${n} p" raneff_formula.txt`
savefile=`sed -n "${n} p" savefile.txt`
arm=`sed -n "${n} p" arm.txt`

Rscript kfoldCV.R $raneff_formula $savefile $arm
