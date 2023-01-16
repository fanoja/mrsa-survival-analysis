# K-fold CV

This directory contains scripts and a notebook example for running 10-fold CV for fitted models. See `kfold.ipynb` for a tutorial that can be run locally.

Scripts:
- `kfold_triton.sh`: A batch script for running k-fold on a computing cluster.
- `kfoldCV.R`: An R script for running K-fold CV. Same as in `kfold.ipynb`.

Other files:
- `folds_strain_d.RData` and `folds_strain_e.RData`: Pregenerated folds for running k-fold CV.
- `arm.txt`, `raneff.txt`, `savefile.txt`: Helper files for using an array job on a computing cluster (running kfold on multiple models, for example with and without random effects, at the same time).

Directories:
- `res/`: Results are saved here.
- `strain_surv_data/`: Survival data files.