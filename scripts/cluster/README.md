# README

If you wish to save the previous models and create a new folder for new analysis, create the following directories `saved_models`
- saved_models/interval_censored/model2/
- saved_models/interval_censored/model2/site_models
  
To run on a computing cluster in parallel, input the following command:

For `sbatch mod2_raneff_parallel.sh`

For site-specific models: `sbatch mod2_raneff_parallel_sites  .sh`

Note that the scripts might not work as is on your computing cluster, but can be used as inspiration for adaptation.

`outputs/` is used to store log files.

`data/strain_surv_data/` has the same data files as the main directory. This way the `cluster/` directory can be used independently.

