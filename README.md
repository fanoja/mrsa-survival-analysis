# MRSA Bayesian Survival Analysis

Analysis code for the manuscript *Impact of Antibiotic Resistance on MRSA Cleareance Using Bayesian Survival Analysis.*


## Requirements

R4.0.5 

Packages:
- rstanarm
- bayesplot
- igraph
- network
- reshape2
- ggplot2
- gridExtra
- cowplot

To gain access to stan_surv function for survival analysis, install the developlment version of [rstanarm from GitHub](https://github.com/stan-dev/rstanarm). For this, packages `rstan` and `remotes` are also required. If you have installation issues, revert to R4.0.5 as suggested in this [issue](https://github.com/stan-dev/rstanarm/issues/500#issuecomment-1203904085).

For package management, [renv](https://rstudio.github.io/renv/index.html) is used. Check an introduction to renv [here](https://rstudio.github.io/renv/articles/renv.html).

## Structure

### Directories
- `data/`: Distance matrices, MRSA metadata with resistances and clearance
- `results/`: Figures and tables from the article as generated by this code
- `saved_models/`: Save fitted models here for visualization. Note! These models are not in the git repository at the moment due to file size limitations. Please contact us if you wish to obtain these files.
- `kfold/`: Results for K-fold Cross-Validation
- `scripts/`: R scripts necessary for the analysis.
- `scripts/cluster/`: Scripts to fit the Bayesian model on a computing cluster.
- `strains/`: Files associated with assigning strains to the isolates (see [BaeMBac](https://github.com/mjarvenpaa/bacterial-colonization-model))


### Data
- Data in survival analysis compatible format should be included in the directory `data/surv_data`. The following files are required:
    - `data.RData`: data "pooled" over sites
    - `data_ast.RData`: same as above, but with phenotypic resistance for mupirocin.
    - `data_sites.RData`: site-specific survival data
    - `treatments.RData`: contains the names of the antimicrobials of interest

### Jupyter notebooks
- `MRSA Clearance.ipynb`: Generation of all figures in the manuscript, instructions to running the analysis with examples.
- `revision.ipynb`: Notebook containing additional results for revision.

### Source code
- surv_preprocessing_strains.R: Preprocess the data to a survival analysis format
- mdata_preprocessing.R: Preprocessing the data prior to creating a survival analysis compatible dataset
- mcmc_plotting_functions.R: Various functions to visualize the data and results
- get_model_functions.R: Various functions for running the models 
- utilities.R: Utility functions

## Recommended reading

The `stan_surv` function used is part of the `rstanarm` package. A detailed explanation can be found in "Bayesian Survival Analysis using the rstanarm R package" (Brilleman et al. 2020),  	
https://doi.org/10.48550/arXiv.2002.09633

`BaeMBac` software used to assign strains to isolates within-host. A detailed explanation is provided in https://doi.org/10.1371/journal.pcbi.1006534

