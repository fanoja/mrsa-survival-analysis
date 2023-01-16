# Run model 2 with CV
# Note: works only for the pooled over sites model!

library("rstanarm")

# First, load the data and premade folds:
data <- readRDS("strain_surv_data/data.RData")
data_sites <- readRDS("strain_surv_data/data_sites.RData")
treatments <- readRDS("strain_surv_data/treatments.RData")

get_new_folds <- T # Do you wish to create new folds using get_folds function?

print("Loaded data.")

get_folds <- function(data, seed = 02122021){
  #' Create and save new folds
    
  set.seed(seed)
  
  folds_strain_d <- loo::kfold_split_stratified(K = 10, x = data[which(data$ARM == 1),"strain"])
  folds_strain_e <- loo::kfold_split_stratified(K = 10, x = data[which(data$ARM == 0),"strain"])
  
  saveRDS(folds_strain_d, file = "folds_strain_d.RData")
  saveRDS(folds_strain_e, file = "folds_strain_e.RData")
  
  print("Folds created and saved.")
}

if (get_new_folds){
  get_folds(data)
}


folds_strain_d <- readRDS("folds_strain_d.RData")
folds_strain_e <- readRDS("folds_strain_e.RData")

options(mc.cores = parallel::detectCores())

# Read raneff_formula, savefile and arm from command line (saved to txt files)
args = commandArgs(trailingOnly=TRUE)
print(args)
print(length(args))

raneff <- unlist(strsplit(args[1], split = ","))

raneff_formula <-gsub(",", "",toString(paste(paste('+ (1', raneff, sep = '|'), ')', sep = '')))

if(raneff == "fixed"){
  raneff_formula <- ''
}

savefile <- args[2]
arm <- as.numeric(args[3])

print(args[2])
print(args[3])

data <- data[which(data$ARM == arm),]


# 1 = decolonization, 0 = education
if (arm == 1){
    folds_strain <- folds_strain_d
}

if(arm == 0){
    folds_strain <- folds_strain_e    
}

pooled_formula <- paste(paste('Surv(time = t0, time2 = delta_t, event = y, type = "interval") ~', gsub(',', ' + ', toString(treatments)), sep = ''), raneff_formula, sep = '')

#pooled_formula <- paste(pooled_formula, raneff_formula, sep = "")  
  
print(pooled_formula)

# Run model 2 (all antibiotics). Iter and adapt_delta higher, because with both strain and host raneff,
# there would be divergent transitions issues and ESS tail etc problems.
print(paste("Generating model: ", raneff_formula))
print(paste("ARM", arm))

all <- stan_surv(data = data[which(data$ARM == arm),], formula = as.formula(pooled_formula),
                 basehaz = "exp", refresh = 0,
                 prior_intercept = normal(0,20), prior = normal(0,2.5),
            iter = 7500, adapt_delta = 0.999) #iter = 10000, adapt_delta = 0.99)

# Kfold CV trick for stan_surv:
# Courtesy of John Baums, https://github.com/stan-dev/rstanarm/issues/473
f <- eval(parse(
  text=sub(
    'fit_k_call$subset <- eval(fit_k_call$subset)', 
    'fit_k_call$subset <- if(is.stansurv(x)) NULL else eval(fit_k_call$subset)', 
    deparse(rstanarm:::kfold.stanreg), fixed=TRUE)
), envir=environment(rstanarm:::kfold.stanreg))

assignInNamespace('kfold.stanreg', value=f, ns='rstanarm')

# Run kfold:
print("Running 10-fold CV...")
kfold_data <- kfold(all, K=10, folds = folds_strain)
saveRDS(kfold_data, file = savefile)
print("Done.")

       
 