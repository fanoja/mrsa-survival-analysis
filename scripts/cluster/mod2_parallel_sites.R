## Parallellization by random effect: each random effect model pair (d/e) runs in parallel.
## Site-specific models
## Probably slower, since each batch contains now 8 models.
options(mc.cores = parallel::detectCores())
source(paste(getwd(), "get_model_functions.R", sep = "/"))

model_id <- "site_mod2"
savedir <- paste(getwd(), "/saved_models/interval_censored/model2/site_models/", sep = "")
modelfun <- get_model2

args = commandArgs(trailingOnly=TRUE)

raneff <- unlist(strsplit(args[1], split = ",")) # Can now handle more than one random effect, if separated by ",".
savefile <- paste(model_id, paste(raneff, collapse = ""), sep = "_")

if (raneff == "fixed"){
  #savefile <- paste(model_id, "fixed", sep = "_") # Reasonable file name for fixed effects.
  raneff <- c()
}

adapt_delta <- as.numeric(args[2]) # Second term is the adapt_delta; same for both models.
print(paste("Raneff:", raneff))
print(paste("Adapt_delta:", adapt_delta))

mod_d <- get_model_for_all_sites("D", modelfun = modelfun, params = treatments, raneff = c(raneff), adapt_delta = c(0.9999, 0.9999999, 0.9999, 0.9999), iter =  c(10000, 10000, 10000, 10000))
mod_e <- get_model_for_all_sites("E", modelfun = modelfun, params = treatments, raneff = c(raneff), adapt_delta = c(0.9999, 0.9999999, 0.9999, 0.9999), iter =  c(10000, 10000, 10000, 10000))

saveRDS(mod_d, file=paste(savedir, paste(paste(savefile, "D", sep = ""), ".Rdata", sep = ""), sep = ""))
saveRDS(mod_e, file=paste(savedir, paste(paste(savefile, "E", sep = ""), ".Rdata", sep = ""), sep = ""))