## Parallellization by random effect: each random effect model pair (d/e) runs in parallel.
options(mc.cores = parallel::detectCores())
source(paste(getwd(), "get_model_functions.R", sep = "/"))

model_id <- "mod2"
savedir <- paste(getwd(), "/saved_models/interval_censored/model2/", sep = "")
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

mod_d <- modelfun(data[which(data$ARM == 1),], params = treatments, raneff = c(raneff), adapt_delta = adapt_delta, iter = 7500)
mod_e <- modelfun(data[which(data$ARM == 0),], params = treatments, raneff = c(raneff), adapt_delta = adapt_delta, iter = 7500)

saveRDS(mod_d, file=paste(savedir, paste(paste(savefile, "D", sep = ""), ".Rdata", sep = ""), sep = ""))
saveRDS(mod_e, file=paste(savedir, paste(paste(savefile, "E", sep = ""), ".Rdata", sep = ""), sep = ""))