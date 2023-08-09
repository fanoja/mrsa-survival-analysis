## Functions for running survival models ##
# By default, uses interval censoring.

# Model 0: no covariates
# Model 1: each resistance in a separate model
# Model 2: All resistances in a single model. Used in the manuscript.

# Load library via require
load_library <- function(package_name, repos = "https://cloud.r-project.org/"){
    print(paste("Required", package_name))
    if (!require(package_name, character.only = TRUE)) install.packages(package_name, repos = repos, character.only = TRUE)
    
}

load_library("rstanarm")

set.seed(150545)

strain_based <- TRUE

if(strain_based){
    # Read the data required for the model in question
    data <- readRDS("data/surv_data/data.RData")
    data_sites <- readRDS("data/surv_data/data_sites.RData")
    treatments <- readRDS("data/surv_data/treatments.RData")
}

check_no_resistant_obs <- function(data, covariates){
    #' Check for cases where there are no resistant isolates.
    #'
    #' For example, in site-specific data, decolonization group + throat has no Tetracycline or Rifampicin resistant 
    #' observations, which in turn halts stan_surv. Exclude these from the covariates.
    #'
    #' @param data Survival data frame.
    #' @param covariates Covariates to check, here antibiotics of interest.
    #'
    #' @return covariates The names of antibiotics that have at least one resistant observation.
    
    for (treatment in treatments){
        if (sum(na.omit(data[,treatment])) == 0){ # no resistant cases
            print(paste("Non-resistant cases, antibiotic:", treatment))
            covariates <- covariates[covariates != treatment]
        }
    }
    
    return(covariates)
}


get_model0 <- function(data, params = NULL, raneff = c(), formula_start = "Surv(time = t0, time2 = delta_t, event=y, type = 'interval') ~", adapt_delta = 0.95, prior_intercept = normal(0,20), prior = normal(0,2.5), iter = 5000){
    #' Get a model with only the intercept term, no antibiotics.
    #'
    #' @param params Parameters of interest.
    #' @param raneff A vector containing the names of random effect variables, such as "strain".
    #' @param formula_start The beginning of the formula that is later passed to stan_surv. Decide the type of censoring here.
    #'
    #' @return mod0 Fitted model.
    
    formula <- paste(formula_start, "1", sep = "") #Surv(delta_t, y) ~ 1"
    # Note that 1 + (1|host_id) = (1|host_id) in lme4 formula syntax
    
    if (length(raneff) != 0){
        re <-gsub(",", "", toString(paste(paste("+(1|", raneff, sep = ""), ")", sep = "")))
        formula <- paste(formula, re)
    }
    
    mod0 <- stan_surv(formula = as.formula(formula), data = as.data.frame(data), basehaz = "exp", refresh = 0, iter = iter, adapt_delta = adapt_delta, prior_intercept = prior_intercept)
    
    return(mod0)
 }

get_model1 <- function(data, params, raneff = c(), formula_start = "Surv(time = t0, time2 = delta_t, event=y, type = 'interval') ~", adapt_delta = c(0.95), prior_intercept = normal(0,20), prior = normal(0,2.5), iter = 5000){
    #' Get separate models (Model 1) for each parameter/antibiotic in params.
    
    params <- check_no_resistant_obs(data, params)
    
    if (length(adapt_delta) != length(params)){
        adapt_delta <- rep(adapt_delta[1], length(params)) # Repeat the first item in adapt_delta for length(params) times.
       }
    
    mod1s <- list()
    
    # Process random effects:
    re <- ""
    if (length(raneff) != 0){
        re <-gsub(",", "", toString(paste(paste("+(1|", raneff, sep = ""), ")", sep = ""))) # add random effects to formula
    }
             
    for (t in 1:length(params)){
        formula <- as.formula(paste(paste(formula_start, params[t], sep = ""), re, sep = ""))
        print(formula)
        mod1s[[paste("mod", t, sep = "")]] <- stan_surv(formula = as.formula(formula), data = as.data.frame(data), basehaz = "exp", refresh = 0, iter = iter, adapt_delta = adapt_delta[t])

    }
    
    return(mod1s)
    
 }

get_model2 <- function(data, params, raneff = c(), formula_start = "Surv(time = t0, time2 = delta_t, event=y, type = 'interval') ~", adapt_delta = 0.95, prior = normal(0, 2.5), prior_intercept = normal(0,20), iter = 5000){
    #' Get a joint model for all params (antibiotics). This is the model used in the manuscript.

    # Process random effects:
    re <- ""
    if (length(raneff) != 0){
    re <-gsub(",", "", toString(paste(paste("+(1|", raneff, sep = ""), ")", sep = "")))
    }
    
    params <- check_no_resistant_obs(data, params)
    
    formula <- paste(paste(formula_start, gsub(",", " + ", toString(params)), sep = ""),re, sep = "")
    print(formula)
    mod2 <- stan_surv(formula = as.formula(formula), data = as.data.frame(data), basehaz = "exp", refresh = 0, iter = iter, adapt_delta = adapt_delta, prior_intercept = prior_intercept, prior = prior)
    
    return(mod2)
}

get_model_for_all_sites <- function(arm, modelfun, params = c(), raneff = c(), formula_start = "Surv(time = t0, time2 = delta_t, event=y, type = 'interval') ~", prior = normal(0, 2.5), prior_intercept = normal(0,20), adapt_delta = c(0.95, 0.95, 0.95, 0.95), iter = c(5000, 5000, 5000, 5000)){
    #' Get a model for all sites (nares, throat, skin, wound).
    #'
    #' @param arm "D" or "E", decolonization or education.
    #' @param modelfun Model function of interest (get_model0, get_model1, get_model2)
    #'
    #' @return A list of site-specific models fitted using modelfun.
 
    # Note: Input four adapt_deltas, one for each site. Order: nares, throat, skin, wound
    # Note! Assumes the same prior for all sites.
    
    if (arm == "D"){
        arm  <- 1
    }
    if (arm == "E"){
        arm <- 0
    }
    
    data_nares <- data_sites[which(data_sites$site == "nares" & data_sites$ARM == arm),]
    data_throat <- data_sites[which(data_sites$site == "throat" & data_sites$ARM == arm),]
    data_skin <- data_sites[which(data_sites$site == "skin" & data_sites$ARM == arm),]
    data_wound <- data_sites[which(data_sites$site == "wound" & data_sites$ARM == arm),]
           
    
    print("Nares")
    params_n <- check_no_resistant_obs(data_nares, params)
    mod_n <- modelfun(data_nares, params = params_n, raneff = raneff, formula_start = formula_start, prior = prior, prior_intercept = prior_intercept, adapt_delta = adapt_delta[1], iter = iter[1])
    print("Throat")
    params_t <- check_no_resistant_obs(data_throat, params)
    mod_t <- modelfun(data_throat, params = params_t, raneff = raneff, formula_start = formula_start, prior = prior, prior_intercept = prior_intercept, adapt_delta = adapt_delta[2], iter = iter[2])
    print("Skin")
    params_s <- check_no_resistant_obs(data_skin, params)
    mod_s <- modelfun(data_skin, params = params_s, raneff = raneff, formula_start = formula_start, prior = prior, prior_intercept = prior_intercept, adapt_delta = adapt_delta[3], iter = iter[3])
    print("Wound")
    params_w <- check_no_resistant_obs(data_wound, params)
    mod_w <- modelfun(data_wound, params = params_w, raneff = raneff, formula_start = formula_start, prior = prior, prior_intercept = prior_intercept, adapt_delta = adapt_delta[4], iter = iter[4])
    
    return(list("mod_n" = mod_n, "mod_t"= mod_t, "mod_s" = mod_s, "mod_w" = mod_w))
    
}

