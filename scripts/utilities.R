## General utility functions
# For example testing functions, debugging utilities, loading models...


load_library <- function(package_name, repos = "https://cloud.r-project.org/"){
    #' Load an R package via require()
    
    print(paste("Required", package_name))
    if (!require(package_name, character.only = TRUE)) install.packages(package_name, repos = repos, character.only = TRUE)
    
}
    
save_pdf <- function(fig, fig_name, fig_savedir, h = 4, w = 4){
    #' Save a figure as pdf, with easy adjustment of height etc. 
    #'
    #' @description
    #' Adds .pdf automatically at the end of fig_name. 
    #' Note: if you can't see a pdf figure in the file, try to use fig = print(fig) as input.
    
    pdf(paste(fig_savedir, paste(fig_name, ".pdf", sep = ""), sep = ""), height = h, width = w)
    print(fig)
    dev.off()
}
    
# Save as tiff
# Load/install required packages
#load_library("png")
#load_library("ggplot2")
#load_library("extrafont")
#font_import()
#loadfonts(device="pdf")
    
save_tiff <- function(fig, fig_name, fig_savedir, w=789, h=789, res = 300){
    #' Save a figure as a .tiff file
    #'
    #' Save as .tiff according to PLOS requirements at the time of writing. capabilities(): test if tiff TRUE, if not this does not work.
    #'
    #' @param fig A ggplot2 figure to be saved.
    #' @param fig_name Name of the figure, without .tif at the end.
    #' @param w Width of the figure in pixels
    #' @param h Height of the figure in pixels
    #' @param res Resolution
    
    tiff(paste(fig_savedir, paste(fig_name, ".tif", sep = ""), sep = ""), height = h, width = w, res = res, type = "cairo")
    print(fig + theme_cowplot(12, font_family = "Arial", font_size = 11))
    dev.off()
   
}
    
test_dims <- function(len1, len2, test_name = ""){
    #' Compare two lengths (unit test style)
    #'
    #' Input two lengths len1 and len2, compare them. If they are equal, pass.
    #'
    #' @param len1 Length 1
    #' @param len2 Length 2
    #' @return pass True or false, depending on whether the test passed
    #'
    #' @examples
    #' test_dims(dim(test)[1] - dim(mykrobe_data)[1], length(setdiff(blast_data$Isolate, mykrobe_data$Isolate)), test_name = "Mykrobe and BLAST merge")

    pass <- F
    print(paste("Testing", test_name))
    if (len1 == len2){
        print("PASSED: Same length.")
        pass <- T
    }
    else{
        print("FAILED: Lengths do not agree.")
        print(paste("len1:", len1))
        print(paste("len2:", len2))
    }
    return(pass)
}
    
fupper <- function(str) {
    #' Convert the first letter of a string/character str to uppercase.
    #'
    #' @param str string to modify
    substr(str, 1, 1) <- toupper(substr(str, 1, 1))
    return(str)
}

unit_test <- function(func_output, correct_output){
    #' Compare two outputs
    #'
    #' Compares the output of the function with the expected output and assesses the results
    #'
    #' @param func_output The output to be assessed.
    #' @param correct_output The correct/expected output.

    rownames(func_output) <- NULL
    rownames(correct_output) <- NULL 
    if (identical(func_output, correct_output)){ # NOTE: Identical is sensitive to the order of the output.
        print("PASSED.")
    }
    else{
        print("The output of your function does not match the correct output.")
        print("Correct output:")
        print(correct_output)
        print("Function output:")
        print(func_output)
    }
    
}

check_data <- function(df){
    #' Explore a data frame
    #'
    #' Check the number of NAs, draw histograms of numeric values, display counts of categorical values. A good first check for any data frame.
    #'
    #' @param df Data frame to be explored.
    
    for (col in colnames(df)){
        print(table(is.na(df[,col])))
        if(length(table(df[,col])) < 5){
            print(col)
            print(table(df[,col]))
        }else if(typeof(df[,col]) == "double" | typeof(df[,col]) == "integer"){
            hist(df[,col], main = col)
        }
    }
        
}

add_trailing_zeros <- function(char, n){
    #' Add zeros to the end of a string.
    #'
    #' @examples df <- apply(df, c(1,2), add_trailing_zeros, n = 1)
    #'
    #' @param char String to modify.
    #' @param n Number of zeros to add.
    #' @return char Modified character vector char.
    
    zeros <- paste(rep("0", n), collapse = "")
    
    if (!grepl("\\.", char)){
        char <- paste(char,zeros, sep = ".")
    }else{
        num_after_dot <- unlist(strsplit(gsub(x = char, pattern = "[0-9]*\\.", replacement = ""), ""))
        if (length(num_after_dot) < n){
            zeros <- paste(rep("0", n - length(num_after_dot)), collapse = "")
            char <- paste(char,zeros, sep = "")
        }
    } 
    
    char
}

#add_trailing_zeros("100.2") # should return 100.2
#add_trailing_zeros("100") # should return 100.0
    
    
# new and improved add trailing zeros
add_trailing_zeros_to_vec <- function(vec, n_extra_zeros = 0){
    #' Add trailing zeros to vector
    #'
    #' @desc Balances the number of digits after "." with trailing zeros such that all numbers have the same amount of digits.
    #'
    #' @param vec Vector of numbers
    #' @param n_extra_zeros Number of extra zeros. For example, one could want to balance the number of digits after dot for an entire data frame.
    #'
    #' @return vec Character vector with the added zeros. Note! Using as.numeric will lose the trailing zeros.
    
    vec <- as.character(vec)
    # find the digit with most characters AFTER THE DOT.
    ls <- sapply(strsplit(gsub(x = vec, pattern = "[-0-9]*\\.", replacement = ""), ""), length) # n digits after dot
    n <- max(ls)
    for (digit in vec){
        l <- length(strsplit(gsub(x = digit, pattern = "[-0-9]*\\.", replacement = ""), "")[[1]])
        if(l < n){ # less digits after dot than the max length
            zeros <- paste(rep("0",n - l), collapse = "") # find the number of zeros to add
            if (grepl("\\.", digit)){ 
                vec[vec == digit] <- paste(digit, zeros, sep = "")
            }
            else{ # if no dots in the digit, add . and the n zeros
                vec[vec == digit] <- paste(digit, paste(rep("0",n), collapse = ""), sep = ".")
            }
        }
        
        #vec[vec == digit] <- paste(digit, paste(rep("0", n_extra_zeros), collapse = "")
    }

    vec

}

## CI table

summary_table <- function(mod, pars, remove_lu_ci = TRUE){
    #' Create summary table with median, exp(median), CI
    #' Note: rounded up to two digits because signif produces too many trailing zeros.

    fit <- as.matrix(mod)

    qs <- apply(fit[,pars], 2, quantile, probs = c((1-0.95)/2, (1-(1-0.95)/2)))
    medians <- apply(fit[,pars], 2, median)

    df <- data.frame("median" = round(medians,2),
                     "hr" = round(exp(medians),2),
                     "l95" = round(qs[1,],2),
                     "u95" = round(qs[2,],2))

    df$median <- add_trailing_zeros_to_vec(df$median)
    df$hr <- add_trailing_zeros_to_vec(df$hr)
    df$l95 <- add_trailing_zeros_to_vec(df$l95)
    df$u95 <- add_trailing_zeros_to_vec(df$u95)
    df$CI <- paste(df$l95, df$u95, sep = ";")
    
    if (remove_lu_ci){
        df[, c("l95", "u95")] <- NULL
    }
    
    df
    
}

    
## Utilities for MRSA data
    
convert_v0_to_v1 <- function(v0, v0s = c("R0", "V1", "V2", "V3", "V4"), v1s = c("V1", "V2", "V3", "V4", NA)){
    #' Convert v0 visits (R0, V1, V2, V3, V4) to the corresponding next visit v1 (V1, V2, V3, V4, NA).
    #'
    #' Match v0 with the corresponding second visit. V4 is paired with an NA, since V4 is the last visit.
    #'
    #' @param v0 A vector of v0 visits of interest, for example c("R0", "V2").
    #' @param v0 A vector of all v0 visits.
    #' @param v1s A vector of all v1 visits
    #' @return Returns the next visits v1 corresponding to v0.

    return(v1s[match(v0,v0s)])   
}
    
convert_v1_to_v0 <- function(v1, v0s = c(NA, "R0", "V1", "V2", "V3", "V4"), v1s = c("R0","V1", "V2", "V3", "V4")){
    #' Same as `convert_v0_to_v1` but convert v1 to v0
     v0s[match(v1, v1s)]
}

check_no_resistant_obs <- function(data, covariates){
    #' Check for cases where there are no resistant isolates.
    #'
    #' For example, in site-specific data, decolonization group + throat has no Tetracycline or Rifampicin resistant 
    #'   observations, which in turn halts stan_surv. Exclude these from covariates.
    #'
    #' @param data Survival data frame.
    #' @param covariates Covariates for the survival model of interest.
    #' @return covariates Covariates without missing resistances.
    
    for (treatment in treatments){
        if (sum(na.omit(data[,treatment])) == 0){ # no resistant cases
            print(paste("No resistant cases, antibiotic:", treatment))
            covariates <- covariates[covariates != treatment]
        }
    }
    
    return(covariates)
}
    
    
check_same_visit_date_for_all_sites <- function(mykrobe_data){
    #' Check same visit date for all sites assumption
    #'
    #' @description
    #' Can it be assumed, based on mykrobe_data, that during a given visit (for example, R0),
    #'  all sites were swabbed during the same day? This function warns about cases where this does not hold and
    #'  returns a vector of the hosts that do not satisfy the assumption.
    #'
    #' @param mykrobe_data Unpreprocessed MRSA dataset
    #' @return A vector of hosts that do not satisfy the assumption
    
    hosts_to_check <- c()
    
    for (host in unique(mykrobe_data$ID)){
        hostdata <- mykrobe_data[which(mykrobe_data$ID == host),]
        visits <- unique(hostdata$Visit..)

        for (v in visits){
            tab <- table(hostdata[which(hostdata$Visit.. == v),"Cx.Site"], hostdata[which(hostdata$Visit.. == v),"Cx.Date"])

            if (dim(tab)[2] != 1){
                print(paste("WARNING! Different sites have different dates. Check host:", host))
                print(paste("Visit:", v, sep = " "))
                print(tab)
                hosts_to_check <- c(host, hosts_to_check)
            }
        }
    }
    
    hosts_to_check
}

    
read_model_set <- function(savedir, mod_id, mods_of_interest, is_site = FALSE){
    #' Load the pretrained models in savedir
    #'
    #' @param savedir Save directory (as a string) where the pretrained models are saved.
    #' @param mods_of_interest Character vector of names for the models of interest, like host_idstrain (consists of random effects)
    #' @param mod_id 0 (no covariates), 1 (separate covariates) or 2 (all covariates), depending on the model type.
    #' @return mods Loaded models.
    
    mod_id <- paste("mod", mod_id, sep = "")
    if (is_site){
        mod_id <- paste("site_", mod_id, sep = "")
    }
    
    # Construct model filenames based on mod_id and mods_of_interest
    fd <- paste(paste(paste(mod_id, mods_of_interest, sep = "_"), "D", sep = ""), ".Rdata", sep = "")
    fe <- paste(paste(paste(mod_id, mods_of_interest, sep = "_"), "E", sep = ""), ".Rdata", sep = "")
    filenames <- paste(savedir, c(fd, fe), sep = "") # Add save directory
    
    # Load the models:
    mods <- lapply(filenames, readRDS)
    
    # Remove model id from the name (load only 1 model type at a time):
    names(mods) <- gsub("2","",gsub(".Rdata", "", substr(filenames, nchar(savedir) + 1, max(nchar(filenames)))))
    names(mods) <- gsub("1", "", names(mods))
    names(mods) <- gsub("0", "", names(mods))
    
    return(mods)
    

}
    
