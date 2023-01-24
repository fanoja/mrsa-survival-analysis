## Functions for running BaeMBac on the data

# Note that you will need the BaeMBac code (https://github.com/mjarvenpaa/bacterial-colonization-model) in this directory for running this.

if (!require("latex2exp")) install.packages("latex2exp", repos = "https://cloud.r-project.org/")
if (!require("coda")) install.packages("coda", repos = "https://cloud.r-project.org/")
if (!require("BayesianTools")) install.packages("BayesianTools", repos = "https://cloud.r-project.org/")
if (!require("MASS")) install.packages("MASS", repos = "https://cloud.r-project.org/")


run_baembac <- function(arm, root = "BaeMBac/", p_sample = 0.1){
    #' Run BaeMBac on D (decolonization) or E (education) arm. 
    #'
    #' Requires that dfst5 and dfst8 are loaded. Assigns same strain indicator based on the probability of being same strain (p > 0.5 -> same strain). Saves the result as a .csv file.
    #'
    #' @param arm D for decolonization and E for education.
    #' @param root Folder where BaeMBac source code is.
    
    
    source(paste(root, "code/mixture_model_demos.R", sep = ""))
    #arm <- "E"
    #p_sample <- 0.1

    df5 <- dfst5[which(dfst5$host1 == dfst5$host2 & dfst5$arm1 == arm & dfst5$arm2 == arm),] # One observation had a distance of 0.
    df8 <- dfst8[which(dfst8$host1 == dfst8$host2 & dfst8$arm1 == arm & dfst8$arm2 == arm ),] # & dfst8$arm1 == arm & dfst8$arm2 == arm 

    df <- rbind(df5, df8)

    #n <- round(dim(df)[1]*p_sample)
    #df <- df[sample(nrow(df), n), ]

    # Get the SNP distance d_i and time difference in bacterial generations t_i
    d_is <- df$value
    t_is <- as.integer(df$timediff)*24*60/90 # Time in generations
    t_is[t_is == 0] <- 1

    save(d_is, t_is, file = paste(root, "data_CLEAR_minimal/mrsa_clear_testdata_armE.RData", sep = "")) # Replace the original data with the new data.

    compute.true.data(root) # Uncomment to run. Warning: will take time.

    # Add a column to df:
    load(paste(root, "simulation_outputs_mixture_model/mixture_true_data_E_final_pr/mixture_fitting_test4.RData", sep = ""))
    mean_z1 <- rowMeans(zs[,1,]) # posterior means of same strain cases
    samestrain <- ifelse(mean_z1 > 0.5,1,0) # if mean > 0.5, consider same strain
    df$colgroup <- samestrain
    df$baembac_ps <- mean_z1

    write.csv(df, paste(paste(paste(root, "bmbdf", sep = ""), arm, sep = ""), ".csv", sep = ""))
    
}

add_fitted_colgroups_baembac <- function(dfst, bmbdf_file = "/u/50/ojalaf2/unix/Dropbox/Mrsa_clear/results/BaeMBac/bmbdf.csv"){
    #' Add same/different strain indicators.
    #'
    #' If BaeMBac was run only on part of the data, use this function to add the rest of the same-different strain indicators via logistic regression. Data (bmbdf_file) should be within-host and include columns distance (SNP), timediff and colgroup.
    #'
    #' @param dfst The data frame we wish to add same strain indicator to (the consecutive distances data frame).
    #' @param bmbdf_file Path to the BaeMBac data created by run_baembac function. Used as the training data for logreg model.
    #' @return dfst Consecutive distances data frame with same strain indicator included.

    bmbdf <- read.csv(bmbdf_file)
    bmbdf$X <- NULL

    data <- data.frame(cbind(bmbdf$value, bmbdf$timediff, bmbdf$colgroup))
    names(data) <- c("distance", "timediff", "colgroup")
    data$colgroup <- as.factor(data$colgroup)

    dfst$colgroup <- rep(0, dim(dfst)[1])

    # Use logistic regression to assing same strain - different strain boolean. 1 = same strain, 0 = different strain.
    # Create the model and asses its accuracy:

    train <- data[0:round(dim(data)[1]/2),]
    start <- round(dim(data)[1]/2) + 1
    test <- data[start:dim(data)[1],]

    model <- glm(colgroup~., data = train, family=binomial)
    print(summary(model))

    # First, assess the model accuracy:

    fitted_colgroups <- predict(model, newdata = test, type='response')
    fitted_colgroups <- ifelse(fitted_colgroups > 0.5,1,0) # Decision boundary of 0.5. If above this, 1, if below, 0.
    #fitted_colgroups <- round(fitted_colgroups, digits = 0)

    #mean(fitted_colgroups == test$colgroup)

    print(paste("Model accuracy:", dim(test[which(test$colgroup == fitted_colgroups),])[1]/dim(test)[1]))

    # Then do the actual prediction:

    #test_st5 <- intersect(dfst5, bmbdf[which(bmbdf$ST == 5),])
    newdata <- dfst
    newdata <- data.frame(cbind(newdata$value, newdata$timediff, newdata$colgroup))
    names(newdata) <- c("distance", "timediff", "colgroup")

    fitted_colgroups <- predict(model, newdata = newdata, type='response')
    fitted_colgroups <- ifelse(fitted_colgroups > 0.5,1,0) # Decision boundary of 0.5. If above this, 1, if below, 0.

    dfst$colgroup <- fitted_colgroups

    return(dfst)
  
}
