# Functions for creating consecutive distance data frames.
# Test functions included for assessing input. Set debug = T to assess.

debug <- F # print out test results for debugging

# Load mykrobe data
#source("mdata_preprocessing.R")

get_all_cons_combos <- function(sites, visits, hosts, add.eos = F){
    #' Get all possible combinations of consecutive/simultaneous visits for each site1-site2 pair. 
    #'
    #' Returns a data frame of these combinations for each host.
    #'
    #' @param sites A vector containing unique body sites in the data (here, nares, skin, throat, wound).
    #' @param visits A vector containing visit names (R0, V1, etc..).
    #' @param hosts A vector containing host ids of interest.
    #' @param add.eos If true, include an end-of-study indicator after V4.
    
    if (add.eos){
        visits <- c(visits, "EOS")
    }
    
    cons_sites <- expand.grid(sites, sites, stringsAsFactors = FALSE)
    cons_visits <- expand.grid(visits, visits, stringsAsFactors = FALSE)
    cons_visits <- cons_visits[which(cons_visits$Var1 <= cons_visits$Var2),] # for example, V1->R0 is not allowed.
    
    if (add.eos){
        cons_visits <- cons_visits[which(cons_visits$Var1 != "EOS"),]
    }
    
    n_hosts <- length(hosts)
    n_visit_combos <- dim(cons_visits)[1]
    n_site_combos <- dim(cons_sites)[1] # number of site1-site2 combinations

    df <- data.frame(v0 = rep(cons_visits$Var1, n_site_combos*n_hosts),
                v1 = rep(cons_visits$Var2, n_site_combos*n_hosts),
                site0 = rep(cons_sites$Var1, n_visit_combos*n_hosts),
                site1 = rep(cons_sites$Var2, n_visit_combos*n_hosts))

    df$host <- rep(hosts, each = n_site_combos*n_visit_combos)
    
    return(df)
    
}

test_get_all_cons_combos <- function(sites, visits, hosts){
    #' Test function for get_all_cons_combos.
    #'
    #' @param sites A vector of sites of interest
    #' @param visits A vector of visits of interest
    #' @param hosts A vector of hosts of interest
    
    print(paste("Sites:", toString(sites)))
    print(paste("Visits:", toString(visits)))
    print(paste("Hosts:", toString(hosts)))
    print(get_all_cons_combos(sites,visits,hosts))
 
}

if (debug){
    test_get_all_cons_combos(c("nares"), c("R0"), c(100))
    test_get_all_cons_combos(c("nares", "throat"), c("R0"), c(100))
    test_get_all_cons_combos(c("nares"),  c("R0"), c(100, 23))
    test_get_all_cons_combos(c("nares", "throat"), c("R0", "V1"), c(100))
}

load_SNP_distances <- function(st_dist_file){
    #' Load and preprocess SNP distance matrix from file st_dist_file.
    #'
    #' @param st_dist_file File name of the SNP distance matrix of interest. There are separate files for ST5 and ST8 isolates.
    #' @return st The SNP distance matrix as data frame.
    
    st <- read.csv(st_dist_file, sep = "\t")
    
    # The first row is the taxon. Remove it.
    st <- st[,-1] # ASSUMPTION: colnames is in the same order as rownames.
    
    # Remove the "X" from the names of isolates.
    colnames(st) <- gsub("X","",colnames(st))
    
    # Some values in colnames of st5 and st8 have a '.' instead of '_'. Fix by replacing.
    colnames(st) <- gsub("\\.", "_", colnames(st))
    
    # Set row names.
    rownames(st) <- colnames(st)
    
    return(st)
    
}

get_most_recent_distances_df <- function(mykrobe_data, st, st_dist_file, add.eos = F){
    #' Get consecutive SNP distances within-host.
    #'
    #' Find SNP distances between two isolates that were sampled at consecutive visits.
    #'
    #' @param mykrobe_data Preprocessed Mykrobe data frame.
    #' @param st ST of interest, either 5 or 8.
    #' @param st_dist_file File name of the SNP distance matrix file.
    #' @param add.eos If true, add end-of-study indicator.
    #' @return df Data frame with SNP distances between consecutive visits within-host.
    
    # Load the SNP distance matrix:
    st_distmat <- load_SNP_distances(st_dist_file)
    
    visits <- c("R0", "V1", "V2", "V3", "V4")
    
    sites <- c("nares", "throat", "skin", "wound")

    mdata <- mykrobe_data[which(mykrobe_data$ST == st & mykrobe_data$Isolate %in% c(colnames(st_distmat))),] # Limit to known SNP distances.

    # Find all site1-site2, v0-v1 and host combinations
    df <- get_all_cons_combos(sites, visits, unique(mdata$ID), add.eos = add.eos)

    # Merge v1
    df <- merge(df, mdata[,c("Isolate", "Cx.Site", "ID", "Visit..")],
                    by.x = c("v1", "site1", "host"),
                    by.y = c("Visit..", "Cx.Site", "ID"))
    colnames(df)[colnames(df) == "Isolate"] <- "Var2"

    # Merge v0
    df <- merge(df, mdata[,c("Isolate", "Cx.Site", "ID", "Visit..", "ARM")],
                    by.x = c("v0", "site0", "host"),
                    by.y = c("Visit..", "Cx.Site", "ID"))
    colnames(df)[colnames(df) == "Isolate"] <- "Var1"

    ## SNP distances, ARM and timediff
    
    # Add snp distances and the time difference in days

    df$value <- rep(NA, dim(df)[1])
    df$timediff <- rep(NA, dim(df)[1])

    for (ri in 1:dim(df)[1]){

        if (df$Var1[ri] %in% colnames(st_distmat) & df$Var2[ri] %in% colnames(st_distmat)){ # check that the SNP dist is not missing
            df[which(df$Var1 == df$Var1[ri] & df$Var2 == df$Var2[ri]), "value"] <- st_distmat[df$Var1[ri], df$Var2[ri]]
            df[which(df$Var1 == df$Var1[ri] & df$Var2 == df$Var2[ri]), "timediff"] <- difftime(mdata[which(mdata$Isolate == df$Var2[ri]), "Cx.Date"],
                                                                                               mdata[which(mdata$Isolate == df$Var1[ri]), "Cx.Date"],
                                                                                               units = "days")
        }
    }

    # Final additions: ST
    df$ST <- rep(st, dim(df)[1])
   
    # Reorder the columns of df (essential for assigning strain ids)
    col_order <- c("Var1", "Var2", "value", "timediff", "host", "v0", "site0", "v1", "site1", "ST", "ARM")
    df <- df[,col_order]
    
    # Remove duplicates: same visit, different site. 
    # - For example skin-throat == throat-skin at V4-V4. 
    to_keep <- paste(c("nares", "nares", "nares", "skin", "skin", "throat"),c("skin", "throat", "wound", "throat", "wound", "wound"), sep = "_") # These duplicates will be removed
    df <- df[-which(paste(df$site0, df$site1, sep = "_") %in% to_keep & df$v0 == df$v1),]

    if(length(df[which(df$timediff < 0),"timediff"]) > 0){
        warning("Negative timediff. Multiply with -1")
        print(df[which(df$timediff < 0),])
        df[which(df$timediff < 0), "timediff"] <- -1*df[which(df$timediff < 0), "timediff"]
    }
    
    return(df)
    
}

test_random_cons_isolates <- function(df, mdata, st_distmat){
    #' Test a pair of random isolates: compare consecutive distances data frame and mykrobe_data
    #'
    #' @param df Consecutive SNP distances data frame.
    #' @param mdata Mykrobe data frame
    #' @param st_distmat Preprocessed SNP distance matrix.
    
    iso1 <- sample(df$Var1, 1)
    iso2 <- sample(df[which(df$Var1 == iso1), "Var2"],1)
    print(paste("Testing isolates:", toString(c(iso1, iso2))))
    print("Mykrobe data:")
    print(mdata[which(mdata$Isolate %in% c(iso1, iso2)),])
    print(difftime(mdata[which(mdata$Isolate == iso2), "Cx.Date"], mdata[which(mdata$Isolate == iso1), "Cx.Date"], units = "days"))
    print("Df (consecutive distances):")
    print(df[which(df$Var1 == iso1 & df$Var2 == iso2),])
    print("SNP distance matrix:")
    print(st_distmat[iso1, iso2])
}

if(debug){
    print("Printing a pair of random consecutive isolates:")
    test_random_cons_isolates(df, mdata, st_distmat)
}


filter_dfst <- function(dfst){
    #' Filter out distances between the isolate and itself.
    #'
    #' Note: will remove one isolate hosts completely! These are added back later on as isolates of a unique strain.
    #'   Different hosts always have different strains according to our definition.
    
    dfst_filtered <- dfst[which(dfst$Var1 != dfst$Var2),]

    return(dfst_filtered)
    
}

