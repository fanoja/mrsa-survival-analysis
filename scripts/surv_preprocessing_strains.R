# Formatting the data for survival analysis
# - Strain-based clearance

data_root <- "data/"
THRESHOLD <- 45 # from BaeMBac output figure z_vs_t_and_d_pred.pdf, 10 %, E arm

## Required scripts ##

source("scripts/utilities.R")
source("scripts/mdata_preprocessing.R")
source("scripts/strain_functions.R")

## Functions ##

# Helper functions for testing:

check_strain_data <- function(strain, mdata, bmbdf, surv_data){
    #' Plot a network figure of the host with strain 'strain'.
    #'
    #' In addition, print strain df, consecutive SNP distance data frame for the strain and restructured survival data.
    #'
    #' @param strain Strain id of interest.
    #' @param mdata Mykrobe data set with resistances etc metadata.
    #' @param bmbdf Consecutive SNP distances data frame with BaeMBac results included.
    #' @param surv_data Data in survival analysis compatible format (time difference, event indicator by intervals).
    
    print(mdata[which(mdata$substrain == strain),])
    print(surv_data[which(surv_data$strain == strain),])
    host <- unique(mdata[which(mdata$substrain == strain),"ID"])
    print(host)
    t <- vis_mrc_df_6(add_substrains(bmbdf[which(bmbdf$host1 == host & bmbdf$host2 == host),]))

}

# Utility functions:

compress_delta_t <- function(delta_t, min_or_max = max){
    #' Removes duplicate rows from delta_t. 
    #'
    #' If there are any duplicate Visit..s after this with different delta_t:s
    #'   corresponding to them, chooses the min or max aof these delta_t:s to assign to the visit according to the 
    #'   function specified in min_or_max argument.
    #'
    #' @param min_or_max Min or max function; used to choose the delta_t if the same visit has more than one different delta_t associated with it. By default max.
    #' @param delta_t An nx2 data frame with columns Visit.. and delta_t.
    #' @return delta_t Modified delta_t data frame.
    
    delta_t <- unique(delta_t)
    
    # More than one same Visit.. with different delta_ts after unique:
    max_visits <- names(table(delta_t$Visit..)[table(delta_t$Visit..) > 1])

    # For each of these max_visits, take only the max delta_t into account:
    for (v in max_visits){
        
        max_delta_t <- min_or_max(delta_t[which(delta_t$Visit.. == v),"delta_t"])

        delta_t <- delta_t[-which(delta_t$Visit.. == v),]

        delta_t <- rbind(delta_t, data.frame("Visit.." = v, "delta_t" = max_delta_t))

    }

    return(delta_t)
    
}



check_strain_switch_clearance <- function(clear, mdata, strain){
    #' Check for a clearance event associated with a strain switch.
    #'
    #' If strain A is followed by an MRSA isolate of strain B, this should be considered a clearance of strain A, not persistent colonization.
    #'
    #' @param clear A named vector with visit.Swab.site as names and 0, 1 and 6/NA as values.
    #' @param mdata Mykrobe data with strain identifiers.
    #' @return clear
    
    # Check if the non-cleared cases are of the same strain:
    not_cleared_cases <- names(clear[clear == "0"])
    
    # Extract visit and site from these (KEEP THE SAME ORDER!)
    extracted <- strsplit(not_cleared_cases, "\\.Swab\\.")
    for (ex in extracted){
        #print(mdata[which(mdata$Visit.. == ex[1] & mdata$Cx.Site == tolower(ex[2]) & mdata$substrain == strain),]) # If the not_cleared cases are not here, then they are of a different strain/ST.
        
        if (dim(mdata[which(mdata$Visit.. == ex[1] & mdata$Cx.Site == tolower(ex[2]) & mdata$strain == strain),])[1] == 0){
            
            clear[paste(ex[1], fupper(ex[2]), sep = ".Swab.")] <- 1 # Strain switch == clearance!
        }
    }
    
    return(clear)
}


remove_non_adh_intervals <- function(data, data_root = "data/"){
    #' Remove non-adherent observations from survival data.
    #'
    #' Note: excludes V4 from D (decolonization) arm, since the decolonization treatment was stopped at V3 -> hosts are no longer adherent in V4.
    #'
    #' @param data Survival data frame.
    #' @return surv_data Survival data after removing non-adherent observations from the decolonization arm.
    
    
    # Load adherence data
    adherence_data <- load_adherence_data(adherence_data_file)
    print(colnames(adherence_data))
    
    adherence_data[which(adherence_data$CHG.Bathing %in% c("Full", "HyperAdherent")),"CHG.Bathing"] <- 1
    adherence_data[which(adherence_data$CHG.Oral.Rinse %in% c("Full", "HyperAdherent")),"CHG.Oral.Rinse"] <- 1
    adherence_data[which(adherence_data$Mupirocin.Adherence %in% c("Full", "HyperAdherent")),"Mupirocin.Adherence"] <- 1
    adherence_data[which(adherence_data$CHG.Bathing %in% c("Partial", "None")),"CHG.Bathing"] <- 0
    adherence_data[which(adherence_data$CHG.Oral.Rinse %in% c("Partial", "None")),"CHG.Oral.Rinse"] <- 0
    adherence_data[which(adherence_data$Mupirocin.Adherence %in% c("Partial", "None")),"Mupirocin.Adherence"] <- 0
    
    colnames(adherence_data)[colnames(adherence_data) == "ID"] <- "host_id"

    # Merge survival data and adherence data
    surv_data <- merge(data, adherence_data[,c("Visit..", "host_id","CHG.Bathing", "CHG.Oral.Rinse", "Mupirocin.Adherence")], by = c("Visit..", "host_id"), all = T)
    surv_data <- surv_data[which(!is.na(surv_data$strain)),] # Remove NA cases
    orig_surv_data <- surv_data
    # Check: How many/are there D hosts do not have adherence data available:
    #print(paste("No adherence data for n_hosts=", dim(surv_data[which(is.na(surv_data$Mup.Adherence) & surv_data$ARM == "D"),])[1])) # all D hosts have adherence data available

    # Remove host-visit that is not adherent to all treatments at visit v1 (end of the interval):
    surv_data_D <- surv_data[which(surv_data$ARM == 1),]
    surv_data_E <- surv_data[which(surv_data$ARM == 0),]

    surv_data_D <- surv_data_D[which(surv_data_D$Mupirocin.Adherence == 1 & surv_data_D$CHG.Bathing == 1 & surv_data_D$CHG.Oral.Rinse == 1),]
    surv_data <- rbind(surv_data_D, surv_data_E)
    
    print(paste("Removed non-adherent cases, D arm n_removed_hosts=",length(unique(data[which(data$ARM == 1),"host_id"])) - length(unique(surv_data[which(surv_data$ARM == 1),]$host_id))))
    print(paste("Removed non-adherent cases, D arm n_removed_isolates=",dim(data[which(data$ARM == 1),])[1] - dim(surv_data[which(surv_data$ARM == 1),])[1]))
    
    print(paste("dim(surv_data) before removal:", dim(data)[1]))
    print(paste("n_hosts before removal:", length(unique(data$host_id))))
    print(paste("dim(surv_data) after removal:", dim(surv_data)[1]))
    print(paste("n_hosts after removal:", length(unique(surv_data$host_id))))
    
    return(surv_data)
}

add_missing_resistance <- function(mdata, debug = F, treatments = c('Ciprofloxacin','Clindamycin','Erythromycin','Gentamicin','Mupirocin','Rifampicin','Tetracycline','Trimethoprim','Chlorhexidine')){
    #' Add missing resistances based on the resistance of isolates of the same strain, on the same site.
    #'
    #' @param mdata Mykrobe data with strains added
    #' @param debug If TRUE, print some debug information
    #'
    #' @return mdata Modified mykrobe data
    
    #mdata[,treatments][!mdata[,treatments] %in% c("R", "S", "r")] <- "N"
    mdata[,treatments][mdata[,treatments] == "R"] <- 1
    mdata[,treatments][mdata[,treatments] == "S" | mdata[,treatments] == "r"] <- 0
    mdata[,treatments][mdata[,treatments] == "" | is.na(mdata[,treatments]) | mdata[,treatments] == "NA"] <- "N"

    #print(paste("Number of missing resistances:", length(mdata[,treatments][mdata[,treatments] == "N"])))
    print(paste("Adding resistance... Number of rows with missing resistances:", dim(mdata[unique(which(mdata[,treatments] == "N", arr.ind = T)[,1]),c("Isolate", "strain", treatments)])[1]))
    n_added_res <- 0


    strains <- unique(mdata$strain)
    #strains <- 142
    #strains <- c(245, 126, 13)
    
    for (strain in strains){ # assess on a strain-by-strain basis
    
        if(debug){print(strain)}
    
        for (treatment in treatments){

            if ("N" %in% mdata[which(mdata$strain == strain),treatment]){ # there is a missing observation for this antibiotic
                method <- paste("most common res status in strain;", treatment)
                tab <- table(mdata[which(mdata$strain == strain),treatment]) # counts of 0 and 1 and N per antibiotic per strain
                #tab <- tab[names(tab) != "N"] # remove N counts
                if(debug){print(tab)}

                if (length(tab) == 1){ # only N obs
                    if (names(tab) == "N"){ # only one isolate in the strain, and that is marked as NA. 
                        if(debug){
                           print("ONLY N")
                           print(tab) 
                        }
                        replace_NA_with <- "N" # unknown, if only one isolate/strain
                        method <- paste("uknown resistance (1 isolate only in the strain);", treatment)

                        if(debug){print(paste("Replace NA with, only N:", replace_NA_with))}
                    }

                }
                #print(names(tab)[tab == max(tab[names(tab) != "N"])]) # which has most obs, 0 or 1?
                else{
                    max_obs <- names(tab)[tab == max(tab[names(tab) != "N"])]
                    max_obs <- max_obs[max_obs != "N"] # don't count the number of "N" observations
                    if (length(max_obs) > 1){ # equal number of 0 and 1 cases
                        replace_NA_with <- "N" # RULE FOR ADDING MISSING DATA
                        method <- paste("unknown antibiotic resistance (equal 0 and 1 obs);", treatment)
                    }
                    else{
                        replace_NA_with <- max_obs
                    }
                }

                if(debug){
                    print(paste("Replaced NA with:", replace_NA_with))
                    print(paste("Method:", method))
                }

                #print(dim(mdata[which(mdata$strain == strain),treatment][mdata[which(mdata$strain == strain),treatment] == "N"])[1])
                n_added_res <- n_added_res + length(mdata[which(mdata$strain == strain),treatment][mdata[which(mdata$strain == strain),treatment] == "N"])[1]
                mdata[which(mdata$strain == strain),treatment][mdata[which(mdata$strain == strain),treatment] == "N"] <- replace_NA_with

                method <- paste("most common res status in strain;", treatment) # reset method desc to default   

            }
        }
    }

    print(paste("DONE. Number of rows with missing resistance after adding:",
               dim(mdata[unique(which(mdata[,treatments] == "N", arr.ind = T)[,1]),c("Isolate", "strain", treatments)])[1]))
    
    
    mdata[,treatments][mdata[,treatments] == "N"] <- NA # "N" to "official" NA
    
    mdata
}

# Helper function to create a dummy data frame:
get_empty_df <- function(colnames, dummy_value = 2){
    #' Get a 1xlength(colnames) data frame
    
    vec <- rep(2,length(colnames))
    names(vec) <- colnames
    df <- t(data.frame(vec))
    rownames(df) <- NULL
    
    as.data.frame(df)
}


# Main function to create the strain-based survival data
get_surv_data_strain <- function(mdata, 
                                 treatments = c('Ciprofloxacin','Clindamycin','Erythromycin','Gentamicin','Mupirocin','Rifampicin','Tetracycline','Trimethoprim','Chlorhexidine'),
                                 by_site = FALSE, site = "Nares"){
    #' Restructure the original data (mdata) into  a survival analysis compatible form.
    #'
    #' @param mdata Original dataset (Mykrobe data with strain ids).
    #' @param treatments Character vector of resistances that should be used as covariates.
    #' @param by_site Boolean indicating if the restructuring should be made by site (T) or over sites (F).
    #' @param site If by_site, construct survival data for this site.
    #' @return surv_data Survival analysis compatible data.
    
    # An empty surv_data data frame:
    surv_data <- data.frame(Visit.. = NA, strain = NA, host_id = NA, y = NA, ARM = NA, delta_t = NA, ST = NA)
    surv_data[,treatments] <- NA
    
    if (by_site){
        surv_data$site <- NA
        print("BY SITE!")
        mdata <- mdata[which(mdata$Cx.Site == tolower(site)),]
        print(table(mdata$Cx.Site))
    }
    
    
    print(paste("mdata before surv_data:", dim(mdata)[1]))
    print(paste("mdata before surv_data n_hosts:", length(unique(mdata$ID))))
    mdata <- mdata[which(!is.na(mdata$ID)),] # Exclude rows with NA as ID. Strain 228 will cause problems without this.
    print(paste("mdata, NA hosts removed:", dim(mdata)[1]))
    strains <- unique(mdata$strain)
    
    print(paste("mdata, overall strains:", length(unique(mdata$strain))))
    approx_delta_t <- c(30, 60, 90, 90)
    v1s <- c("V1", "V2", "V3", "V4")
    v0s <- c("R0", "V1", "V2", "V3")
    
    for (strain in strains){

        #print(paste("Strain:", strain))
        straindata <- mdata[which(mdata$strain == strain),]
        v0 <- straindata[which(straindata$Visit.. != "V4"),"Visit.."]
        
        v1 <- convert_v0_to_v1(v0) #v1s[match(v0, v0s)] # t+1 visits corresponding to the visits present at t in the strain. 
        #print(paste("Strain:", strain))
        
        if (length(v0 != 0)){ 
            if (by_site){
                v1_swab <- get_charvec_combo(c("R0", "V1", "V2", "V3", "V4"), fupper(site))
                y <- C[na.omit(as.character(unique(straindata$ID))),v1_swab] # 1 = cleared, 0 = not cleared, 6 = no swab
            }
            else{
                v1_swab <- get_charvec_combo(c("R0", "V1", "V2", "V3", "V4"), fupper(unique(mdata$Cx.Site))) # limit to sites present in mykrobe_data. This way, one can easily omit for example wound.
                y <- C[na.omit(as.character(unique(straindata$ID))), v1_swab]
            }
     
            y_prev <- y
                y <- check_strain_switch_clearance(y, mdata, strain) # Now this is clearance for this strain. Strain switch = clearance.
          
            if (length(y_prev) != length(y)){
                print("WARNING y != y before strain_switch_clearance")
                print(paste("Strain", strain))
                }
            
            names(y) <- substr(names(y), start = 1, stop = 2) # Rename the status vector y with visit names.
            
            clear <- unique(data.frame(names(y), rep(1, length(y)))) # By default, everything is cleared. This means that NAs will be treated as cleared.
            colnames(clear) <- c("Visit..", "y")

            for (v in clear$Visit..){
                

                if (0 %in% y[names(y) == v]){# Just one 0 per visit is enough to indicate a not cleared case (remember that 0 = not cleared)
                    clear[which(clear$Visit.. == v),"y"] <- 0 
                }
                if (6 %in% y[names(y) == v]){ # New 13.12.2021: Jump case fix.
                    clear[which(clear$Visit.. == v),"y"] <- 6 # Even one NA means that the y status is not reliably known
                }
                # In other cases, all 1s -> visit was cleared.
            }

            # Determine the "overall" resistance status 
            # (over all sites: just one resistant isolate during the visit - consider the strain resistant)

            res <- straindata[,c("Visit..",treatments)]

            res[res == "R"] <- 1
            res[res == "S"| res == "r"] <- 0

            #newres <- data.frame("Ciprofloxacin" = c(2), "Clindamycin" = c(2), "Erythromycin"=c(2), "Gentamicin"=c(2), 
                                 #"Mupirocin" = c(2), "Rifampicin" = c(2), "Tetracycline" = c(2), "Trimethoprim"=c(2), 
                                 #"Chlorhexidine"=c(2))
            newres <- get_empty_df(treatments) # 23.11.22 now generalizes to different resistance sets

            t_visits <- unique(res$Visit..)[unique(res$Visit..) != "V4"] # Resistance at t; V4 resistance is not relevant

            for (v in t_visits){
                vres <- res[which(res$Visit.. == v),treatments]
                #print(vres[,treatments])
                vres[,treatments] <- apply(vres[,treatments], 2, as.numeric)
                vres[,treatments] <- apply(vres[,treatments], 2, as.logical)
                #vres <- apply(na.omit(vres[,treatments]), 2, any)
                vres <- apply(vres[,treatments], 2, any) # new 3.11.2022; no omitting NAs
                newres <- rbind(newres, vres)
            }
            
            newres <- newres[-1,]
            newres <- cbind(t_visits, newres) # Assumes the same order for t_visits and res profiles.
            newres$Visit.. <- v1s[match(newres$t_visits,v0s)] # Find v1s corresponding to v0s
            newres <- merge(clear, newres, by = "Visit..") # Add clearance status y

            # Add the other relevant metadata: host_id, strain, ST and ARM:
            newres$ARM  <- rep(unique(straindata$ARM), dim(newres)[1])
            newres$host_id <- rep(unique(straindata$ID), dim(newres)[1])
            newres$strain <- rep(strain, dim(newres)[1])
            newres$ST <- rep(unique(straindata$ST), dim(newres)[1])
            
            if (by_site){
                newres$site <- rep(unique(straindata$Cx.Site), dim(newres)[1])
            }

            # Convert ARM D and E to 1 and 0:
            newres[which(newres$ARM == "D"),"ARM"] <- 1
            newres[which(newres$ARM == "E"), "ARM"] <- 0

            # Calculate delta_t (time between the consecutive visits; approximated or directly calculated)
            t0 <- straindata[which(straindata$Visit.. %in% v0s), c("Visit..", "Cx.Date")]
            t0$t1 <- v1
            # Find the dates corresponding to v1:
            date2 <- straindata[which(straindata$Visit.. %in% t0$t1), c("Cx.Date", "Visit..")]
            colnames(date2) <- c("date2","t1")
            delta_t <- merge(t0, date2, by = "t1", all = T) # Merge by v1 to get the available v1 dates.
            delta_t$delta_t <- as.numeric(difftime(as.Date(delta_t$date2, "%Y-%m-%d"), as.Date(delta_t$Cx.Date, "%Y-%m-%d"), units = "days"))
             # Timediff corresponding to the missing visits
            delta_t[which(is.na(delta_t$delta_t)), "delta_t"] <- approx_delta_t[match(delta_t[which(is.na(delta_t$delta_t)),"t1"], v1s)] # pick max values for delta_t if multiple different delta_ts for the same visit
            delta_t$Visit.. <- convert_v0_to_v1(delta_t$Visit..) #v1s[match(delta_t$Visit..,v0s)] # Convert Visit.. to v1.
            delta_t <- unique(delta_t[,c("Visit..", "delta_t")])

            # Add delta_t to newdata
            newres <- merge(compress_delta_t(delta_t), newres, by = "Visit..")

            # Remove excess columns:
            newres$t_visits <- NULL
            
            # Combine this resturcutred strain data (newres) with already structured data (surv_data)
            surv_data <- rbind(surv_data, newres)

        }
                
    }
    
    surv_data <- surv_data[-1,]
    
    return(surv_data)

}

## Functions for final preprocessing ##

modify_data <- function(data){
    #' Adds t0, turns delta_t to months, interval censoring fix (all event indicator 1s to 3s), omit NAs, omit V3-V4 observations.
    #'
    #' @param data Survival data frame.
    #' @return data Modified survival data frame.
    
    t0 <- as.matrix(rep(0.01, dim(data)[1]))
    colnames(t0) <- "t0"
    data <- cbind(data, t0)

    # Process delta_t by months
    data[,"delta_t"] <- data[,"delta_t"]/30

    # Interval censoring fix:
    data[which(data[,"y"] == 1),"y"] <- 3
    data[,"y"] <- as.numeric(data[,"y"])

    
    print(paste("surv_data before omitting missed swabs, n obs:", dim(data)[1]))
    print(paste("surv_data before omitting missed swabs, n hosts:", length(unique(data$host_id))))
    # For now, turn number 6 status values (missed swab) to NA and omit them:
    print(paste("Missing values:", length(data[which(data$y == 6), "y"])))
    data[which(data[,"y"] == 6),"y"] <- NA
    
    # Omit NA values:
    print(paste("Removing missing swabs, n =", dim(data[which(is.na(data[,"y"])),])[1]))
    data <- data[which(!is.na(data[,"y"])),]
    print(paste("surv_data after omitting missed swabs, n obs:", dim(data)[1]))
    print(paste("surv_data after omitting missed swabs, n hosts:", length(unique(data$host_id))))

    # Interval censoring with stan_surv:
    data[which(data[,"y"] == 0),"t0"] <- data[which(data[,"y"] == 0),"delta_t"]

    # Omit V4 observations: V3-V4 is not comparable with the first three intervals
    # Note that Visit.. indicates the end of the interval, thus remove V4 not V3.
    print(paste("Removing V3-V4 observations (t0 == V3), n =", dim(data[which(data[,"Visit.."] == "V4"),])[1]))
    data <- data[which(data[,"Visit.."] != "V4"),]
    print(paste("surv_data after omitting V4", dim(data)[1]))
    print(paste("surv_data after omitting V4, n hosts:", length(unique(data$host_id))))

    return(data)
}


## MAIN ## 

visit_data <- load_visit_data(visit_data_file)
mykrobe_data <- preprocess_mykrobe_data()

C <- get_clearance_data(visit_data)
C <- modify_clearmat_wound_6_to_1(C)

dfst5 <- readRDS(paste(data_root, "dfst5.RData", sep = ""))
dfst8 <- readRDS(paste(data_root, "dfst8.RData", sep = ""))

dfst <- add_threshold_strains(rbind(dfst5, dfst8), threshold = THRESHOLD)

# Add strain identifiers to mykrobe_data
mdata <- add_substrains_to_mykrobe_data(mykrobe_data, dfst)
mdata <- add_missing_resistance(mdata)  # add missing resistances if applicable

# Pooled over sites survival data
print("Constructing survival data...")
data <- modify_data(get_surv_data_strain(mdata))
data <- remove_non_adh_intervals(data) # Pooled sites


get_site_data <- function(mdata){
    #' Utility function to get site-specific data.
    #'
    #' @param mdata Mykrobe data frame with strains.
    #' @return data_sites Site-specific data.
    
    data_nares <- get_surv_data_strain(mdata, by_site = T, site = "Nares")
    data_throat <- get_surv_data_strain(mdata, by_site = T, site = "Throat")
    data_skin <- get_surv_data_strain(mdata, by_site = T, site = "Skin")
    data_wound <- get_surv_data_strain(mdata, by_site = T, site = "Wound")

    data_sites <- rbind(data_nares, data_throat, data_skin, data_wound)

    return(data_sites)
}

print("Constructing site-specific data...")
data_sites <- get_site_data(mdata)
data_sites <- remove_non_adh_intervals(modify_data(data_sites))

treatments <- c('Ciprofloxacin','Clindamycin','Erythromycin','Gentamicin','Mupirocin','Rifampicin','Tetracycline','Trimethoprim','Chlorhexidine')
treatments_td <- treatments[treatments != "Gentamicin"]


