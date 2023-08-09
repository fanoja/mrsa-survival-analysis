## PREPROCESSING ##

# Required scripts: utilities.R
source("scripts/utilities.R")
#source("scripts/utilities.R")


# Location of raw data files
data_root <- "/u/50/ojalaf2/unix/Dropbox (Aalto)/Mrsa_clear/jupyter_notebooks/Publish_saved_models_and_data/"

# Raw data file names
mykrobe_data_file <- "PCLEAR_Mykrobe_R_predictions.tsv"
blast_data_file <- "PCLEAR_genotype_prediction_and_blast_results.tsv"
visit_data_file <- "Visit_and_Swab_Records.csv"
adherence_data_file <- "Adherence.csv"

output <- T # If T, will print the number of isolates/hosts after certain preprocessing steps
include_adherence <- F # If T, adds adherence data to mykrobe_data
call_preprocess_mykrobe_data <- T # If T, this script calls preprocess_mykrobe_data function to perform all preprocessing steps



## Loading relevant datasets: mykrobe data, visit data, adherence data

load_visit_data <- function(visit_data_file){
    #' Load visit data.
    #'
    #' Visit data contains information on completed visits by host and MRSA positive swabs by body site.
    #'
    #' @param visit_data_file The name of the file containing raw visit data.
    #' @return visit_data Visit data as a data frame.
    
    visit_data <- read.csv(paste(data_root, visit_data_file, sep = ""), sep = ";")
    
    return(visit_data)
}

get_clearance_data <- function(vdata){
    #' Construct a clearance matrix using mykrobe data and visit data.
    #'
    #' Construct a clearance matrix with actual clearances, not cleared and not swabbed cases.
    #' @param vdata A data set with completion indicators for swabs and visits by host.
    #' @return C A 5xn_host clearance matrix with MRSA positive visits denoted by 0, cleared visits by 1 and missing observations by 6
    #'   Each row denotes one host.
    
    print("Constructing clearance matrix...")
    visit <- c("R0", "V1", "V2", "V3", "V4")
    sites <- c("Nares", "Throat", "Skin", "Wound")

    # Reconstruct the column names from visit_data
    swab_colnames <- c()
    MRSA_colnames <- c()
    for (site in c("Nares", "Throat", "Skin", "Wound")){
        swab_colnames <- c(swab_colnames, paste(visit, paste("Swab", site, sep = "."), sep = "."))
        MRSA_colnames <- c(MRSA_colnames, paste(visit, paste("MRSA.", site, sep = "."), sep = "."))
    }
    
    swabs <- vdata[,swab_colnames]
    mrsa <- vdata[,MRSA_colnames]

    # Trues are true clearances. Falses are either skips or true positives.
    C <- swabs == 1 & mrsa == 0
    
    C[C == TRUE] <- "cleared"
    

    # Mark missing swabs:
    MS <- swabs == 0 & mrsa == 0
    C[MS == TRUE & C == FALSE] <- "missing_swab"
    table(C)

    # Remaining Falses are true positives. Make sure it is the case:
    TP <- swabs == 1 & mrsa == 1
    C[TP == TRUE & C == FALSE] <- "mrsa_pos"
    table(C) # Here, no FALSEs should remain.
    head(C)

    #  cleared == 1, missing_swab == 6 (to avoid stan_surv confusion), mrsa_pos == 0.
    C[C == "missing_swab"] <- 6
    C[C == "cleared"] <- 1
    C[C == "mrsa_pos"] <- 0

    rownames(C) <- vdata$ID # Rows by host
    
    return(C)
}

get_charvec_combo <- function(vec1, vec2, sep=".Swab."){
    #' Combine the contents of the character vectors vec1 and vec2 into one string separated by sep.
    #'
    #' Takes into account all combinations. In our case, vec1 = visit, vec2 = site, sep = ".Swab."
    #'
    #' @param vec1 First character vector to combine.
    #' @param vec2 Second character vector to combine.
    #' @return comb All combinations of vec1 and vec2 separated by sep.
    
    comb <- c()

    for (s1 in vec1){
        comb <- c(comb, paste(s1, vec2, sep = sep))
    }
    
    return(comb)
}

modify_clearmat_wound_6_to_1 <- function(C){
    #' Modify values for wound in clearance matrix C.
    #'
    #' Transform 6s of wound (missing observations) to 1 (clearance), since NAs in wound indicate that
    #'   the wound was not found -> wound was healed -> no MRSA in wound.
    #'   Tested by comparing with table() the distributions of 1, 0 and 6 with modified and unmodified C.
    #'
    #' @param C The clearance matrix.
    #' @return clearmat The modified clearance matrix
    
    wound_cols <- get_charvec_combo(c("R0", "V1", "V2", "V3", "V4"), fupper("wound")) # In wound, NA == no wound == clearance
    clearmat <- as.data.frame(C)
    
    woundmat <- clearmat[,wound_cols]
    woundmat[woundmat == 6] <- 1
    clearmat <- clearmat[!(colnames(clearmat) %in% wound_cols)]

    clearmat <- cbind(clearmat, woundmat)
    clearmat <- as.matrix(clearmat)
    
    return(clearmat)
   
}

load_adherence_data <- function(adherence_data_file){
    #' Load and preprocess adherence data.
    #'
    #' Adherence data contains information on the adherence to the decolonization protocol (mupirocin and chlorhexidine) by host.
    #'   Note: no adherence data available for R0 or V4, because the hosts are no longer treated at V4 and R0 is enrollment.
    #'
    #' @param adherence_data_file Name of the raw adherence data file.
    #' @return adherence_data Loaded and preprocessed adherence data as a data frame.
    
    adherence_data <- read.csv(paste(data_root, adherence_data_file, sep = ""), sep = ";")
    adherence_data <- plyr::rename(adherence_data, c("Study.ID" = "ARM", "X"="ID", "X.1" = "Unknown")) # rename some columns
    
    adherence_data <- read.csv(paste(data_root, "Adherence.csv", sep = ""), sep = ";")
    adherence_data[which(adherence_data$Visit == "1-month"),"Visit"] <- "V1"
    adherence_data[which(adherence_data$Visit == "3-month"),"Visit"] <- "V2"
    adherence_data[which(adherence_data$Visit == "6-month"),"Visit"] <- "V3"
    adherence_data <- plyr::rename(adherence_data, c("Study.ID" = "ARM", "X"="ID", "X.1" = "group", "Visit" = "Visit..", "Mupirocin" = "Mupirocin.Adherence"))
    
    return(adherence_data)
}

load_blast_data <- function(blast_data_file){
    #' Load and preprocess blast_data
    #'
    #' BLAST data contains the resistance indicators for mupirocin and chlorhexidine as obtained using BLAST in addition to the metadata. It also includes additional data provided by BLAST (such as presence of mupA gene etc).
    #'
    #' @param blast_data_file Name of the data file containing raw BLAST data and metadata.
    #' @return blast_data Loaded and preprocessed BLAST data.
    
    blast_data <- unique(read.csv(paste(data_root, blast_data_file, sep= ""), sep = "\t"))
    blast_data <- blast_data[which(!blast_data$Isolate %in% c("5024_45")),] # this isolate id has twin observations
    blast_data <- plyr::rename(blast_data, c("prediction_Chlorhexidine" = "Chlorhexidine"))
    
    return(blast_data)
}

load_mykrobe_data <- function(mykrobe_data_file, output = F){
    #' Load and preprocess mykrobe data.
    #'
    #' Mykrobe data contains resistance predictions for different antibiotics as predicted by the Mykrobe predictor software, in addition to metadata. The metadata in mykrobe data is the same as in BLAST data.
    #'
    #' @param mykrobe_data_file File name of mykrobe data.
    #' @return mykrobe_data Preprocessed mykrobe data frame.
    
    mykrobe_data <- read.csv(paste(data_root, mykrobe_data_file, sep = ""), sep = "\t") # Load the original data set with Mykrobe resistance predictions
    mykrobe_data$Cx.Site <- tolower(mykrobe_data$Cx.Site) # body site to lowercase
    mykrobe_data[which(mykrobe_data$Isolate == "5045_58"), "Cx.Date"] <- "6/19/13" # One isolate has clearly the wrong year associated with them. (Host 1510, V2 has a smaller date than R0)
    mykrobe_data  <- mykrobe_data[which(mykrobe_data$ID != 3058),] # This host is assigned to both D and E arms simultaneously
    
    # Fix broken dates: All hosts expect these had the same date associated with a certain visit for all bodysites.
    #  These are assumed to be typos and will be fixed accordingly.
    mykrobe_data[which(mykrobe_data$ID == 1077 & mykrobe_data$Cx.Site == "wound" & mykrobe_data$Visit.. == "R0"), "Cx.Date"] <-  "4/22/11" # Choose nares swab date (although, impossible to know which date is the correct one)
    mykrobe_data[which(mykrobe_data$ID == 1703 & mykrobe_data$Cx.Site == "skin" & mykrobe_data$Visit.. == "V4"), "Cx.Date"] <- "1/28/12" # This skin isolate seems like it should be at V3 and a year earlier.
    mykrobe_data[which(mykrobe_data$ID == 1703 & mykrobe_data$Cx.Site == "skin" & mykrobe_data$Visit.. == "V4"), "Visit.."] <-  "V3"
    # These three have a difference of 1 day. Choosing the nares swab as the actual swab (no effect on the results, since the difference is small)
    mykrobe_data[which(mykrobe_data$ID == 1963 & mykrobe_data$Cx.Site == "wound" & mykrobe_data$Visit.. == "V3"), "Cx.Date"] <- "8/9/13" 
    mykrobe_data[which(mykrobe_data$ID == 6199 & mykrobe_data$Cx.Site == "skin" & mykrobe_data$Visit.. == "R0"), "Cx.Date"] <- "7/30/13"
    mykrobe_data[which(mykrobe_data$ID == 2592 & mykrobe_data$Cx.Site == "wound" & mykrobe_data$Visit.. == "R0"), "Cx.Date"] <- "8/20/13"

    mykrobe_data$Cx.Date <- as.Date(mykrobe_data$Cx.Date, "%m/%d/%Y") # date to Date object
    
    mykrobe_data[which(mykrobe_data$Isolate == "-"),]$Isolate <- mykrobe_data[which(mykrobe_data$Isolate == "-"),]$Accession.. # Some isolates lack a name - replace with Accession id
    # Replace all dashes with underscores in mykrobe data (because they are marked as such in the distance matrices)
    mykrobe_data$Isolate <- gsub("-", "_", mykrobe_data$Isolate) # There should be 33 of these
    mykrobe_data[which(mykrobe_data$Isolate %in% c('T12805','W24219')),"Isolate"] <- mykrobe_data[which(mykrobe_data$Isolate %in% c('T12805','W24219')),"Accession.."] # Two isolate ids are identical, replace their names with the  corresponding accession id

    if (output){
        print(paste("Isolates in mykrobe_data (unfiltered):", dim(mykrobe_data)[1]))
        print(paste("Hosts in mykrobe_data (unfiltered):", length(unique(mykrobe_data$ID))))      
    }
    
    return(mykrobe_data)
    
}

## ADDITIONAL MYKROBE DATA PREPROCESSING ##

add_blast_results_to_mykrobe_data <- function(mykrobe_data, blast_data_file, output = F){
   #' Add ALL BLAST results to mykrobe_data by a merge.
    #'
    #' Since Mykrobe predictor does not include resistance prediction for chlorhexidine, add this to mykrobe data from BLAST data.
    #'
    #' @param mykrobe_data Loaded mykrobe_data as returned by the `load_mykrobe_data` function.
    #' @param blast_data_file Name of the file containing BLAST data.
    #' @return mykrobe_data Mykrobe data with chlorhexidine resistance included.
    
    print("Adding Chlorhexidine from BLAST...")
    
    blast_data <- load_blast_data(blast_data_file)
    
    mykrobe_data_blast <- merge(mykrobe_data, blast_data, by = "Isolate", all = T)
    test_dims(dim(mykrobe_data_blast)[1] - dim(mykrobe_data)[1], length(setdiff(blast_data$Isolate, mykrobe_data$Isolate)), test_name = "Mykrobe and BLAST merge")
    
    mykrobe_data <- mykrobe_data_blast[which(!is.na(mykrobe_data_blast$ID)),]

    if (output){print(paste("nrow after adding BLAST:", dim(mykrobe_data)[1]))} 
    
    return(mykrobe_data)
}

## Optional: Add adherence data to mykrobe data ##

add_adherence_info_to_mykrobe_data <- function(mykrobe_data, adherence_data_file, output = F){
    #' Combine adherence data and mykrobe data. Optional.
    #'
    #' @param mykrobe_data Mykrobe data frame.
    #' @param adherence_data_file Name of the adherence data file.
    #' @return mykrobe_data_adh Mykrobe data with adherence included.
    
    print("Adding adherence information...")
    
    adherence_data <- load_adherence_data(adherence_data_file)
    mykrobe_data_adh <- merge(adherence_data[,c("ARM", "ID", "Visit..", "CHG.Bathing", "CHG.Oral.Rinse", "Mupirocin.Adherence")], mykrobe_data, by = c("ID", "Visit..", "ARM"), all = T)
    
    return(mykrobe_data_adh)
    
}

remove_incompleted_hosts <- function(mykrobe_data, visit_data_file, output = F){
    #' Remove hosts that did not complete all visits from Mykrobe data
    #'
    #' @param mykrobe_data Mykrobe data frame
    #' @param visit_data_file Name of the visit data file
    #' @param output If true, prints additional output from this function for debugging.
    #' @return mykrobe_data Mykrobe data without hosts that did not complete all visits.
    
    print("Adding visit data...")
    
    visit_data <- load_visit_data(visit_data_file)
    
    # Find the hosts that completed all visits: 
    all_visits_completed_hosts <- unique(visit_data[which(visit_data$FINAL.2121 == 1 & visit_data$R0.Visit.Completed == 1 & visit_data$V1.Visit.Completed == 1 &visit_data$V2.Visit.Completed == 1 & visit_data$V3.Visit.Completed == 1 & visit_data$V4.Visit.Completed == 1),"ID"])
    sus_hosts <- setdiff(mykrobe_data$ID, all_visits_completed_hosts)

    if(output){    
        print(paste("Removing hosts that did not complete all visits, n_hosts removed =",length(unique(mykrobe_data$ID)) - length(unique(intersect(mykrobe_data$ID, all_visits_completed_hosts)))))
        print(paste("Removing hosts that did not complete all visits, n_isolates removed = ", dim(mykrobe_data[which(mykrobe_data$ID %in% sus_hosts),])[1]))
    }

    # Filter out hosts that did not complete all visits:
    mykrobe_data_visits <- mykrobe_data[which(mykrobe_data$ID %in% all_visits_completed_hosts),] # This row excludes NA hosts also
    mykrobe_data <- mykrobe_data_visits
    test_dims(length(unique(mykrobe_data$ID)), length(intersect(all_visits_completed_hosts, mykrobe_data$ID)), test_name = "Removed hosts that did not complete all visits")
    if (output){print(paste("nrow after removing hosts that did not attend all visits:", dim(mykrobe_data)[1]))}    
    
    return(mykrobe_data)
    
 }


## MISC PREPROCESSING ##

remove_contaminated_isolates <- function(mykrobe_data, contaminated = c("5007_56", "5025_02", "5016_39", "H30734-6", "H30734-7", "5061_37", "5009_01", "5012_44", "5025_20", "5009_53", "5050_65", "5016_47", "5016_40"), output = F){
    #' Remove contaminated isolates from mykrobe_data
    #'
    #' Some isolates were known to be contaminated. These are excluded from the data.
    #'
    #' @param mykrobe_data Mykrobe data frame.
    #' @param contaminated A vector containing the isolate ids of the contaminated isolates.
    #' @param output If true, prints additional output for debugging.
    #' @return mykrobe_data Mykrobe data without contaminated isolates.
    
    print("Removing contaminated isolates..")
    
    mykrobe_data_cont <- mykrobe_data[which(!mykrobe_data$Isolate %in% contaminated),]
    test_dims(length(intersect(contaminated, mykrobe_data$Isolate)), dim(mykrobe_data)[1] - dim(mykrobe_data_cont)[1])

    if(output){
        print(paste("Removed contaminated isolates, n_isolates removed:", dim(mykrobe_data)[1] - dim(mykrobe_data_cont)[1]))
        print(paste("Removed contaminated isolates, n_hosts removed:", length(unique(mykrobe_data$ID)) - length(unique(mykrobe_data_cont$ID))))
    }

    mykrobe_data <- mykrobe_data_cont
    
}


# Remove Methicillin == S cases

remove_MSSA <- function(mykrobe_data, output = F){
    #' Remove isolates that are marked as MSSA in Mykrobe data.
    #'
    #' Note that the MSSA cases are most likely due to error related to Mykrobe predictor predictions.
    #'
    #' @param mykrobe_data Mykrobe data frame.
    #' @param output If true, print additional output for debugging.
    #' @return mykrobe_data_mssa Mykrobe data without isolates marked as MSSA.
    
    print("Removing MSSA cases...")
    mssa_isolates <- mykrobe_data[which(mykrobe_data$Methicillin == "S"),"Isolate"]
    mykrobe_data_mssa <- mykrobe_data[which(!mykrobe_data$Isolate %in% mssa_isolates),]
    test_dims(length(unique(mssa_isolates)), length(unique(mykrobe_data$Isolate)) - length(unique(mykrobe_data[which(!mykrobe_data$Isolate %in% mssa_isolates),"Isolate"])), test_name = "MSSA isolates")

    if (output){
        print(paste("Removed MSSA cases, n_isolates removed:", dim(mykrobe_data)[1] - dim(mykrobe_data_mssa)[1]))
        print(paste("Removed MSSA cases, n_hosts removed:", length(unique(mykrobe_data$ID)) - length(unique(mykrobe_data_mssa$ID))))
    }

    return(mykrobe_data_mssa)
}

remove_duplicate_swabs <- function(mykrobe_data, duplicates_to_remove = c("5016_56", "5052_05", "5052_09", "5007_15", "5011_25", "5016_54","T24515", "5036_51", "H30734_7"), output = F){
    #' Remove duplicate swabs
    #'
    #' Some body sites were swabbed more than one time and recorded more than once. Remove the latter observation of these duplicates.
    #'   For example, host 1724 has 2 isolates from nares at V1 - remove the latter of these.
    #'
    #' @param mykrobe_data Mykrobe data frame.
    #' @param duplicates_to_remove A vector of duplicate isolate ids to remove. These were observed manually from the data.
    #' @param output If true, prints additional output.
    #' @return mykrobe_data_dupl Mykrobe data without duplicate swabs.
    
     print("Removing duplicates: same host, same visit, same site.")
    duplicates_to_remove <- c("5016_56", "5052_05", "5052_09", "5007_15", "5011_25", "5016_54","T24515", "5036_51", "H30734_7") #
    mykrobe_data_dupl <- mykrobe_data[which(!mykrobe_data$Isolate %in% duplicates_to_remove),]

    if (output){
        print(paste("Removed duplicates, n_isolates=", length(mykrobe_data$Isolate) - length(mykrobe_data_dupl$Isolate), sep = ""))
    }

    return(mykrobe_data_dupl)
}


### MAIN FUNCTION: call this to perform all preprocessing steps on mykrobe_data ###

preprocess_mykrobe_data <- function(output = F){
    #' Perform all preprocessing steps on Mykrobe data (see above functions).
    #'
    #' @param output True for additional input, by default False.
    #' @return mykrobe_data Preprocessed Mykrobe data.
    
    mykrobe_data <- load_mykrobe_data(mykrobe_data_file, output = output) # unpreprocessed data
    mykrobe_data <- add_blast_results_to_mykrobe_data(mykrobe_data, blast_data_file, output = output) # add chlorhexidine from BLAST data
    if (include_adherence){ # add adherence info, optional
        mykrobe_data <- add_adherence_info_to_mykrobe_data(mykrobe_data, adherence_data_file, output = output)
    }
    mykrobe_data <- remove_incompleted_hosts(mykrobe_data, visit_data_file, output = output) # remove hosts that did not complete all visits
    mykrobe_data <- remove_contaminated_isolates(mykrobe_data, output = output) # remove contaminated isolates
    mykrobe_data <- remove_MSSA(mykrobe_data, output = output) # remove isolates marked as mssa
    mykrobe_data <- remove_duplicate_swabs(mykrobe_data, output = output) # remove duplicate swabs (same host, same site, same visit)
    
    if(output){
        print(paste("Hosts in mykrobe_data (all_visits_completed, contaminated removed, methicillin != S, no duplicates):", length(unique(mykrobe_data$ID))))
        print(paste("Isolates in mykrobe_data (all visits completed, contaminated removed, methicillin != S, no duplicates):", dim(mykrobe_data)[1]))
    }
    
    return(mykrobe_data)
    
}


if (call_preprocess_mykrobe_data){
    mykrobe_data <- preprocess_mykrobe_data(output = output)
}