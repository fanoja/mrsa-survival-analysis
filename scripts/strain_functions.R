# Utility functions for assigning strain ids etc to consecutive SNP distance data frames and mykrobe_data
source("scripts/utilities.R")

add_substrains <- function(df, to_add = 0){
    #' Add substrain id based on same strain strain indicator (stored in colgroup column, values are 1 if same strain and 0 if different strain)
    #'
    #' Creates a graph and assigns the same id to all connected isolates (same strain). Required packages: igraph.
    #'
    #' @param df The consecutive SNP distance data frame with same strain indicator included.
    #' @return df Consecutive SNP distance data frame with strain id.

    load_library("igraph")

    G <- graph_from_data_frame(df, directed = T)
    G <- G - E(G)[!E(G)$colgroup] # Remove nonexistent edges
    clstrs <- clusters(G) # Get the clusters

    ms <- clstrs$membership

    df$substrain1 <- rep(0, dim(df)[1])
    df$substrain2 <- rep(0, dim(df)[1])

    df$substrain1 <- ms[as.character(df$Var1)] + to_add #rep(0, dim(df)[1])
    df$substrain2 <- ms[as.character(df$Var2)] + to_add

    #df <- get_substrain_index(clstrs, df)
    print(paste("# of different substrains:", length(unique(c(df$substrain1, df$substrain2)))))

    return(df)
}


find_unique_isolates_within_ST <- function(mdata){
    #' Find hosts with only one isolate per ST and assign unique strain identifiers for these.
    #'
    #'
    #' @param mdata Mykrobe data frame.
    
    add_unique_strain <- c()
    
    for (st in c("5", "8")){

        tab <- table(mdata[which(mdata$ST == st),]$ID)
        tab <- tab[tab == 1] # Hosts with one isolate of the st of interest.

        isodf <- mdata[which(mdata$ST == st & mdata$ID %in% names(tab)),]
        #test_dims(dim(isodf)[1], length(unique(isodf$ID)), "Unique ID vs isolates") # each host has only one obs in this ST

        add_unique_strain <- c(add_unique_strain, isodf$Isolate) # these isolates should be unique within the host/ST group

    }
    
    return(unique(add_unique_strain))
 
}

add_substrains_to_mykrobe_data <- function(mykrobe_data, df){
    #' Add strain identifiers to Mykrobe data
    #'
    #' Limits Mykrobe data to ST5 and ST8, since these sequence types are clearly most prominent
    #'
    #' @param mykrobe data Mykrobe data frame.
    #' @param df Consecutive distances data frame with strain identifiers included.
    #' @return mdata_strains Mykrobe data with strain identifiers included. NAs assigned to isolates that do not have SNP distances available.
    
    mdata <- mykrobe_data[which(mykrobe_data$ST %in% c(5,8)),]
    
    max_strain <- max(c(df$substrain1, df$substrain2)) # current maximum strain id
    
    single_isolates <- find_unique_isolates_within_ST(mdata) # find single isolates within ST
    single_isolate_strains <- seq(max_strain + 1, max_strain + length(single_isolates), 1) # assign a unique strain id to each of these
    
    isostrain <- unique(data.frame("Isolate" = c(df$Var1, df$Var2, single_isolates), "strain" = c(df$substrain1, df$substrain2, single_isolate_strains))) # combine all available strain id information
    
    mdata_strains <- merge(mdata, isostrain, by = "Isolate", all = T) # add strain id to mykrobe_data
    
    return(mdata_strains)
   
}

## TEST FUNCTIONS ##

add_threshold_strains <- function(df, threshold=45, to_add = 0){
    #' For testing purposes, assign threshold-based strain ids to consecutive distances data frame df
    df$colgroup <- rep(0, dim(df)[1])
    df[which(df$value <= threshold), "colgroup"] <- 1 # replace this step with bmb results later
    df <- add_substrains(df, to_add)
    return(df)
}

test_single_isolate_st_host <- function(mykrobe_data, testhost){
    #' Check if an isolate marked as single occurence of the ST within the host is actually the single isolate.
    #'   If not, outputs a warning.
    #'   To check all hosts, iterate over mykrobe_data$ID
    
    mdata <- mykrobe_data[which(mykrobe_data$ST %in% c(5,8)),]
    add_unique_strain <- find_unique_isolates_within_ST(mdata)
    test <- mdata[which(mdata$ID %in% mdata[which(mdata$Isolate %in% add_unique_strain), "ID"]),]
    
    
    if (testhost %in% test$ID){ # otherwise, the host did not have a single isolate/ST.
        #test[which(test$ID == testhost),]
        tab <- table(test[which(test$ID == testhost),"ST"])

        candidate_single_st <- as.numeric(test[which(test$ID == testhost & test$Isolate %in% add_unique_strain),"ST"])
        #print(tab)
        #print(candidate_single_st)
        
        for (st in candidate_single_st){ # sometimes, the host has 1 isolate of ST5 and 1 isolate of ST8
            if(unname(tab[names(tab) == st]) != 1){ # Check that this ST is only there once
                warning("FAILED. There is a problem with single isolate of ST.")
                print(paste("Check host:", testhost))}
        }
    }
    
    
}

# Check ALL of the single isolate hosts. If this outputs nothing, all good. SLOW.
#print("Testing...")
#for (host in mykrobe_data$ID){ # Warning: Takes some time.
#    test_single_isolate_st_host(mykrobe_data, host)
#}
#print("Done.")

