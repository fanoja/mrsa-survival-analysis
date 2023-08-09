## Preprocessing of the raw phenotypic data file

debug <- FALSE

# Load original data before changing to surival data format:

source("scripts/mdata_preprocessing.R")
source("scripts/strain_functions.R")

dfst5 <- readRDS("data/dfst5.RData")
dfst8 <- readRDS("data/dfst8.RData")

dfst <- add_threshold_strains(rbind(dfst5, dfst8))

mdata <- add_substrains_to_mykrobe_data(mykrobe_data, dfst)

# Load AST data
install.packages("readxl")
library("readxl")
ast_data_uci <- read_excel("/u/50/ojalaf2/unix/Dropbox (Aalto)/Mrsa_clear/jupyter_notebooks/Spring2022/data/PCLEAR_Raw_AST_Phenotypes_modified.xlsx", sheet = "Site UCI") # Loads only UCI data.
ast_data_ucla <-  read_excel("/u/50/ojalaf2/unix/Dropbox (Aalto)/Mrsa_clear/jupyter_notebooks/Spring2022/data/PCLEAR_Raw_AST_Phenotypes_modified.xlsx", sheet = "Site UCLA")

# AST data preprocessing:

ast_data <- rbind(as.data.frame(ast_data_uci), as.data.frame(ast_data_ucla))
ast_data <- plyr::rename(ast_data, c("Accession #" = "Accession..", "Culture Site" = "Cx.Site", "Visit" = "Visit.."))
ast_data[which(ast_data$Visit.. == "RO"),"Visit.."] <- "R0"
print(dim(ast_data))

# Remove excess columns:
ast_data <- ast_data[,!(names(ast_data) %in% c("MUP MIC (ug/ml)", "CHG MIC (ug/ml)", "QAC A/B PCR RESULT"))]
ast_data <- unique(ast_data) # Gets rid of some duplicate observations

# Fix typos:
ast_data[which(ast_data[,"CHG Susceptibility designation"] == "Suceptible"),"CHG Susceptibility designation"] <- "Susceptible"
ast_data[which(ast_data[,"CHG Susceptibility designation"] == "Reduced-Suceptibility"), "CHG Susceptibility designation"] <- "Reduced-Susceptibility"

# "NA" to NA:
ast_data[which(ast_data[,"MUP Susceptibility designation"] == "NA"),"MUP Susceptibility designation"] <- NA
ast_data[which(ast_data[,"CHG Susceptibility designation"] == "NA"),"CHG Susceptibility designation"] <- NA

# Add an identifier to AST data which can be used to merge this data with Mykrobe data.
# Format: visit_site_host <- these should be unique, given that only one sample was taken per time point.

ast_data$merge_id <- paste(paste(ast_data$Visit.., ast_data$Cx.Site, sep = "_"), ast_data$ID, sep = "_")
#head(ast_data[,c("ID", "Visit..", "Cx.Site", "merge_id")])

tab <- table(ast_data$merge_id)
print(paste("Number of duplicates:", length(tab[tab > 1]))) 

# Remove duplicates:
remove_duplicates_from_ast_data <- function(ast_data){
    #' Remove duplicate rows from ast_data (multiple samples from the same visit, site and host)
    tab <- table(ast_data$merge_id)
    sus_des <- "MUP Susceptibility designation"
    for (merge_id in names(tab[tab > 1])){

        #print(ast_data[which(ast_data$merge_id == merge_id & !is.na(ast_data[, sus_des])),])

        if(dim(ast_data[which(ast_data$merge_id == merge_id & !is.na(ast_data[, sus_des])),])[1] == 1){ # One designation available, the other is NA
            ast_data[which(ast_data$merge_id == merge_id & is.na(ast_data[, sus_des])),] <- NA # Delete the NA row, keep the one with susceptibility designation available.
        }else if(dim(ast_data[which(ast_data$merge_id == merge_id & !is.na(ast_data[, sus_des])),])[1] == 2){
            ast_data[which(ast_data$merge_id == merge_id &  grepl("-1", ast_data[,"Accession.."])),] <- NA # Delete the first accession, assume accession_nr-2 is the newest sample 
        }

    }
    
    return(ast_data)
    
}
ast_data <- remove_duplicates_from_ast_data(ast_data)

 # For the remaining 15 duplicates, drop the first case. This is quite arbitrary, but necessary for later analysis and affects only a few samples.
for (vs in names(tab[tab > 1])){
    ast_data[which(ast_data$merge_id == vs),][1,] <- NA

}
tab <- table(ast_data$merge_id)
print(paste("Number of remaining duplicates:", length(tab[tab > 1]))) 

# Remove NA merge_ids:
ast_data <- ast_data[which(!is.na(ast_data$merge_id)),]
dim(ast_data)

# These should match:
length(ast_data$merge_id)
length(unique(ast_data$merge_id))


head(ast_data[,c("ID", "Visit..", "Cx.Site", "merge_id")])

## Adding missing phenotypic resistances ##

#Phenotypic susceptibility testing was not done on all isolates. If the first and last isolate of the host in question had the same phenotypic resistance designation, the isolates in between these were assumed to be of the same susceptibility designation. Here, we add resistances based on this assumption while also checking that the isolates are of the same strain. If some isolates are not the same strain as the first and last tested isolates, their susceptibility designation is assumed unknown.

# Number of missing susceptibility designations

if (debug){
    print("Number of observations with missing MUP susceptibility desgination:")
    length(ast_data[which(is.na(ast_data[,"MUP Susceptibility designation"]) | ast_data[,"MUP Susceptibility designation"] == "NA"),"ID"])
    print("Number of hosts with missing MUP susceptibility desgination:")
    length(unique(ast_data[which(is.na(ast_data[,"MUP Susceptibility designation"]) | ast_data[,"MUP Susceptibility designation"] == "NA"),"ID"]))
}
res <- "Mupirocin"
sus_designation <- "MUP Susceptibility designation"


add_genetic_res_to_ast_data <- function(ast_data, mdata, res = "Mupirocin", sus_designation = "MUP Susceptibility designation"){
    #' Add genetic resistance to AST data where phenotypic resistance was not tested. See above for a detailed explanation.
    
    tab <- table(mdata[,c("strain", res)])
    heterores_strains <- names(which(tab[,"R"] > 0 & tab[,"S"] > 0)) # Strains that have both resistant and non-resistant isolates in the genetic data.
    # Ignore these.
    new_ast_data <- ast_data

    new_ast_data$visitsite <- paste(new_ast_data$Visit.., new_ast_data$Cx.Site, sep = "_")

    strains <- unique(mdata$strain)
    n_added_resistances <- 0

    for (strain in strains){

        if (!strain %in% heterores_strains){ # Do not count strains that have both S and R isolates

            straindata <- mdata[which(mdata$strain == strain),]
            straindata$visitsite <- paste(straindata$Visit.., fupper(straindata$Cx.Site), sep = "_") 
            genetic_res <- unique(straindata[,res])

            v0 <- min(straindata$Visit..)
            v1 <- max(straindata$Visit..)

            host <- unique(straindata$ID)
            
            host_ast_data <- ast_data[which(ast_data$ID == host),]
            host_ast_data$visitsite <- paste(host_ast_data$Visit.., host_ast_data$Cx.Site, sep = "_")
            host_ast_data <- host_ast_data[which(host_ast_data$visitsite %in% straindata$visitsite),]
            
            if (dim(host_ast_data)[1] > 0){ # If the strain has phenotypic resistance data available
                  
                if (length(genetic_res) == 1){
                        to_add <- NA
                        if(is.na(genetic_res)){
                            to_add <- NA
                        }
                        else if (genetic_res == "S"){
                            to_add <- "Susceptible"
                        }
                        else if (genetic_res == "R"){
                            to_add <- "High-level R"
                        }


                        # Find the missing visits from AST data:

                        missing_res_data <- host_ast_data[which(is.na(host_ast_data[,sus_designation]) | host_ast_data[,sus_designation] == "NA"),]

                        if(dim(missing_res_data)[1] > 0){ # If there are missing sus designations in ast_data for this strain
                            for (i in 1:nrow(missing_res_data)) {
                                row <- missing_res_data[i, ]

                                v_row <- row$Visit..
                                if(v_row > v0 & v_row < v1){
                                   #print(paste("Can insert resistance:", to_add))

                                   new_ast_data[which(new_ast_data$ID == host & 
                                                      new_ast_data$visitsite == row$visitsite), sus_designation] <- to_add

                                   n_added_resistances <- n_added_resistances + 1

                                }
                            }
                        }

                }#else{print(genetic_res)}
                
            }

        }

    }
    
    print(paste("Number of susceptibility designations added:", n_added_resistances))
    
    return(new_ast_data)

}
# How this works:
# - If there are missing genetic resistances in the same strain, no resistances are added to the phenotypic resistances.
# - If there are differing sus designations associated with the same strain, no resistances added to the phenotypic data.

## Finally, merge into a data frame with phenotypic and genetic resistances ##

# First, add some missing resistances to AST data:

ast_data <- add_genetic_res_to_ast_data(ast_data, mdata, res = "Mupirocin", sus_designation = "MUP Susceptibility designation")
ast_data <- add_genetic_res_to_ast_data(ast_data, mdata, res = "Chlorhexidine", sus_designation = "CHG Susceptibility designation")

# Then merge mdata and ast_data using inner join on merge_id:
mdata$merge_id <- paste(paste(mdata$Visit.., fupper(mdata$Cx.Site), sep = "_"), mdata$ID, sep = "_")
merged_ast_data <- merge(mdata, ast_data, by = "merge_id")

# Fix some column names:
colnames(merged_ast_data)[colnames(merged_ast_data) == "ID.x"] <- "ID"
colnames(merged_ast_data)[colnames(merged_ast_data) == "Visit...x"] <- "Visit.."
colnames(merged_ast_data)[colnames(merged_ast_data) == "Cx.Site.x"] <- "Cx.Site"
colnames(merged_ast_data)[colnames(merged_ast_data) == "Accession...x"] <- "Accession.."

# Dimensions of merged and original data sets:
dim(merged_ast_data)
dim(ast_data)
dim(mdata)

# Percent of mdata that also has phenotypic mupirocin resistance available:
# - This is directly the dimensions after the merge.
# - Mykrobe data = genetic resistances

if (debug){
    print("Proportion of Mykrobe data that also has phenotypic resistance data (NA's included):")
    dim(merged_ast_data)[1]/dim(mdata)[1] 
    print("Proportion of Mykrobe data that has non-NA phenotypic resistance available for mupirocin:")
    dim(merged_ast_data[which(!is.na(merged_ast_data[,"MUP Susceptibility designation"])),])[1]/dim(merged_ast_data)[1] # 81% 
    dim(merged_ast_data[which(!is.na(merged_ast_data[,"MUP Susceptibility designation"])),])[1]/dim(mdata)[1] # 80%
}




## Add HLR and LLR columns to mdata:
# High-level resistance:
merged_ast_data$mupirocin_hlr <- merged_ast_data[,"MUP Susceptibility designation"]
merged_ast_data[which(merged_ast_data$mupirocin_hlr == "High-level R"),"mupirocin_hlr"] <- 1
merged_ast_data[which(merged_ast_data$mupirocin_hlr == "Low-level R" | merged_ast_data$mupirocin_hlr == "Susceptible"), "mupirocin_hlr"] <- 0

# Low-level resistance:
merged_ast_data$mupirocin_llr <- merged_ast_data[,"MUP Susceptibility designation"]
merged_ast_data[which(merged_ast_data$mupirocin_llr == "Low-level R"),"mupirocin_llr"] <- 1
merged_ast_data[which(merged_ast_data$mupirocin_llr == "High-level R" | merged_ast_data$mupirocin_llr == "Susceptible"), "mupirocin_llr"] <- 0

mdata_ast <- merged_ast_data

# Only include unnecessary columns
preprocessed_ast_data <- merged_ast_data[,c("ID", "Visit..", "Cx.Site", "MUP Susceptibility designation",
                                      "CHG Susceptibility designation", "strain",
                                      "mupirocin_llr", "mupirocin_hlr", "Mupirocin", "Chlorhexidine")]

# Save the data
saveRDS(preprocessed_ast_data, file = "data/ast_data.RData")


### Survival data from phenotypic resistances ###

source(paste(root, "scripts/surv_preprocessing_strains.R", sep = ""))

ts <- c(treatments, "mupirocin_hlr", "mupirocin_llr")

mdata_ast[which(mdata_ast$mupirocin_hlr == "NA"), "mupirocin_hlr"] <- NA
mdata_ast[which(mdata_ast$mupirocin_llr == "NA"), "mupirocin_llr"] <- NA

mdata_ast <- add_missing_resistance(mdata_ast)  # add missing resistances if applicable. Only for genetic resistance!

visit_data <- load_visit_data(visit_data_file)
mykrobe_data <- preprocess_mykrobe_data()

C <- get_clearance_data(visit_data)
C <- modify_clearmat_wound_6_to_1(C)

dfst5 <- readRDS(paste(data_root, "dfst5.RData", sep = ""))
dfst8 <- readRDS(paste(data_root, "dfst8.RData", sep = ""))

dfst <- add_threshold_strains(rbind(dfst5, dfst8), threshold = THRESHOLD)

data_ast <- modify_data(get_surv_data_strain(mdata_ast, treatments = ts)) 

data_ast <- remove_non_adh_intervals(data_ast, data_root)

# Save the survival data

saveRDS(data_ast, file = "data/surv_data/data_ast.RData")

