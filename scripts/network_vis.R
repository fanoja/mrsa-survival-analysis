# Visualize the colonization of a single host with a network figure
source("scripts/utilities.R")

# Network figure of within-host colonization

n_rand_hex_cols <- function(n){
    #' Get n random hex colors.
    #'
    #' @param n Number of colors to generate.
    #' @return hex_cols A vector of n hex values.

    hex_cols <- c()
    for (i in 1:n){
        char <- paste("#", paste(sample(c(0:9,letters[1:6]), 6), collapse = ""), sep="")
        while(T){ # No duplicate colors allowed
            if (char %in% hex_cols){
            char <- paste("#", paste(sample(c(0:9,letters[1:6]), 6), collapse = ""), sep="")
            }
            else {
            hex_cols <- c(hex_cols, char)
            break
            }
        }
    }
    
    return(hex_cols)
}

vec_to_col <- function(vec){
    #' Creates random color for each unique value of the vector.
    #'
    #' Converts each value in vec to a color (same values = same color) and returns this color vector.
    #'
    #' @param vec A vector.
    #' @return cols A vector containing colors associated with each unique elemnt in vec.
    
    vec <- as.character(vec)
    colmap <- n_rand_hex_cols(length(unique(vec)))
    names(colmap) <- as.character(unique(vec))
    
    cols <- c()
    for (v in vec){
        if (v %in% names(colmap)){
          cols <- c(cols, colmap[[v]])
        }
    }

    return(cols)
}

charvec_to_letter <- function(charvec){
    #' Assign a letter from a to z to elements of charvec. One letter per unique element.
    #'
    #' Returns a vector with letters. Length(vector) == length(charvec)
    #'
    #' @param charvec A character vector.
    
    dict <- list()

    i <- 1
    for (unique_char in unique(charvec)){
        dict[unique_char] <- letters[as.integer(i)]
        i <- i + 1
    }

    for (char in charvec){
        charvec[charvec == char] <- dict[char]
    }

    return(as.character(charvec))

}

get_layout <- function(visits, sites){
    #' Get new coordinates for the layout of a host colonization network figure.
    #'
    #' @param visits Visits of interest.
    #' @param sites Sites of interest.
    #' @return new_layout Coordinates for nodes in host colonization network figure.
    
    x <- c()
    for (v in visits){
        if (v == "R0"){
          x <- c(x,0)
        }
        if (v == "V1"){
          x <- c(x,1)
        }
        if (v == "V2"){
          x <- c(x,2)
        }
        if (v == "V3"){
          x <- c(x,3)
        }
        if (v == "V4"){
          x <- c(x,4)
        }
    }

    y <- c()

    for (s in sites){
        if (s == "wound"){
          y <- c(y,0)
        }
        if (s == "skin"){
          y <- c(y,1)
        }
        if (s == "throat"){
          y <- c(y,2)
        }
        if (s == "nares"){
          y <- c(y,3)
        }

    }

    new_layout <- matrix(data = c(x,y), nrow = length(x), ncol = 2)

    return(new_layout)
}


vis_mrc_df_6 <- function(df, rmv_edges = T, title_addon = "default", vs = 20){
    #' Create a network visualization of the colonization of a host over the study period. 
    #'
    #' A node represents an MRSA positive sample (isolate) at a visit. Node color and letter denote the strain of the MRSA isolate. y-axis is the body site the MRSA isolate was sampled from. Pass "default" to title_addon if you wish to have the id number as title_addon.
    #'
    #' @param df Consecutive distance data frame
    #' @param rmv_edges If true, removes edges that have an SNP distance below a certain threshold.
    #' @param title_addon Add a string to the title of the figure, such as host id.
    #' @param vs Vertex size.
     
    load_library("igraph")
    
    G <- graph_from_data_frame(df, directed = T)

    if (rmv_edges){
     G <- G - E(G)[!E(G)$colgroup] # Remove nonexistent edges (edges with a distance below the SNP threshold. By default 20)
    }

    clstrs <- clusters(G) # Get the clusters)
    nodestrains <- unique(rbind(cbind(df$Var1, df$substrain1), cbind(df$Var2, df$substrain2)))[,2]
    nodecols <- vec_to_col(nodestrains)

    visits <- unique(rbind(cbind(df$site0, df$v0), cbind(df$site1, df$v1)))[,2]#unique(c(df$visit1,df$visit2))#c()
    sites <- unique(rbind(cbind(df$site0, df$v0), cbind(df$site1, df$v1)))[,1] # Not sure if this works :) 
    
    L <- get_layout(visits, sites) # Get a new layout for the network graph
    
    vertex_labels <- charvec_to_letter(nodestrains)

    if(title_addon == "default"){
        title_addon <- paste(", Host", unique(df$host), sep = " ")
     }

    # Color edge and vertex frame by strain for clarity:
    unique_nodes <- unique(nodestrains)
    cols <- c("black", "cadetblue3", "brown2", "blueviolet", "darkgoldenrod2")#brewer.pal(length(unique_nodes), "Dark2")#c("red", "blue", "black", "grey", "green", "pink", "darkyellow")
    nodecols <- cols[match(nodestrains, unique_nodes)]

    V(G)$color <- nodecols
    E(G)$color <- tail_of(G, E(G))$color
    

    #par(font = 1, family = "Times", ps = 12)
    #par(ps = 12)
    # vertex.size = 20 to adjust node size
    plot(G, asp = 0, xlim = c(0, max(L[,1])), ylim = c(0, max(L[,2])), main = paste("Within-host Colonization", title_addon, sep = ""), rescale = F, vertex.label.color = "black", vertex.color = "white", vertex.frame.color = V(G)$color, vertex.size = vs, edge.color = E(G)$color, vertex.label = toupper(vertex_labels), layout = get_layout(visits, sites), edge.arrow.mode = 0)


    # Match visits with labels (0-month, 1-month etc):
    visits <- unique(visits)
    labels <- c("0-month", "1-month", "3-month", "6-month", "9-month")
    x_axis_labels <- labels[match(visits,c("R0", "V1", "V2", "V3", "V4"))]

    #axis(1, at = unique(L[,1]), labels = x_axis_labels) # Use this if only relevant visits should be visible. Can be messy!
    axis(1, at = c(0,1,2,3,4), labels = labels) # Fixed x-axis for all hosts.

    axis(2, at = c(3,2,1,0), labels = c("nares", "throat", "skin", "wound"))
  
}
