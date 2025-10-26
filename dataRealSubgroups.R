# -----------------------------------------------------------
# Script to obtain subgroup labels for the real data, based on hierarchical clustering. 
# 
# OUTPUTS: .RData file, with subgroup labels coinciding with groups labels
#          .RData file, with subgroup labels obtained from hierarchical clustering
#          .pdf file with the dendograms obtained from hierarchical clustering (optional, see plotDendogram option)
#          .pdf file with the graphic representation of subgroup division of cells.
#
# INPUT: 
# - typeSIM:            Type of real data (Pancreas or DentateGyrus)
# - typeProcessing:     Pre-processing of scVelo functions that are still performed (: filter, filter_and_normalize, filter_and_normalize_noLog, moments)
# - plotDendogram:      Specify if we want to plot the dendograms obtained from the hierarchical clustering or not
# 
# DEPENDENCIES:
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: cluster
#             ggplot2
# -----------------------------------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
# PACKAGES
# -----------------------------
library(cluster)
library(ggplot2)

# -----------------------------
# INPUT
# -----------------------------
typeSIM <- "Pancreas"                  # type of data (Pancreas, DentateGyrus)
typeProcessing <- "filter"             # Specify which scVelo pre-processing steps are still performed (options: filter, filter_and_normalize, filter_and_normalize_noLog, moments)
plotDendogram <- TRUE                  # Boolean, deciding if we want or not to plot the dendogram from the hierarchical clustering
# -----------------------------
#  PATH 
# -----------------------------
# Set working directory ad load auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
chrTypeSIM <- ifelse(typeSIM == "sim", "simulations", "real data")
source(paste0(pathToYourDirectory, "/functions.R"))

# -----------------------------
# LOAD THE DATA 
# -----------------------------
res <- loadRealData(typeProcessing = typeProcessing, pathData = paste0(pathToYourDirectory, "/", chrTypeSIM, "/", typeSIM)) # import real data, where we have just filtered out data that are not enough expressed
Y_u <- res$unspliced
Y_s <- res$spliced
typeCell <- res$typeCell

n_genes <- dim(Y_u)[2]

# -----------------------------
# PCA PROJECTION
# -----------------------------
set.seed(seed)

data.pca <- prcomp(cbind(Y_u, Y_s), scale = TRUE) # Apply PCA considering together unspliced and spliced data 
pcaYus_df <- data.frame(data.pca$x)
pcaYus_df$typeCell <- as.factor(typeCell)
pcaYus_df$typeCellInt <- as.integer(pcaYus_df$typeCell)


pcaYus_dfFINAL <- pcaYus_df
pcaYus_dfFINAL$subtypeCell <- pcaYus_dfFINAL$typeCellInt

typeCell <- pcaYus_dfFINAL$typeCellInt
subtypeCell <- pcaYus_dfFINAL$typeCellInt
typeCellT0_off <- pcaYus_dfFINAL$typeCellInt

# groups coincides with subgroups
save.image(paste0(pathToYourDirectory, "/", chrTypeSIM, "/", typeSIM, "/", typeProcessing, "/", typeSIM, n_genes, "genes_1subgr.RData"))

# -----------------------------
# HIERARCHICAL CLUSTERING
# -----------------------------
set.seed(seed)

# considering different clustering methods to be tested
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute the agglomerative coefficient for a given distance matrix and method
ac <- function(dist_mat, x) {
  return(agnes(dist_mat, diss = TRUE, method = x)$ac)
}
ac_vec <- rep(NA, length(m))


if(plotDendogram){
    pdf(paste0(pathToYourDirectory, "/", chrTypeSIM, "/", typeSIM, "/", typeProcessing, "/hierClust_dendogram.pdf"), height = 7*1.5, width = 7*1.2)
    par(mfrow = c(1,1))
}

last <- 0  
for(ty in 1:length(unique(pcaYus_dfFINAL$typeCellInt))){     # loop over each group
    cellTy <- which(pcaYus_dfFINAL$typeCellInt == ty)

    # compute Euclidean distance matrix using the first two PCA components
    dist_mat <- dist(pcaYus_dfFINAL[cellTy,1:2], method = 'euclidean')
    
    # Evaluate each grouping method by its agglomerative coefficient
    for(i in 1:length(m)){
        ac_vec[i] <- ac(dist_mat, m[i])
    }
    # Select the method with the maximum agglomerative coefficient
    meth <- m[which.max(ac_vec)]
    
    # Perform agglomerative clustering with the selected method
    aggClust <- agnes(dist_mat, diss = TRUE, method = meth)

    # Plot the dendogram
    if(plotDendogram){
        pltree(aggClust, main = pcaYus_dfFINAL$typeCell[pcaYus_dfFINAL$typeCellInt == ty][1]) 
    }

    # Set the initial number of subgroups depending on sample size (we ensure there are enough cells in each subgroup)
    if(length(cellTy) > 30){
        n_clusters <- 10
        
    }else{
        n_clusters <- 1
    }
    
    # Iteratively adjust the number of subgroups 
    # Rule: all subgroups must have at least 30 cells
    stop <- FALSE   
    while(!stop){
        cut_min <- cutree(aggClust, k = n_clusters)
        
        if(n_clusters == 1){
            stop <- TRUE
        }else if(any(table(cut_min) < 30)){ # if any subgroup has fewer than 30 cells reduce the number of clusters
            n_clusters <- n_clusters - 1
        }else{
            stop <- TRUE
        }
    }
        
    # Assign subgroups labels 
    pcaYus_dfFINAL[cellTy,]$subtypeCell <- (cut_min) + last
    last <- last + n_clusters
}         
dev.off()


pdf(paste0(pathToYourDirectory, "/", chrTypeSIM, "/", typeSIM, "/", typeProcessing, "/hierClust.pdf"), height = 7*1.5, width = 7*1.2)
par(mfrow = c(1, 1))
for(ty in unique(pcaYus_dfFINAL$typeCellInt)){
    cellTy <- which(pcaYus_dfFINAL$typeCellInt == ty)
    gg <- ggplot(pcaYus_dfFINAL[cellTy, ], aes(x = PC1, y = PC2, col = as.factor(subtypeCell))) + geom_point(size = 3) + labs(title = pcaYus_dfFINAL$typeCell[pcaYus_dfFINAL$typeCellInt == ty][1], col = "Subgroups") + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), plot.title = element_text(hjust = 0.5, size = 30), legend.title = element_text(size = 25), legend.text = element_text(size = 25), legend.position = "bottom")
    plot(gg)
}
dev.off()

save.image(paste0(pathToYourDirectory, "/", chrTypeSIM, "/", typeSIM, "/", typeProcessing, "/", typeSIM, n_genes, "genes_multsubgr.RData"))
