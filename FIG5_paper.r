# ---------------------------------------------
# Script to produce Fig. 5 (left and central panel) of BayVel paper with comparison
# scVelo results on real data, projected with PCA
# ----------------------------------------------
#
# OUTPUTS
# - 2 PDF files, for left and central panel in Figure 5
#
# INPUTS (set at the beginning):
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
# ---------------------------------------------


rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set directory with BayVel results and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
pathToResults <- paste0(pathToYourDirectory, "/real data")
setwd(paste0(pathToYourDirectory))
source(paste0(pathToYourDirectory, "/functions.R"))
# Set the output path where you will save the images
pathOutput <- paste0(pathToYourDirectory, "/figuresPaper/")

# -----------------------------
#  PACKAGES 
# -----------------------------
library(ggplot2)

# -----------------------------
# TYPE OF SIMULATION WE ARE CONSIDERING
# -----------------------------
n_genes <- 2000                   # number of genes
typeSIM <- "Pancreas"             # specify that we are considering simulated data
addScVelo <- TRUE                 # specify if you want to plot the results of scVelo

# ------------------------------
# LOAD SCVELO DATA
# ------------------------------
nameSim <- "moments"
typeSW <- "SW1"
scv <- new.env(parent = baseenv()) 
source(paste0(pathToYourDirectory, "/loadScVelo_results.r"))


# ------------------------------
# PLOTS
# ------------------------------

pc = 2 # PCA component

# create an unique df for the plots
dfMs <- data.frame(
    pcI = rep(NA, n_cells * (1 + addScVelo)), 
    pcII = NA, 
    pcI_prime_norm = NA, 
    pcII_prime_norm = NA, 
    color = c(as.factor(scv$typeCell), rep("blue", n_cells*addScVelo)), 
    fill = c(as.factor(scv$typeCell), rep(scv$predPCA$typeCell, addScVelo)), 
    type = c(rep("Ms", n_cells), rep("scVelo", n_cells*addScVelo)), 
    typeCell = rep(as.factor(scv$typeCell), 1 + addScVelo), 
    subtypeCell = rep(as.factor(scv$subtypeCell), 1 + addScVelo)
)

# compute PCA on the pre-processed Ms 
scv$pcaMs <- prcomp(scv$Ms, scale = TRUE)
scv$pcaMs_df <- data.frame(scv$pcaMs$x)
scv$pcaMs_df$typeCell <- as.factor(scv$typeCell)
scv$pcaMs_df$subtypeCell <- as.factor(scv$subtypeCell)

# compute S' = S_estim + v*dt (using derivative definition)
scv$v <- sweep(scv$pos_u, 2, scv$beta, "*") - sweep(scv$pos_s, 2, scv$fit_gamma, "*")
dt <- 0.001
scv$s_prime <- scv$pos_s + scv$v*dt

# project s and s' into the PCA projection
colnames(scv$pos_s) <- colnames(scv$Ms)
scv$predPCA <- predict(scv$pcaMs, newdata = scv$pos_s)
scv$predPCA <- data.frame(scv$predPCA)
scv$predPCA$subtypeCell <- as.factor(scv$subtypeCell)
scv$predPCA$typeCell <- as.factor(scv$typeCell)

colnames(scv$s_prime) <- colnames(scv$Ms)
scv$predPCA_prime <- predict(scv$pcaMs, newdata = scv$s_prime)
scv$predPCA_prime <- data.frame(scv$predPCA_prime)
scv$predPCA_prime$subtypeCell <- as.factor(scv$subtypeCell)
scv$predPCA_prime$typeCell <- factor(scv$typeCell)

# normalize the length of the PCA arrows
norm <- 9
scv$predPCA_prime$norm <- sqrt((scv$predPCA$PC1 - scv$predPCA_prime$PC1)^2 + (scv$predPCA$PC2 - scv$predPCA_prime$PC2)^2)
scv$predPCA_prime$PC1norm <- scv$predPCA$PC1 + (scv$predPCA_prime$PC1 - scv$predPCA$PC1)*norm/scv$predPCA_prime$norm
scv$predPCA_prime$PC2norm <- scv$predPCA$PC2 + (scv$predPCA_prime$PC2 - scv$predPCA$PC2)*norm/scv$predPCA_prime$norm

# prepare the dataset
dfMs[which(dfMs$type == "Ms"), c("pcI", "pcII")] <- scv$pcaMs_df[, c("PC1", paste0("PC", pc))]      # M_s
dfMs[which(dfMs$type == "scVelo"), c("pcI", "pcII")] <- scv$predPCA[, c("PC1", paste0("PC", pc))]   # s 
dfMs[which(dfMs$type == "scVelo"), c("pcI_prime_norm", "pcII_prime_norm")] <- scv$predPCA_prime[, c("PC1norm", paste0("PC", pc, "norm"))] # s'


# graphics parameter
dfMs$alphaArrow <- ifelse(dfMs$type == "scVelo", 0.8, 1) # arrows transparency
gg <- ggplot(data = dfMs) +
geom_point(aes(x = pcI, y = pcII, fill = fill, color = color))
plot_build <- ggplot_build(gg)
col <- unique(plot_build$data[[1]]$fill)
dfMs$fill <- rep(col[scv$typeCell], sum(1 + addScVelo))
dfMs$color <- c(col[scv$typeCell], rep(alpha("blue", .8), n_cells*addScVelo))

fillLegendMs <- c("white", rep("blue", 1*addScVelo))

# -------------------------------------------------
# FIG. 5a
# -------------------------------------------------
# create the plot for the positions
ggMs <- ggplot(data = dfMs) + 
        geom_point(data = subset(dfMs, type == "Ms"), aes(x = pcI, y = pcII, shape = type, fill = fill, color = fill), alpha = 0.5, size = 1) + 
        geom_point(data = subset(dfMs, type != "Ms"), aes(x = pcI, y = pcII, shape = type, fill = fill, color = color), alpha = 0.7, size = 3) +
        # internal colors (inner shapes)
        scale_fill_manual(values = sort(unique(dfMs$fill)), labels = unique(scv$typeCellReal)[order(unique(dfMs$fill))]) +
        # external colors (margin)
        scale_color_manual(values = sort(unique(dfMs$color))) +
        # shape
        scale_shape_manual(values = c( "scVelo" = 24, "Real" = 22, "Ms" = 21)) + 
        # legends
        guides(
            shape = guide_legend(override.aes = list(fill = fillLegendMs, size = 7)),
            fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 7)),
            color = "none",  
            alpha = "none",
            size = "none"
        ) + 
        labs(x = "PCA 1", y = paste0("PCA ", pc), fill = "Groups", shape = "Methods", color = "Color Legend") 

ggMs <- ggMs +
        # x-axis
        scale_x_continuous(breaks = c(-40, -20, 0, 20, 40), labels = c("-40", "-20", "0", "20", "40")) + 
        xlim(c(-40, 40)) +
        # graphics parameters        
        theme(
            axis.text = element_text(family = "serif", size = 35), 
            axis.title = element_text(family = "serif", size = 45), 
            legend.text = element_text(family = "serif", size = 35), 
            legend.title = element_text(family = "serif", size = 37), 
            legend.position = c(0.27,0.35), 
            legend.key = element_rect(colour = "transparent", fill = NA), 
            legend.background=element_blank(),
            plot.margin = unit(c(0.2,2.5,0,0.2), "lines"),
            plot.title = element_text(family = "serif", size=50,  hjust = 0.5), 
        ) 

pdf(paste0(pathOutput, "/Fig5_a.pdf"), width = 7*1.3, height = 7*1.5)
    print(ggMs)
dev.off()
        



# -------------------------------------------------
# FIG. 5b
# -------------------------------------------------
# add the velocity to the plot     

dfMs_scvelo <- dfMs[which(dfMs$type != "Ms"), ]
# add the arrow
for(c in 1:n_cells){
    l <- 0.01
    if(dfMs_scvelo$pcI[c] < dfMs_scvelo$pcI_prime_norm[c]){
        while(dfMs_scvelo$pcI[c] + l  < dfMs_scvelo$pcI_prime_norm[c]){
            dfMs_scvelo$pcI_terminalArrow[c] <- dfMs_scvelo$pcI[c] + l
            l <- l + 0.01
        }
    }else{
        while(dfMs_scvelo$pcI[c] - l  > dfMs_scvelo$pcI_prime_norm[c]){
            dfMs_scvelo$pcI_terminalArrow[c] <- dfMs_scvelo$pcI[c] - l
            l <- l + 0.01
        }
    }
}
# sample just 15 cells for each type in order to allow the plot to be understandable
set.seed(123)
sampledCell <- c()
for(ty in 1:max((scv$typeCell))){
    type <- which(scv$typeCell == ty)
    sampledCell <- c(sampledCell, sample(type, min(length(type), 15)))
}
dfMs_scvelo <- dfMs_scvelo[sort(sampledCell), ]

# add the arrows to the plot
k <- 0.7
ggMs <- ggMs + 
        geom_segment(data = dfMs_scvelo, aes(x = pcI, y = pcII, xend = pcI_prime_norm, yend = pcII_prime_norm, color = fill), color = alpha("black", 0.5), arrow = arrow(length = unit(0.7, "cm")), linewidth = 1.3*k) + 
        geom_segment(data = subset(dfMs_scvelo, type != "Ms"), aes(x = pcI, y = pcII, xend = pcI_prime_norm, yend = pcII_prime_norm, color = fill), linewidth = 0.75*k, arrow = arrow(length = unit(0.7, "cm")), alpha = 0.5) 

pdf(paste0(pathOutput, "/Fig5_b.pdf"), width = 7*1.3, height = 7*1.5)
    print(ggMs)
dev.off()



  



    
    