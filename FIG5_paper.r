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
#             latex2exp
#             RColorBrewer
#             grid
# ---------------------------------------------


rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set directory with BayVel results and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
pathToResults <- "pathToYourDirectory/real data"
setwd(paste0(pathToYourDirectory))
source(paste0(pathToYourDirectory, "/functions.R"))
# Set the output path where you will save the images
pathOutput <- paste0(pathToYourDirectory, "/figuresPaper/")

# -----------------------------
#  PACKAGES 
# -----------------------------
library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(grid)


# -----------------------------
# TYPE OF SIMULATION WE ARE CONSIDERING
# -----------------------------
n_genes <- 2000                   # number of genes
typeSIM <- "Pancreas"             # specify that we are considering simulated data
addScVelo <- TRUE                 # specify if you want to plot the results of scVelo

# PCA
pc = 2

pdf(paste0(pathOutput, "/Fig5_a.pdf"))

# create an unique df

    dfMs <- data.frame(
        pcI = rep(NA, env$n_cells * (1 + addScVelo + addReal)), 
        pcII = NA, 
        pcI_prime_norm = NA, 
        pcII_prime_norm = NA, 
        color = c(as.factor(env$typeCell), rep("blue", env$n_cells*addScVelo), rep("red", env$n_cells*addReal)), 
        fill = c(as.factor(env$typeCell), rep(scVelo$predPCA$typeCell, addReal), rep(real$predPCA_Ms$typeCell, addReal)), 
        type = c(rep("Ms", env$n_cells), rep("scVelo", env$n_cells*addScVelo), rep("Real", env$n_cells*addReal)), 
        typeCell = rep(as.factor(env$typeCell), 1 + addScVelo + addReal), 
        subtypeCell = rep(as.factor(env$subtypeCell), 1 + addScVelo + addReal)
    )

   
    dfMs[which(dfMs$type == "Ms"), c("pcI", "pcII")] <- env$pcaMs_df[, c("PC1", paste0("PC", pc))]


    if(addScVelo){
        dfMs[which(dfMs$type == "scVelo"), c("pcI", "pcII")] <- scVelo$predPCA[, c("PC1", paste0("PC", pc))]
        dfMs[which(dfMs$type == "scVelo"), c("pcI_prime_norm", "pcII_prime_norm")] <- scVelo$predPCA_prime[, c("PC1norm", paste0("PC", pc, "norm"))]
    }

    
    dfMs$alphaArrow <- ifelse(dfMs$type == "scVelo", 0.8, 1)
    

    # select better colors
    p <- ggplot(data = dfYs) +
    geom_point(aes(x = pcI, y = pcII, fill = fill, color = color))
    plot_build <- ggplot_build(p)
    colorazione <- unique(plot_build$data[[1]]$fill)


    dfMs$fill <- rep(colorazione[env$typeCell], sum(1 + addScVelo + addReal))
    dfMs$color <- c(colorazione[env$typeCell], rep(alpha("blue", .8), n_cells*addScVelo), rep("red", n_cells*addReal))

    fillLegendYs <- c(rep("darkgreen", 1*addBayVel), rep("red", 1*addReal), "white")
    fillLegendMs <- c("white", rep("blue", 1*addScVelo), rep("red", 1*addReal))



    # -------------------------------------------------
    # create the plot for the positions
    ggMs <- ggplot(data = dfMs) +  geom_point(data = subset(dfMs, type == "Ms"), aes(x = pcI, y = pcII, shape = type, fill = fill, color = fill), alpha = 0.6, size = 0.5) +  geom_segment(data = subset(dfMs, type != "Ms"), aes(x = pcI, y = pcII, xend = pcI_prime_norm, yend = pcII_prime_norm, color = color, alpha = alphaArrow), linewidth = 0.3, arrow = arrow(length = unit(0.1, "cm"))) +        geom_point(data = subset(dfMs, type != "Ms"), aes(x = pcI, y = pcII, shape = type, fill = fill, color = color), alpha = 0.9, size = 1) +
        # Scale per i colori di riempimento (interno)
        scale_fill_manual(values = sort(unique(dfMs$fill)), labels = unique(scVelo$typeCellReal)[order(unique(dfMs$fill))]) +
        # Scale per i colori esterni (bordo)
        scale_color_manual(values = sort(unique(dfMs$color))) +
        # Scale per la forma
        scale_shape_manual(values = c( "scVelo" = 24, "Real" = 22, "Ms" = 21)) 

    
    # Legende separate con colori personalizzati per la legenda di type
    ggMs <- ggMs + guides(
                shape = guide_legend(override.aes = list(fill = fillLegendMs, size = 3)),
                fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 3)),
                color = "none",  # Nascondi la leggenda del colore esterno se non necessaria
                alpha = "none",
                size = "none"
    )

    ggMs <- ggMs +
        # Titoli per le legende
        labs(x = "PC1", y = paste0("PCA ", pc), fill = "Types", shape = "Methods", color = "Color Legend", title = TeX(paste0("PCA projection of scVelo results"))) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        themeGGPLOT +
        themeLEGEND 

    # Mostra il plot
    print(ggMs)


    # ---------------------------------------
    # add the velocity to the plot     
    ggMs <- ggMs + geom_segment(data = subset(dfMs, type!= "Ms"), aes(x = pcI, y = pcII, xend = pcI_prime_norm, yend = pcII_prime_norm, color = color, alpha = alphaArrow), arrow = arrow(length = unit(0.2, "cm")))

    pdf(paste0(pathSavePlot, "/",  nome, "/pca1", pc, "-CFRseparated.pdf"))
    
    print(ggYs)
    print(ggMs)

dev.off()
