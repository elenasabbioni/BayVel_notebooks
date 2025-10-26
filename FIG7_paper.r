# ---------------------------------------------
# Script to produce Fig. 7 of BayVel paper with BayVel results on real data
# ----------------------------------------------
#
# OUTPUTS
# - 4 .pdf files with plots in Fig. 7
#
# INPUTS (set at the beginning):
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
#             LaplacesDemon
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
library(LaplacesDemon)
library(ggplot2)

# -----------------------------
# TYPE OF SIMULATION WE ARE CONSIDERING
# -----------------------------
n_genes <- 2000                   # number of genes
typeSIM <- "Pancreas"             # specify that we are considering simulated data
addBayVel <- TRUE                 # specify if you want to plot the results of BayVel

# ------------------------------
# LOAD SCVELO DATA
# ------------------------------
nameSim_vec <- c("SW1-T1-D4", "SW2-T3-D4")
mcmc <- 250000

for(nameSim in nameSim_vec){
    bv <- new.env(parent = baseenv()) 
    maxLog <- TRUE
    source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))   

    source(paste0(pathToYourDirectory, "/functions.R"))

    # ------------------------------
    # PCA projection of s and s' = s + dt*v
    # ------------------------------
    norm <- 5
    set.seed(seed)
    # compute PCA on the Ys/Catt
    bv$pcaYs <- prcomp(bv$Y_s, scale = TRUE)
    bv$pcaYs_df <- data.frame(bv$pcaYs$x)
    bv$pcaYs_df$typeCell <- as.factor(bv$typeCell)
    bv$pcaYs_df$subtypeCell <- as.factor(bv$subtypeCell)

    # compute S' = S_estim + v*dt
    dt <- 0.001
    bv$s_prime <- bv$s_chain + bv$v*dt

    # project s into the PCA projection
    bv$u_chainNEW <- bv$u_chain[bv$subtypeCell, ,]
    bv$s_chainNEW <- bv$s_chain[bv$subtypeCell, ,]
    bv$mean_sNB <- sweep(bv$s_chainNEW, 1, bv$Catt_chain, "*")

    colnames(bv$s_chain) <- colnames(bv$Y_s)
    colnames(bv$mean_sNB) <- colnames(bv$Y_s)

    bv$predPCA <- predict(bv$pcaYs, newdata = bv$mean_sNB)
    bv$predPCA <- data.frame(bv$predPCA)
    bv$predPCA$subtypeCell <- as.factor(bv$subtypeCell)
    bv$predPCA$typeCell <- as.factor(bv$typeCell)

    # project s' into the PCA projection
    bv$mean_sNB_prime <- sweep(bv$s_prime[bv$subtypeCell, ,], 1, bv$Catt_chain, "*")

    colnames(bv$s_prime) <- colnames(bv$Y_s)
    colnames(bv$mean_sNB_prime) <- colnames(bv$Y_s)

    bv$predPCA_prime <- predict(bv$pcaYs, newdata = bv$mean_sNB_prime)
    bv$predPCA_prime <- data.frame(bv$predPCA_prime)
    bv$predPCA_prime$subtypeCell <- as.factor(bv$subtypeCell)
    bv$predPCA_prime$typeCell <- as.factor(bv$typeCell)

    # normalize the length of the velocity arrows
    bv$predPCA_prime$norm <- sqrt((bv$predPCA$PC1 - bv$predPCA_prime$PC1)^2 + (bv$predPCA$PC2 - bv$predPCA_prime$PC2)^2)
    bv$predPCA_prime$PC1norm <- bv$predPCA$PC1 + (bv$predPCA_prime$PC1 - bv$predPCA$PC1)*norm/bv$predPCA_prime$norm
    bv$predPCA_prime$PC2norm <- bv$predPCA$PC2 + (bv$predPCA_prime$PC2 - bv$predPCA$PC2)*norm/bv$predPCA_prime$norm

    # find an application point for the velocity for each subgroup
    bv$applicationPoint <- apply(bv$predPCA[,1:2000], MARGIN = 2, FUN = function(x) tapply(x, bv$subtypeCell, FUN = mean))
    bv$applicationPoint <- data.frame(bv$applicationPoint)
    
    mapping <- c()
    for(sty in 1:bv$n_subtypeC){
        mapping[sty] <- which(bv$subtypeCell == sty)[1]
    }
    bv$applicationPoint$PC1norm <- bv$applicationPoint$PC1 + (bv$predPCA_prime$PC1norm[mapping] - bv$predPCA$PC1[mapping])    
    bv$applicationPoint$PC2norm <- bv$applicationPoint$PC2 + (bv$predPCA_prime$PC2norm[mapping] - bv$predPCA$PC2[mapping])    


    # ------ PLOT FOR POSITIONS    
    pc <- 2
    # create an unique df
    dfYs <- data.frame(
        pcI = NA, 
        pcII = NA, 
        pcI_prime_norm = NA, 
        pcII_prime_norm = NA, 
        color = c(as.factor(bv$typeCell), rep("darkgreen", length(bv$subtypeCell)*addBayVel)),
        fill = c(as.factor(bv$typeCell), rep(bv$predPCA$typeCell, addBayVel)), 
        type = c(rep("Ys", bv$n_cells), rep("BayVel", length(bv$subtypeCell)*addBayVel)), 
        typeCell = c(as.factor(bv$typeCell), rep(bv$predPCA$typeCell, addBayVel)), 
        subtypeCell = c(as.factor(bv$subtypeCell), rep(bv$predPCA$subtypeCell, addBayVel))
    )
    dfYs[which(dfYs$type == "Ys"), c("pcI", "pcII")] <- bv$pcaYs_df[, c("PC1", paste0("PC", pc))]
    dfYs[which(dfYs$type == "BayVel"), c("pcI", "pcII") ] <- bv$predPCA[, c("PC1", paste0("PC", pc))]
    dfYs[which(dfYs$type == "BayVel"), c("pcI_prime_norm", "pcII_prime_norm")] <- bv$predPCA_prime[, c("PC1norm", paste0("PC", pc, "norm"))]

    bv$applicationPoint[, c("pcI", "pcII", "pcI_prime_norm", "pcII_prime_norm")] <- bv$applicationPoint[, c("PC1", paste0("PC", pc), "PC1norm", paste0("PC", pc, "norm"))]
    
    # select better colors
    gg <- ggplot(data = subset(dfYs, type == "Ys")) +
    geom_point(aes(x = pcI, y = pcII, fill = typeCell, color = typeCell))
    plot_build <- ggplot_build(gg)
    color <- unique(plot_build$data[[1]]$fill)
    if(nameSim == "SW1-T2-D4"){
        color <- c(color[1:5], "darkgrey", color[6:8])
    }
    
    # graphical parameters
    dfYs$fill <- c(rep(color[bv$typeCell], 1), rep(color[bv$predPCA$typeCell],addBayVel))
    dfYs$color <- c(color[bv$typeCell], rep("darkgreen", length(bv$subtypeCell)*addBayVel))
    fillLegendYs <- c(rep("darkgreen", 1*addBayVel), "white")
   
   
    # -------------------------------------------------
    # create the plot for the positions
    ggYs <- ggplot(data = dfYs) + 
        geom_point(aes(x = pcI, y = pcII, shape = type, fill = fill, color = ifelse(type == "Ys", fill, color)), alpha = ifelse(dfYs$type == "Ys", 0.5, 0.7),  size = ifelse(dfYs$type == "Ys", 1, 3)) + 
        # Internal colors
        scale_fill_manual(values = sort(unique(dfYs$fill)), labels = unique(bv$typeCellReal)[order(unique(dfYs$fill))]) +
        # External colors
        scale_color_manual(values = sort(unique(dfYs$color))) +
        # Shape
        scale_shape_manual(values = c("Ys" = 21, "Real" = 22, "BayVel" = 23)) +
        # Legend
        guides(
                shape = guide_legend(override.aes = list(fill = fillLegendYs, size = 7)),
                fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 7)),
                color = "none",  
            )    
        
    ggYs <- ggYs +
        # Titles for legend
        labs(x = "PCA 1", y = paste0("PCA ", pc),  fill = "Groups", shape = "Methods", color = "Color Legend") + 
        # graphical parameters
        theme(
            plot.title = element_text(family = "serif", size=50,  hjust = 0.5), 
            axis.text = element_text(family = "serif", size = 35), 
            axis.title = element_text(family = "serif", size = 45), 
            legend.text = element_text(family = "serif", size = 38), legend.title = element_text(family = "serif", size = 40), legend.position = c(0.31,0.33), 
            legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), 
            plot.margin = unit(c(0.2,2.5,0,0.2), "lines")
        )

    pdf(paste0(pathOutput, "/Fig7_", ifelse(nameSim == "SW1-T1-D4", "a", "c"), ".pdf"), width = 7*1.2, height = 7 * 1.5)
    # Print the plot
        print(ggYs)
    dev.off()

    # ------ PLOT FOR VELOCITIES
    dfYs_velocity <- data.frame(
        pcI = NA, 
        pcII = NA, 
        pcI_prime_norm = NA, 
        pcII_prime_norm = NA, 
        color = c(as.factor(bv$typeCell), rep("darkgreen", max(bv$subtypeCell)*addBayVel)),
        fill = c(as.factor(bv$typeCell), rep(bv$predPCA$typeCell[mapping], addBayVel)), 
        type = c(rep("Ys", bv$n_cells), rep("BayVel", max(bv$subtypeCell)*addBayVel)), 
        typeCell = c(as.factor(bv$typeCell), rep(bv$predPCA$typeCell[mapping], addBayVel)), 
        subtypeCell = c(as.factor(bv$subtypeCell), rep(bv$predPCA$subtypeCell[mapping], addBayVel))
    )

    dfYs_velocity[which(dfYs_velocity$type == "Ys"), c("pcI", "pcII")] <- bv$pcaYs_df[, c("PC1", paste0("PC", pc))]

    dfYs_velocity[which(dfYs_velocity$type == "BayVel"), c("pcI", "pcII") ] <- bv$applicationPoint[, c("PC1", paste0("PC", pc))]
    dfYs_velocity[which(dfYs_velocity$type == "BayVel"), c("pcI_prime_norm", "pcII_prime_norm")] <- bv$applicationPoint[, c("PC1norm", paste0("PC", pc, "norm"))]
    
    bv$applicationPoint$color <- color[bv$typeCell[mapping]]   
    
    # graphical parameters
    dfYs_velocity$fill <- c(rep(color[bv$typeCell], 1), rep(color[bv$predPCA$typeCell[mapping]],addBayVel))
    dfYs_velocity$color <- c(color[bv$typeCell], rep("darkgreen", max(bv$subtypeCell)*addBayVel))
    fillLegendYs <- c(rep("darkgreen", 1*addBayVel), "white")

    # plot velocities
    ggYs_velocity  <-  ggplot(data = dfYs_velocity) + 
        geom_point(aes(x = pcI, y = pcII, shape = type, fill = fill, color = ifelse(type == "Ys", fill, color)), alpha = ifelse(dfYs_velocity$type == "Ys", 0.4, 1),  size = ifelse(dfYs_velocity$type == "Ys", 1, 6.5)) +
        # Internal colors
        scale_fill_manual(values = sort(unique(dfYs_velocity$fill)), labels = unique(bv$typeCellReal)[order(unique(dfYs_velocity$fill))]) +
        # External colors
        scale_color_manual(values = sort(unique(dfYs_velocity$color))) +
        # Shape
        scale_shape_manual(values = c("Ys" = 21, "Real" = 22, "BayVel" = 23)) +
        # Legend
        guides(
                shape = guide_legend(override.aes = list(fill = fillLegendYs, size = 7)),
                fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 7)),
                color = "none",  
                alpha = "none"
        )

    
    ggYs_velocity <- ggYs_velocity +
        # Titles for legend
        labs(x = "PCA 1", y = paste0("PCA ", pc), fill = "Groups", shape = "Methods", color = "Color Legend") + 
        # graphical parameters
        theme(
            plot.title = element_text(family = "serif", size=50,  hjust = 0.5), 
            axis.text = element_text(family = "serif", size = 35), 
            axis.title = element_text(family = "serif", size = 45), 
            legend.text = element_text(family = "serif", size = 38), legend.title = element_text(family = "serif", size = 40), legend.position = c(0.31,0.33), 
            legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), 
            plot.margin = unit(c(0.2,2.5,0,0.2), "lines")
        )

    # add arrows
    aspect_ratio <- 1   
    k <- 1.5
    bv$applicationPoint$pcI_prime_norm_corrected <- aspect_ratio * (bv$applicationPoint$pcI_prime_norm - bv$applicationPoint$pcI) + bv$applicationPoint$pcI


    for(sty in 1:bv$n_subtypeC){
        l <- 0.01
        if(bv$applicationPoint$pcI[sty] < bv$applicationPoint$pcI_prime_norm_corrected[sty]){
            while(bv$applicationPoint$pcI[sty] + l  < bv$applicationPoint$pcI_prime_norm_corrected[sty]){
                bv$applicationPoint$pcI_terminalArrow[sty] <- bv$applicationPoint$pcI[sty] + l
                l <- l + 0.01
            }
        }else{
            while(bv$applicationPoint$pcI[sty] - l  > bv$applicationPoint$pcI_prime_norm_corrected[sty]){
                bv$applicationPoint$pcI_terminalArrow[sty] <- bv$applicationPoint$pcI[sty] - l
                l <- l + 0.01
            }
        }
    }

    ggYs_velocity <- ggYs_velocity + 
    geom_segment(data = bv$applicationPoint, aes(x = pcI, y = pcII, xend = pcI_prime_norm_corrected, yend = pcII_prime_norm), arrow = arrow(length = unit(0.5, "cm")), linewidth = 1*k, color = "black") + 
    geom_segment(data = bv$applicationPoint, aes(x = pcI_terminalArrow , y = ((pcI_terminalArrow  - pcI_prime_norm_corrected)/(pcI_prime_norm_corrected - pcI)*(pcII_prime_norm - pcII) + pcII_prime_norm), xend = pcI_prime_norm_corrected, yend = pcII_prime_norm), color = "black", arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.25*k) + 
    geom_segment(data = bv$applicationPoint, aes(x = pcI_terminalArrow , y = ((pcI_terminalArrow  - pcI_prime_norm_corrected)/(pcI_prime_norm_corrected - pcI)*(pcII_prime_norm - pcII) + pcII_prime_norm), xend = pcI_prime_norm_corrected, yend = pcII_prime_norm, color = color), arrow = arrow(length = unit(0.5, "cm")), linewidth = 0.65*k)

    pdf(paste0(pathOutput, "/Fig7_", ifelse(nameSim == "SW1-T1-D4", "b", "d"), ".pdf"), width = 7*1.2, height = 7 * 1.5)
        print(ggYs_velocity )
    dev.off()

}
