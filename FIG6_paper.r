# ---------------------------------------------
# Script to produce Fig. 6 of BayVel paper with representation of subgroups on real data
# ----------------------------------------------
#
# OUTPUTS
# - 3 .pdf files with plots in Fig. 6
#
# INPUTS (set at the beginning):
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
#             LaplacesDemon
#             latex2exp
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
library(LaplacesDemon)
library(ggplot2)
library(latex2exp)

# -----------------------------
# TYPE OF SIMULATION WE ARE CONSIDERING
# -----------------------------
n_genes <- 2000                   # number of genes
typeSIM <- "Pancreas"             # specify that we are considering simulated data
addBayVel <- TRUE                 # specify if you want to plot the results of BayVel

# ------------------------------
# LOAD SCVELO DATA
# ------------------------------
typeSW <- "SW1"
typeT_vec <- c("T1", "T2", "T3")
typeD <- "D4"
nameSim_vec <- paste0(typeSW, "-", typeT_vec, "-", typeD)
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

    pdf(paste0(pathOutput, "/Fig6_", ifelse(nameSim == "SW1-T1-D4", "a", ifelse(nameSim == "SW1-T2-D4", "b", "c")), ".pdf"), width = 7*1.2, height = 7 * 1.5)
    
    pc <- 2
    # create an unique df
    dfYs <- data.frame(
        pcI = NA, 
        pcII = NA, 
        pcI_prime_norm = NA, 
        pcII_prime_norm = NA, 
        color = c(as.factor(bv$typeCell)),
        fill = c(as.factor(bv$typeCell)), 
        type = c(rep("Ys", bv$n_cells)), 
        typeCell = c(as.factor(bv$typeCell)), 
        subtypeCell = c(as.factor(bv$subtypeCell))
    )
    dfYs[which(dfYs$type == "Ys"), c("pcI", "pcII")] <- bv$pcaYs_df[, c("PC1", paste0("PC", pc))]
    
    # select better colors
    gg <- ggplot(data = subset(dfYs, type == "Ys")) +
    geom_point(aes(x = pcI, y = pcII, fill = subtypeCell, color = subtypeCell))
    plot_build <- ggplot_build(gg)
    color <- unique(plot_build$data[[1]]$fill)
    if(nameSim == "SW1-T2-D4"){
        color <- c(color[1:5], "darkgrey", color[6:8])
    }
    
    # graphical parameters
    dfYs$fill <- c(rep(color[bv$typeCell], 1))
    dfYs$color <- c(color[bv$subtypeCell])
    fillLegendYs <- c("white")
   
    # -------------------------------------------------
    # create the plot for the positions
    dfYs_Ys <- subset(dfYs, type == "Ys")

    ggYs <- ggplot(data = dfYs_Ys) + 
        geom_point(aes(x = pcI, y = pcII, shape = type, fill = color, color = color), alpha = 1,  size = 1) +
        # Internal colors
        scale_fill_manual(values = sort(unique(dfYs_Ys$color)), labels = unique(bv$typeCellReal)[order(unique(dfYs_Ys$fill))]) +
        # External colors
        scale_color_manual(values = sort(unique(dfYs_Ys$color))) +
        # Shape
        scale_shape_manual(values = c("Ys" = 21, "Real" = 22, "BayVel" = 23)) +
        # Legend
        guides(
            shape = "none",
            color = "none",  
            alpha = "none"
        )    

    # Modify legends depending on the type of panel
    if(nameSim == "SW1-T1-D4"){
        ggYs <- ggYs + 
                guides(
                    fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 7)),
                )    
    }else{
        ggYs <- ggYs + 
                guides(
                    fill = "none",
                )
    }
    
        
    ggYs <- ggYs +
        # Titles for legend
        labs(x = "PCA 1", y = paste0("PCA ", pc), title = TeX(paste0("R = ", length(unique(bv$subtypeCell)))), fill = "Groups") + 
        # graphical parameters
        theme(
            plot.title = element_text(family = "serif", size=50,  hjust = 0.5), 
            axis.text = element_text(family = "serif", size = 35), 
            axis.title = element_text(family = "serif", size = 45), 
            legend.text = element_text(family = "serif", size = 38), legend.title = element_text(family = "serif", size = 40), legend.position = c(0.31,0.3), 
            legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), 
            plot.margin = unit(c(0.2,2.5,0,0.2), "lines")
        )

    # Print the plot
    print(ggYs)

    dev.off()
}


