# -----------------------------------
# Plot simulated data: generates plot with the real gene dynamic for a specific simulation settings and selected genes.

#
# INPUTS (set at the beginning or controlled in loops):
# - pathToYourDirectory: path to the working directory
# - type of simulation parameters:
#    - typeSW: common ("SW1") vs. cluster-specific switching times ("SW2")
#    - typeT: number of subgroups ("T1", "T2", "T3")
#    - typeD: type of data distribution (Poisson without capture efficiency: 
#          "D1", Negative Binomial with capture efficiency: "D4")
# - n_genes: total number of simulated genes (2000 by default)
# - genesToPlot: vector with indexes of the genes to plot
# - plotPosition: boolean, to decide if yu want or not to add the subgroup positions in the plot

# DEPENDENCIES:
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
#             RcolorBrewer
# -----------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set working directory and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
setwd(paste0(pathToYourDirectory, "/simulations"))
source(paste0(pathToYourDirectory, "/functions.R"))


# -----------------------------
#  PACKAGES 
# -----------------------------
library(ggplot2)
library(RColorBrewer)

# -----------------------------
# SIMULATION PARAMETERS
# -----------------------------
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D1"                     # type of data distribution
n_genes <- 2000                   # total number of simulated genes 
genesToPlot <- seq(1, 3)          # indexes of the genes that we want to plot          
plotPosition <- TRUE              # boolean to decide if you want to plot or not the subgroup position in the plot

# -----------------------------
# OUTPUT PATH
# -----------------------------
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
pathOutput <- paste0(pathToYourDirectory, "/simulations/", nameSim, "/")

# -----------------------------
# LOAD SIMULATED DATA
# -----------------------------
load(file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/",nameSim, "_", n_genes, ".RData"))
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)                   


# -----------------------------
# PLOT GENE DYNAMIC 
# -----------------------------
chrPos <- ifelse(plotPosition, "posSubGr", "")
pdf(file = paste(pathOutput, "realDynamic_", chrPos ,"_", nameSim, "_", n_genes, ".pdf", sep = ""))

list_plot <- list()
k <- 1
for(g in genesToPlot){
  for (tyT0_off in 1:max(typeCellT0_off)){ # iterate over the different switching times (we plot each switching dynamic separately, otherwise the plot is too messy)

    cellTyT0_off <- which(typeCellT0_off == tyT0_off)
    subTy <- unique(subtypeCell[cellTyT0_off])
    if(plotPosition){
        # cells associated with the selected switching time
        plotU <- unique(pos_u_real[cellTyT0_off, g])
        plotS <- unique(pos_s_real[cellTyT0_off, g])
    }else{
        plotU <- NA
        plotS <- NA
    }
    
    # gene dynamic for a single switching time, without subgroup-specific positions
    gg_dynamic <- plot_sVSu(t0_off = t0_off_real[tyT0_off, g], t0_on = 0, alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g], pos_u = plotU, pos_s = plotS, g = g, subGrLabels = subTy, add = FALSE, colCell = NA, , xlim = NA, ylim = NA, axisTitle.size = 40, axisText.size = 20, title.size = 40, colDyn = "black", sizePoint = 10, lineSize = 1, gg = NA)      

    list_plot[[k]] <- gg_dynamic   
    k <- k + 1
  }
}

# plot the different dynamics
for(k in seq(1, length(list_plot))){
  plot(list_plot[[k]])
} 
dev.off()
