# ---------------------------------------------
# Script to produce Fig. 3 of BayVel paper with comparison
# between the estimated gene-dynamics between scVelo and BayVel, on simulated data
# ----------------------------------------------
#
# OUTPUTS
# - 3 PDF files, one for each panel in Figure 3
#
# INPUTS (set at the beginning or controlled in loops):
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
pathToResults <- "pathToYourDirectory/simulations"
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
mcmc <- 150                       # number of iterations of the mcmc we have run
typeSIM <- "sim"                  # specify that we are considering simulated data
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D1"                     # type of data distribution
modelScVel <- "Demings"           # type of simulation method for scVelo results
genesToPlot <- c(1, 3, 12)        # indexes of the genes that we want to plot        

addBayVel <- TRUE                 # specify if you want to plot the results of BayVel
addScVelo <- TRUE                 # specify if you want to plot the results of scVelo


# -----------------------------
# LOAD THE DATA 
# -----------------------------
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)

# environment where we will load the results from the different methods
bv <- new.env(parent = baseenv()) 
scv <- new.env(parent = baseenv())

# --- BayVel results
if(addBayVel){
  maxLog <- FALSE
  source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))
  iterToPlot <- bv$loglikTOT[which.max(bv$loglikTOT[, 2]), 1]
}

# --- real parameters
if(typeSIM == "sim" & addReal){  
  real <- new.env(parent = baseenv())
  pathReal <- paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData")
}   

# --- scVelo results
if(addScVelo){
  source(paste0(pathToYourDirectory, "/loadScVelo_results.r"))
}    

# upload the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
source(paste0(pathToYourDirectory, "/functions.R"))

# 
list_plot <- list()
k <- 1
for(g in genesTo_plot){
  for (tyT0_off in 1:max(bv$typeCellT0_off)){  
    # real
    gg <- plot_sVSu(t0_off_real = real$t0_off_real[tyT0_off, g], t0_on_real = 0, alpha_real = real$alpha_real[g,], beta_real = real$beta_real[g], gamma_real = real$gamma_real[g], pos_u_real = NA, pos_s_real = NA, g = g, tipoCellula = bv$subtypeCell, add = FALSE, colCell = NA, , xlim = NA, ylim = NA, axisTitle.size = 20, axisText.size = 10, title.size = 20, colReal = "red", lineSize = 1, gg = NA)      
      
    # scVelo
    gg <- plot_sVSu_scVelo(fit_t0_off = scv$fit_t_[g], fit_t0_on = 0, fit_alpha = scv$fit_alpha[g], fit_beta = scv$fit_beta[g], fit_gamma = scv$fit_gamma[g], fit_scaling = scv$fit_scaling[g], fit_u0_offset = scv$fit_u0[g], fit_s0_offset = scv$fit_s0[g], fit_t = NA, g = g, tipoCellula = bv$subtypeCell, add = TRUE, colCell = NA, colReal = "blue", lineSize = 1, gg = gg)

    # BayVel
    gg <- plot_sVSu(t0_off_real = bv$T0_off_chain[tyT0_off, g, iterToPlot], t0_on_real = 0, alpha_real = bv$alpha_chain[g, , iterToPlot], beta_real = bv$beta_chain[g, iterToPlot], gamma_real = bv$gamma_chain[g, iterToPlot], pos_u_real = NA, pos_s_real = NA, g = g, tipoCellula = bv$subtypeCell, add = TRUE, colCell = NA, , xlim = 10, ylim = 10, axisTitle.size = 20, axisText.size = 20, title.size = 40, colReal = "darkgreen", lineSize = 1, gg = gg)

    # add legend
    gg <- gg + 
          geom_point(data = data.frame(x = bv$alpha_chain[g, 1, iterToPlot]/bv$gamma_chain[g, iterToPlot], y = bv$alpha_chain[g, 1, iterToPlot]/bv$beta_chain[g, iterToPlot], group = c("Real", "scVelo", "BayVel")), aes(x = x, y = y, color = group), size = 0.1) + 
          scale_color_manual(values = c("Real" = "red", "scVelo" = "blue", "BayVel" = "darkgreen")) + 
          labs(color = "") + 
          guides(color = guide_legend(override.aes = list(size = 15, shape = 15))) + 
          theme(legend.text = element_text("serif", size = 38),  legend.title = element_text(family = "serif", size = 40), legend.position="bottom", legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), plot.title = element_text(family = "serif", size=50,  hjust = 0.5), axis.text = element_text(family = "serif", size = 35),  axis.title = element_text(family = "serif", size = 45), plot.margin = unit(c(0.2,0.5,0,0), "lines")) 
 
    pdf(paste0(pathOutput, "Fig3_gene", g, "_tyT0_off", tyT0, ".pdf"), width = 7*1.2, height = 7*1.5)
      plot(gg)
    dev.off()
  }
}
