# ---------------------------------------------
# Script to produce Fig. 8 of BayVel paper with representation of subgroups on real data
# ----------------------------------------------
#
# OUTPUTS
# - 2 .pdf files with plots in Fig. 8c and 8d
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


pathToYourDirectory <- "C:/Users/elena/Dropbox (Politecnico Di Torino Studenti)/sabbioni/Code Elena/Reproducibility/BayVel_notebooks"
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
typeSW <- "SW1"
typeT <- "T1"
typeD <- "D4"
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
mcmc <- 250000


bv <- new.env(parent = baseenv()) 
maxLog <- TRUE
source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))   
source(paste0(pathToYourDirectory, "/functions.R"))

# -------------------------------
# PLOTS
# -------------------------------
# extract colors
df <- data.frame(s = bv$Y_s[,1], u = bv$Y_u[,2], typeCell = as.factor(bv$typeCell), subtypeCell = as.factor(bv$subtypeCell))
gg <- ggplot(data = df) +
  geom_point(aes(x = s, y = u, fill = subtypeCell, color = subtypeCell))
plot_build <- ggplot_build(gg)
color <- unique(plot_build$data[[1]]$fill)

# genes to plot
genesName <- c("Cpe", "Fkbp2")
geniNum <- which(bv$nameG %in% genesName)

for(g in geniNum){
  
  dfYs <- data.frame(
    x = bv$s_chain[, g, 1], 
    y = bv$u_chain[, g, 1],
    fill = as.factor(seq(1, max(bv$subtypeCell)))
  )
    
  gg <- ggplot(data = dfYs) +
      geom_point(aes(x = x, y = y, fill = fill))
  plot_build <- ggplot_build(gg)
  color <- unique(plot_build$data[[1]]$fill)
  dfYs$fill <- sort(color)

  for (tyT0_off in 1:max(bv$typeCellT0_off)){  
    #  BayVel dynamic without  positions
    gg <- plot_sVSu(t0_off = bv$T0_off_chain[tyT0_off, g, 1], 
                    t0_on = 0, 
                    alpha = bv$alpha_chain[g, , 1], 
                    beta = bv$beta_chain[g, 1], 
                    gamma = bv$gamma_chain[g, 1], 
                    pos_u = NA, 
                    pos_s = NA, 
                    g = g, 
                    subGrLabels = bv$subtypeCell, 
                    add = FALSE, 
                    colCell = NA,
                    colDyn = "darkgreen", 
                    lineSize = 1
                    )

    # add position of subgroups
    gg <- gg + geom_point(data = dfYs, aes(x = x, y = y, fill = fill), size = 8, color = "darkgreen", shape = 21) +  
    # add legend and graphical elements
      scale_fill_manual(values = sort(color), labels = c("Ductal", "Beta", "Pre-endocrine", "Epsilon", "Ngn3 low EP", "Ngn3 high EP", "Delta", "Alpha")) +
      labs(fill = "Groups", x = "s", y = "u", title = paste0(bv$nameG[g])) + 
      guides(fill = guide_legend(override.aes = list( size = 7))) + 
      theme(
        plot.title = element_text(family = "serif", size=50,  hjust = 0.5), 
        axis.text = element_text(family = "serif", size = 35), 
        axis.title = element_text(family = "serif", size = 45), 
        legend.text = element_text(family = "serif", size = 38), 
        legend.title = element_text(family = "serif", size = 40), 
        legend.position = c(0.65,0.35), 
        legend.key = element_rect(colour = "transparent", fill = NA), 
        legend.background=element_blank()
      ) 
      

    pdf(paste0(pathOutput, "/Fig8_", ifelse(bv$nameG[g] == "Cpe", "b", "c"), ".pdf"), width = 7*1.2, height = 7 * 1.5)
      print(gg)
    dev.off()
  }
}

