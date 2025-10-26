# ---------------------------------------------
# Script to produce Table 4 of BayVel paper, 
# with the median absolute error estimated by BayVel and scVelo on simulated data.
# ---------------------------------------------
#
# OUTPUT: .txt file with the LaTeX code to reproduce table 4
#
# INPUTS:
# - pathToYourDirectory: path to the working directory
# - pathToResults: path where the results of scVelo are
# - pathOutput: path where youu wanto to save the output
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: xtable
#             data.table
# ---------------------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set working directory and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
pathToResults <- paste0(pathToYourDirectory, "/simulations")
setwd(paste0(pathToYourDirectory))
source(paste0(pathToYourDirectory, "/functions.R"))

# Set the output path where you will save the table
pathOutput <- paste0(pathToYourDirectory, "/tablesPaper/")

# -----------------------------
#  PACKAGES 
# -----------------------------
library(xtable)
library(data.table)
library(LaplacesDemon)

# -----------------------------
# TABLE 4
# -----------------------------
maxLog <- FALSE
n_genes <- 2000
typeSIM <- "sim"
typeSW <- "SW1"
typeT <- "T1"
typeD <- "D4"
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
modelScVel <- "Demings"
mcmc <- 250000

table4 <- data.frame(nameMethod = c("BayVel", "scVelo"), u_SS_OFF = NA, s_SS_OFF = NA, u_SS_ON = NA, s_SS_ON = NA, vel = NA)

# environment where we will load the results from BayVel
bv  <- new.env(parent = baseenv())
# load BayVel results
source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))


# environment where we will load the results from scVelo
scv  <- new.env(parent = baseenv())
# load scVelo results
source(paste0(pathToYourDirectory, "/loadScVelo_results.r"))

# load real results
real <- new.env(parent = baseenv())
load(paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData"), envir = real)

# upload again the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
source(paste0(pathToYourDirectory, "/functions.R"))


# uniform the dimension of the parameters in real and BayVel environment
real$new_pos_u_real <- matrix(NA, nrow = real$n_subtypeC, ncol = real$n_genes)
real$new_pos_s_real <- matrix(NA, nrow = real$n_subtypeC, ncol = real$n_genes)

for(sty in 1:real$n_subtypeC){
    subcellTy <- which(real$subtypeCell == sty)
    real$new_pos_u_real[sty, ] <- real$pos_u_real[subcellTy[1], ]
    real$new_pos_s_real[sty, ] <- real$pos_s_real[subcellTy[1], ]
}
real$pos_u_real <- real$new_pos_u_real
real$pos_s_real <- real$new_pos_s_real

real$v_real <- matrix(NA, nrow = real$n_subtypeC, ncol = real$n_genes)
for(sty in 1:real$n_subtypeC){
    real$v_real[sty, ] = real$beta_real*real$pos_u_real[sty,] - real$gamma_real*real$pos_s_real[sty,]
}   


for(j in colnames(table4)[-1]){
    print(j)

    if(j == "u_SS_ON"){
        chain <- bv$uSS_on_chain
        realPar <- real$alpha_real[, 2]/real$beta_real

        estimScVelo <- scv$fit_s0
        estimScVelo <- estimScVelo[1,]
    }else if(j == "u_SS_OFF"){
        chain <- bv$uSS_off_chain
        realPar <- real$alpha_real[, 1]/real$beta_real

        estimScVelo <- scv$fit_u0
        estimScVelo <- estimScVelo[1,]
    }else if(j == "s_SS_ON"){
        chain <- bv$sSS_on_chain
        realPar <- real$alpha_real[, 2]/real$gamma_real
        
        estimScVelo <- scv$fit_alpha/scv$fit_gamma + scv$fit_s0
        estimScVelo <- estimScVelo[1,]
    
    }else if(j == "s_SS_OFF"){
        chain <- bv$sSS_off_chain
        realPar <- real$alpha_real[, 1]/real$beta_real

        estimScVelo <- scv$fit_alpha/scv$fit_beta + scv$fit_u0
        estimScVelo <- estimScVelo[1,]
    }else if(j == "vel"){
        chain <- bv$v
        realPar <- real$v_real

        estimScVelo <- scv$velocity # computed using fit_beta*fit_scaling 
    }

    if(length(dim(chain)) == 2){
        medianChain <- apply(chain, FUN = median, MARGIN = 1)
    }else if(length(dim(chain)) == 3){
        medianChain <- apply(chain, FUN = median, MARGIN = c(1, 2))
    }

    if(j %in% c("pos_u", "pos_s", "vel")){
        medianChain <- medianChain[bv$subtypeCell,]
    }
        
    include <- which(realPar != 0, arr.ind = TRUE)
    resBayVel <- median(abs(medianChain[include] - realPar[include]))
    resScVelo <- median(abs(estimScVelo[include] - realPar[include]))

    table4[which(table4$nameMethod == "BayVel"), j] <- resBayVel
    table4[which(table4$nameMethod == "scVelo"), j]  <- resScVelo
}


# save the code to generate the Table in TeX
sink(paste0(pathOutput, "/Tab4.txt"))
print(xtable(table4))
sink()