# ---------------------------------------------
# Script to produce Table 1 of BayVel paper, 
# with the median relative error for the parameters estimated by scVelo, with simulated data
# ---------------------------------------------
#
# OUTPUT: .txt file with the LaTeX code to reproduce table 1
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

# -----------------------------
# FUNCTIONS
# -----------------------------

# Function to compute the relative mean square error
relError <- function(estim, real, vettori = FALSE){
    propZero <- NA
    resZero <- NA
    if(!vettori){
        if(length(estim) == length(real)){
            include <- which(real != 0, arr.ind = TRUE)
            resPos <- median(abs((estim[include] - real[include])/real[include]))
        }else{
            error("Error in dimensions of the input objects")
        }
    }else{
        norm1 <- apply(estim - real_par, FUN = norm, MARGIN = 1, type = "2")
        norm2 <- apply(real_par, FUN = norm, MARGIN = 1, type = "2")
        resPos <- sum(norm1[which(norm2 > 0)]/norm2[which(norm2 > 0)])/length(which(norm2 > 0))
    }
    names(resPos) <- c("relPositivi")
    return(resPos)
}


# -----------------------------
# TABLE 1
# -----------------------------

n_genes <- 2000
typeSIM <- "sim"
SW_vec <- c("SW1")
T_vec <- c("T1", "T2", "T3")
g <- expand.grid(SW_vec, T_vec)
g$Var3 <- "D4"
combinations <- paste(g$Var1, g$Var2, g$Var3, sep = "-")

# output
table1 <- data.frame("u_SS_ON" = rep(NA, 6),
        "u_SS_OFF" = NA,
        "s_SS_ON" = NA,
        "s_SS_OFF" = NA,
        "u0_off" = NA, 
        "s0_off" = NA, 
        "pos_u" = NA,
        "pos_s" = NA, 
        "vel" = NA)
rownames(table1) <- c("IN, K = 1, R = 10", "IN, K = 1, R = 30", "IN, K = 1, R = 100", "Dem, K = 1, R = 10", "Dem, K = 1, R = 30", "Dem, K = 1, R = 100")

# iterate over the different models and compute the relative error for each setting
k <- 1
for(modelScVel in c("NormInd", "Demings")){
    for(nameSim in combinations){    
        real <- new.env(parent = baseenv())
        scv  <- new.env(parent = baseenv())

        # load scVelo results
        source(paste0(pathToYourDirectory, "/loadScVelo_results.r"))
        
        # load real results
        load(paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData"), envir = real)

        # upload again the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
        source(paste0(pathToYourDirectory, "/functions.R"))
    
        # compute the error for the different parameters of the model
        for(j in c("u_SS_ON", "u_SS_OFF", "s_SS_ON", "s_SS_OFF", "u0_off", "s0_off", "pos_u", "pos_s", "vel")){
            if(j =="s_SS_ON"){
                estim <- scv$fit_alpha/scv$fit_gamma + scv$fit_s0
                estim <- estim[1,]
                real_par <- real$alpha_real[2,]/real$gamma_real
            }else if(j =="s_SS_OFF"){
                estim <- scv$fit_s0
                estim <- estim[1,]
                real_par <- real$alpha_real[1,]/real$gamma_real
            }else if(j =="u_SS_ON"){
                estim <- scv$fit_alpha/scv$fit_beta + scv$fit_u0
                estim <- estim[1,]
                real_par <- real$alpha_real[2,]/real$beta_real
            }else if(j =="u_SS_OFF"){
                estim <- scv$fit_u0
                estim <- estim[1,]
                real_par <- real$alpha_real[1,]/real$beta_real
            }else if(j == "u0_off"){
                estim <- scv$u0_off[1,,1]
                real_par <- real$u0_off_real
            }else if(j == "s0_off"){
                estim <- scv$s0_off[1,,1]
                real_par <- real$s0_off_real
            }else if(j == "pos_u"){
                estim <- scv$pos_u
                real_par <- real$pos_u_real
            }else if(j == "pos_s"){
                estim <- scv$pos_u
                real_par <- real$pos_s_real
            }else if(j == "vel"){
                estim <- scv$velocity 
                real_par <- matrix(rep(real$beta_real, real$n_cells), nrow = real$n_cells, ncol = real$n_genes, byrow = TRUE)* real$pos_u_real -  matrix(rep(real$gamma_real, real$n_cells), nrow = real$n_cells, ncol = real$n_genes, byrow = TRUE)*real$pos_s_real
            }
            table1[k, j] <- relError(estim, real_par, FALSE)   
        }
        k <- k + 1
    }
}

    
# save the code to generate the Table in TeX
sink(paste0(pathOutput, "/Tab1.txt"))
print(xtable(table1))
sink()


