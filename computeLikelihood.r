# -----------------------------------------------------------
# Script to compute the log-likelihood for BayVel results 
# for each MCMC iteration
# ----------------------------------------------
# OUTPUT: 
# - A .csv file with the log-likelihood of the model for each MCMC iteration

# INPUTS (set at the beginning):
# - pathToYourDirectory: path to the working directory
# - Type of results we are loading:
#    - n_genes: number of genes in the simulated data
#    - typeSim: type of data (simulated "sim", pancreas "Pancreas", dentate gyrus "Dentate Gyrus")
#    - typeSW: common ("SW1") vs. cluster-specific switching times ("SW2")
#    - typeT: number of subgroups ("T1", "T2", "T3")
#    - typeD: type of data distribution (Poisson without capture efficiency: 
#             "D1", Negative Binomial with capture efficiency: "D4")

# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: LaplacesDemon

# -----------------------------------------------------------


rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
# SIMULATION PARAMETERS
# -----------------------------
n_genes <- 2000                   # number of genes
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D4"                     # type of data distribution
typeSIM <- "sim"                  # type of data (simulated, Pancreas, DentateGyrus)
mcmc <- 250000                    # MCMC iterations
# -----------------------------
#  PATH 
# -----------------------------
# Set working directory ad load auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
chrTypeSIM <- ifelse(typeSIM == "sim", "simulations", "real data")
pathToResults <- paste0(pathToYourDirectory, "/", chrTypeSIM)
setwd(pathToBayVel_results)

source(paste0(pathToYourDirectory, "/functions.R"))

# -----------------------------
# OUTPUT PATH
# -----------------------------
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
pathOutput <- paste0(pathToYourDirectory, "/", chrTypeSIM, "/", nameSim, "/output/")

# -----------------------------
# PACKAGES
# -----------------------------
library(LaplacesDemon)


# -----------------------------
# LOAD THE DATA 
# -----------------------------
# --- BayVel results
bv <- new.env(parent = baseenv()) 
maxLog <- FALSE
source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))

# load file with auxiliary functions
source(paste0(pathToYourDirectory, "/functions.R"))


# -----------------------------
# COMPUTE LIKELIHOOD
# -----------------------------
likTOT <- rep(NA, bv$SampleToSave)
    
if(typeD == "D1"){ # Poisson
    for(iter in 1:bv$SampleToSave){
        bv$loglik <- dpois(array(bv$Y_u, dim = c(bv$n_cells, bv$n_genes)), 
                           lambda = bv$u_chain[bv$subtypeCell, , iter], 
                           log = TRUE) + 
                     dpois(array(bv$Y_s, dim =c(bv$n_cells, bv$n_genes)), 
                           lambda = bv$s_chain[bv$subtypeCell, , iter], 
                           log = TRUE)

        bv$likTOT[iter] <- sum(bv$loglik)
    }
}else{             # Negative Binomial
    for(iter in 1:bv$SampleToSave){
        bv$loglik <- dnbinom(bv$Y_u, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), 
                            mu = bv$u_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), 
                            log = TRUE) + 
                     dnbinom(bv$Y_s, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), 
                            mu = bv$s_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), 
                            log = TRUE)

        bv$likTOT[iter] <- sum(bv$loglik)
    }  
}

# -----------------------------
# SAVE RESULTS
# -----------------------------
write.csv(bv$likTOT, file = paste0(pathOutput, "/loglik_", typeSIM, "_", nameSim, "_", n_genes, "genes_", mcmc, ".csv"), row.names = FALSE)
    
