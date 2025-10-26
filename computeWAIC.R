# -----------------------------------------------------------
# Script to compute Watanabe–Akaike information criterion (WAIC) for different BayVel models. 
# WAIC is defined as 
# WAIC = -2*(
#               \sum_{c, g}log(1/MC*\sum_{i = 1}^{MC} f(Y_{u, cg}, Y_{s, cg}| \theta_i)) + 
#            -  \sum_{c, g} Var_{\theta}(log(f(Y_{u, cg}, Y_{s, cg}| \theta)))
#         )
# 
# OUTPUTS: .csv file with WAIC value for the specific simulation setting we are considering
#
# INPUT: 
# - pathToYourDirectory: path to the working directory
# - type of simulation parameters:
#    - typeSW: common ("SW1") vs. cluster-specific switching times ("SW2")
#    - typeT: number of subgroups ("T1", "T2", "T3")
#    - typeD: type of data distribution (Poisson without capture efficiency:
#          "D1", Negative Binomial with capture efficiency: "D4")
# - typeSIM: type of data we are considering (simulated "sim", Pancreas or DentateGyrus)
# - n_genes: total number of simulated genes (2000 by default)
# - mcmc: number of MCMC iterations
#
# DEPENDENCIES:
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: compiler
#             LaplacesDemon
#             
# -----------------------------------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
# PACKAGES
# -----------------------------
library(compiler)
library(LaplacesDemon)

# -----------------------------
# SIMULATION PARAMETERS
# -----------------------------
n_genes <- 2000                   # number of genes
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D4"                     # type of data distribution
typeSIM <- "Pancreas"             # type of data (Pancreas, DentateGyrus)
mcmc <- 250000                    # MCMC iterations

nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
# -----------------------------
#  PATH 
# -----------------------------
# Set working directory ad load auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
chrTypeSIM <- ifelse(typeSIM == "sim", "simulations", "real data")
pathToResults <- paste0(pathToYourDirectory, "/", chrTypeSIM)
setwd(pathToResults)
pathOutput <- paste0(pathToResults, "/", nameSim, "/output/")

source(paste0(pathToYourDirectory, "/functions.R"))

# -----------------------------
#  LOAD THE DATA 
# -----------------------------
# environment where we will load the results from the different methods
bv <- new.env(parent = baseenv()) 
source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))

# upload the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
source(paste0(pathToYourDirectory, "/functions.R"))

   
# -----------------------------
#  COMPUTE WAIC 
# -----------------------------
# function o compute the WAIC
# this function will be previously compiled, thanks to compiler package, and then run, to speed up the cde
computeWAIC <- function(bv, typeD){
    bv$maxLoglik <- matrix(-Inf, nrow = bv$n_cells, ncol = bv$n_genes)            
    # compute, for each cell and for each gene, the iteration with the maximum log-likelihood (to use log-sum-exp trick)
    if(typeD == "D1"){
        for(iter in 1:bv$SampleToSave){
            bv$loglik <- dpois(array(bv$Y_u, dim = c(bv$n_cells, bv$n_genes)), lambda = bv$u_chain[bv$subtypeCell, , iter], log = TRUE) + dpois(array(bv$Y_s, dim =c(bv$n_cells, bv$n_genes)), lambda = bv$s_chain[bv$subtypeCell, , iter], log = TRUE)
            bv$maxLoglik <- pmax(bv$maxLoglik,  bv$loglik)
        }
    }else{
        for(iter in 1:bv$SampleToSave){
            bv$loglik <-  dnbinom(bv$Y_u, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), mu = bv$u_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), log = TRUE) + dnbinom(bv$Y_s, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), mu = bv$s_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), log = TRUE)

            bv$maxLoglik <- pmax(bv$maxLoglik,  bv$loglik)
        }  
    }

    # compute likelihood and loglikehood for each cell,gene and iteration
    bv$likSAMPLE <- matrix(0, nrow = bv$n_cells, ncol = bv$n_genes)
    bv$loglikSAMPLE <- matrix(0, nrow = bv$n_cells, ncol = bv$n_genes)
    pWAIC <- 0
    if(typeD == "D1"){
        for(iter in 1:bv$SampleToSave){
            bv$loglik <- dpois(array(bv$Y_u, dim = c(bv$n_cells, bv$n_genes)), lambda = bv$u_chain[bv$subtypeCell, , iter], log = TRUE) + dpois(array(bv$Y_s, dim =c(bv$n_cells, bv$n_genes)), lambda = bv$s_chain[bv$subtypeCell, , iter], log = TRUE)
            
            bv$likSAMPLE <- bv$likSAMPLE + exp(bv$loglik - bv$maxLoglik)
            bv$loglikSAMPLE <- bv$loglikSAMPLE + bv$loglik
            pWAIC <- pWAIC + bv$loglik^2
        }
    }else{
        for(iter in 1:bv$SampleToSave){
            bv$loglik <-  dnbinom(bv$Y_u, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), mu = bv$u_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), log = TRUE) + dnbinom(bv$Y_s, size = matrix(t(rep(1/bv$Eta_chain[, iter], bv$n_cells)), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE), mu = bv$s_chain[bv$subtypeCell, , iter]*array(rep(bv$Catt_chain[, iter], bv$n_genes), dim = c(bv$n_cells, bv$n_genes)), log = TRUE)

            bv$likSAMPLE <- bv$likSAMPLE + exp(bv$loglik - bv$maxLoglik)
            bv$loglikSAMPLE <- bv$loglikSAMPLE + bv$loglik
            pWAIC <- pWAIC + bv$loglik^2
        }  
    }
    bv$loglikSAMPLE <- bv$loglikSAMPLE/bv$SampleToSave
    pWAIC <- (pWAIC - bv$SampleToSave*bv$loglikSAMPLE^2)/(bv$SampleToSave - 1)
    pWAIC <- sum(pWAIC)
    bv$lppd <- sum(-log(bv$SampleToSave) + bv$maxLoglik + log(bv$likSAMPLE)) # log-sum-exp trick 
    WAIC <- -2*(bv$lppd - pWAIC)

    return(WAIC)
}

# compile the function
compiledWAIC <- cmpfun(computeWAIC)
# compute WAIC
WAIC <- compiledWAIC(bv, typeD)

# save the results       
fwrite(WAIC, file = paste0(pathOutput, "/WAIC.csv"))

  
