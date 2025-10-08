# ---------------------------------------------
# Script to produce Table 2 of BayVel paper, 
# with the median relative error for the parameters estimated by BayVel, with simulated data
# ---------------------------------------------
#
# OUTPUT: .txt file with the LaTeX code to reproduce table 2
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
# TABLE 2
# -----------------------------
maxLog <- FALSE
n_genes <- 2000
typeSIM <- "sim"
SW_vec <- c("SW1", "SW2")
T_vec <- c("T1", "T2", "T3")
g <- expand.grid(SW_vec, T_vec)
g$Var3 <- "D4"
combinations <- paste(g$Var1, g$Var2, g$Var3, sep = "-")
mcmc <- 150

# output
table2 <- data.frame("u_SS_ON" = rep(NA, 6),
        "u_SS_OFF" = NA,
        "s_SS_ON" = NA,
        "s_SS_OFF" = NA,
        "u0_off" = NA, 
        "s0_off" = NA, 
        "pos_u" = NA,
        "pos_s" = NA, 
        "vel" = NA, 
        "eta" = NA,
        "catt" = NA
        )
rownames(table2) <- c("K = 1, R = 10", "K = 1, R = 30", "K = 1, R = 100", "K = 10, R = 10", "K = 10, R = 30", "K = 10, R = 100")

# iterate over the different models and compute the relative error for each setting
k <- 1
for(nameSim in combinations){    
    real <- new.env(parent = baseenv())
    bv  <- new.env(parent = baseenv())

    # load BayVel results
    source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))

    # load real results
    load(paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData"), envir = real)

    # upload again the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
    source(paste0(pathToYourDirectory, "/functions.R"))

    # compute the error for the different parameters of the model
    for(j in c("u_SS_ON", "u_SS_OFF", "s_SS_ON", "s_SS_OFF", "u0_off", "s0_off", "pos_u", "pos_s", "vel", "eta", "catt")){
        if(j =="u_SS_ON"){
            chain <- bv$uSS_on_chain
            real_par <- bv$u_SSon_real
        }else if(j =="u_SS_OFF"){
            chain <- bv$uSS_off_chain
            real_par <- bv$u_SSoff_real
        }else if(j =="s_SS_ON"){
            chain <- bv$sSS_on_chain
            real_par <- bv$s_SSon_real
        }else if(j =="s_SS_OFF"){
            chain <- bv$sSS_off_chain
            real_par <- bv$s_SSoff_real
        }else if(j == "u0_off"){
            chain <- bv$u0_off_chain
            real_par <- bv$u0_off_real
        }else if(j == "s0_off"){
            chain <- bv$s0_off_chain
            real_par <- bv$s0_off_real
        }else if(j == "pos_u"){
            chain <- bv$u_chain
            real_par <- bv$pos_u_real
        }else if(j == "pos_s"){
            chain <- bv$s_chain
            real_par <- bv$pos_s_real
        }else if(j == "vel"){
            chain <- bv$v
            real_par <- matrix(rep(bv$beta_real, bv$n_cells), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE)* bv$pos_u_real -  matrix(rep(bv$gamma_real, bv$n_cells), nrow = bv$n_cells, ncol = bv$n_genes, byrow = TRUE)*bv$pos_s_real
        }else if(j == "eta"){
            chain <- bv$Eta_chain
            real_par <- bv$eta_real
        }else if(j == "catt"){
            chain <- bv$Catt_chain
            real_par <- bv$catt_real
        }

        if(length(dim(chain)) == 2){
            medianChain <- apply(chain, FUN = median, MARGIN = 1)
        }else if(length(dim(chain)) == 3){
            medianChain <- apply(chain, FUN = median, MARGIN = c(1, 2))
        }

        if(j %in% c("pos_u", "pos_s", "vel")){
            medianChain <- medianChain[bv$subtypeCell,]
        }

        table2[k, j] <- relError(medianChain, real_par, FALSE)   
        k <- k + 1
    }
}

    
# save the code to generate the Table in TeX
sink(paste0(pathOutput, "/Tab2.txt"))
print(xtable(table2))
sink()


