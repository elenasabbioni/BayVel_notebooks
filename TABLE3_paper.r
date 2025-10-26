# ---------------------------------------------
# Script to produce Table 3 of BayVel paper, 
# with the percentage of true parameters that fall within the 95% credible interval, as obtained from the posterior distributions estimated by BayVel on simulated data.
# ---------------------------------------------
#
# OUTPUT: .txt file with the LaTeX code to reproduce table 3
#
# INPUTS:
# - pathToYourDirectory: path to the working directory
# - pathToResults: path where the results of scVelo are
# - pathOutput: path where you want to save the output
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
# TABLE 3
# -----------------------------
maxLog <- FALSE
n_genes <- 2000
typeSIM <- "sim"
SW_vec <- c("SW1", "SW2")
T_vec <- c("T1", "T2", "T3")
g <- expand.grid(SW_vec, T_vec)
g$Var3 <- "D4"
combinations <- paste(g$Var1, g$Var2, g$Var3, sep = "-")
mcmc <- 250000


table3 <- data.frame(nameSim = combinations, u_SS_OFF = NA, s_SS_OFF = NA, u_SS_ON = NA, s_SS_ON = NA, u0_off = NA, s0_off = NA, pos_u = NA, pos_s = NA, vel = NA,  eta = NA, catt = NA)


for(nameSim in combinations){
    # environment where we will load the results from BayVel
    bv  <- new.env(parent = baseenv())
    # load scVelo results
    source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))

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

    for(j in colnames(table3)[-1]){
        print(j)

        if(j == "u_SS_ON"){
            chain <- bv$uSS_on_chain
            realPar <- real$alpha_real[, 2]/real$beta_real
        }else if(j == "u_SS_OFF"){
            chain <- bv$uSS_off_chain
            realPar <- real$alpha_real[, 1]/real$beta_real
        }else if(j =="s_SS_ON"){
            chain <- bv$sSS_on_chain
            realPar <- real$alpha_real[, 2]/real$gamma_real
        }else if(j =="s_SS_OFF"){
            chain <- bv$sSS_off_chain
            realPar <- real$alpha_real[, 1]/real$beta_real
        }else if(j == "pos_u"){
            chain <- bv$u_chain
            realPar <- real$pos_u_real
        }else if(j == "pos_s"){
            chain <- bv$s_chain
            realPar <- real$pos_s_real
        }else if(j == "vel"){
            chain <- bv$v
            realPar <- real$v_real
        }else if(j == "u0_off"){
            chain <- bv$u0_off_chain
            realPar <- real$u0_off_real
        }else if(j == "s0_off"){
            chain <- bv$s0_off_chain
            realPar <- real$s0_off_real
        }else if(j == "eta"){
            chain <- bv$Eta_chain
            realPar <- real$eta_real
        }else if(j == "catt"){
            chain <- bv$Catt_chain
            realPar <- real$catt_real
        }

        # q0.025
        if(length(dim(chain)) == 2){
            q1<- apply(chain, FUN = quantile, MARGIN = 1, 0.025)
            q2 <- apply(chain, FUN = quantile, MARGIN = 1, 0.975)
        }else if(length(dim(chain)) == 3){
            q1 <- as.vector(t(apply(chain, FUN = quantile, MARGIN = c(1, 2), 0.025)))
            q2 <- as.vector(t(apply(chain, FUN = quantile, MARGIN = c(1, 2), 0.975)))
            realPar <- as.vector(t(realPar))
        }

        # check if the real parameter is in the credible interval 
        inCI <- (round(realPar, 5) >= round(q1, 5)) & (round(realPar, 5) <= round(q2, 5))

        table3[which(table3$nameSim == nameSim),which(colnames(table3)== j)] <- mean(inCI)
    }
}

View(table3)

# save the code to generate the Table in TeX
sink(paste0(pathOutput, "/Tab3.txt"))
print(xtable(table3))
sink()

