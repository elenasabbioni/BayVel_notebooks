# -----------------------------------------------------------
# Simulation script: generates synthetic single-cell RNA data
# with different gene dynamics, switching times, and data distribution, 
# according to scVelo model (Demings Residuals and Independent Normal distributions). 
# The parameters of the dynamic will be the same of the corresponding BayVel simulation.
#
# The user can decide the simulation settings modifying the 
# parameters (type of switching times, number of time clusters, 
# data distribution, etc.) at the beginning of the code. 
# The code produces simulated counts for unspliced (Y_u) and 
# spliced (Y_s) RNA molecules.
# 
# OUTPUT:
# - ..._Mu.csv and ..._Ms.csv files saved for the chosen simulation setting, containing the unspliced and spliced values according to Independent Normal distributions model
# - ..._DEMINGS_Mu.csv and ..._DEMINGS_Mu.csv files saved for the chosen simulation setting, containing the unspliced and spliced values according to Deming residuals model
#
# INPUTS (set at the beginning):
# - pathToYourDirectory: path to the working directory
# - n_cells: number of cells
# - n_genes: number of genes
# - n_typeC: number of groups
# - typeSW: common ("SW1") vs. cluster-specific switching times ("SW2")
# - typeT: number of subgroups ("T1", "T2", "T3")
# - typeD: type of data distribution (Poisson without capture efficiency: 
#          "D1", Negative Binomial with capture efficiency: "D4")
#
# DEPENDENCIES:
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())


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
# SIMULATION PARAMETERS
# -----------------------------
n_cells <- 3000                   # Number of cells
n_genes <- 2000                   # number of genes
n_typeC <- 10                     # number of types of cells
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D4"                     # type of data distribution

# -----------------------------
# LOAD BAYVEL SIMULATED DATA (TO USE THE SAME SIMULATED PARAMTERES)
# -----------------------------
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
load(file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData"))

# -----------------------------
# SIMULATE ACCORDING TO INDEPENDENT NORMAL DISTRIBUTIONS MODEL
# The simulation strategy follows scVelo's simulation procedure, implemented in scvelo.datasets.simulation
# -----------------------------
set.seed(seed)
# compute variance of the normal distributions (as in scvelo.datasets.simulation)
noise.level.u <- 0.8*quantile(pos_u_real, probs = 0.99)/10
noise.level.s <- 0.8*quantile(pos_s_real, probs = 0.99)/10
# Sample Mu and Ms according to Independent Normal distributions model
Mu <- matrix(rnorm(length(pos_u_real), mean = pos_u_real, sd = noise.level.u), nrow = dim(pos_u_real)[1])
Ms <- matrix(rnorm(length(pos_s_real), mean = pos_s_real, sd =  noise.level.s), nrow = dim(pos_s_real)[1])
Mu[which(Mu < 0)] <- 0
Ms[which(Ms < 0)] <- 0

# Save the data 
write.csv(Mu, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_NormInd_Mu.csv"), row.names = FALSE)
write.csv(Ms, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_NormInd_Ms.csv"), row.names = FALSE)

# -----------------------------
# SIMULATE ACCORDING TO DEMING RESIDUAL MODEL
# The simulation strategy follows scVelo's likelihood
# -----------------------------
set.seed(seed)
# compute variance of the normal distributions (as in scvelo.datasets.simulation)
noise.level.u <- 0.8*quantile(pos_u_real, probs = 0.99)/10
noise.level.s <- 0.8*quantile(pos_s_real, probs = 0.99)/10
# Sample Mu and Ms according to Deming residuals model
theta <- runif(length(pos_u_real), pi/2, 3*pi/2)
rho <- rnorm(length(pos_u_real), 0, mean(c(noise.level.u, noise.level.s)))
Mu <- rho*sin(theta) + pos_u_real
Ms <- rho*cos(theta) + pos_s_real
Mu[which(Mu < 0)] <- 0
Ms[which(Ms < 0)] <- 0

# Save the data
write.csv(Mu, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_Demings_Mu.csv"), row.names = FALSE)
write.csv(Ms, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_Demings_Ms.csv"), row.names = FALSE)



# --- SAVE AUXILIARY FILES
obs <- cbind(typeCell, subtypeCell)
colnames(obs) <- c("clusters_coarse", "clusters")
write.csv(obs, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_obs.csv"), row.names = FALSE)

var <- seq(1, n_genes)
var <- cbind(var, "NaN")
colnames(var) <- c("index", "highly_variable_genes")
write.csv(var, file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, "_var.csv"), row.names = FALSE)

