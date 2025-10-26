# -----------------------------------------------------------
# Simulation script: generates synthetic single-cell RNA data
# with different gene dynamics, switching times, and data distribution, 
# according to BayVel model.
#
# The user can decide the simulation settings modifying the 
# parameters (type of switching times, number of time clusters, 
# data distribution, etc.) at the beginning of the code.
# The code produces simulated counts for unspliced (Y_u) and 
# spliced (Y_s) RNA molecules.
# 
# OUTPUT:
# - An .RData file saved for the chosen simulation setting, containing:
#   - the unspliced and spliced counts Y_u, Y_s
#   - the parameters of the ODE alpha_real, beta_real, gamma_real, k_real, t_real, ...
#   - the parameters of the data distribution: catt_real, eta_real, ...
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
setwd(paste0(pathToYourDirectory))
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
# OUTPUT PATH
# -----------------------------
path <- paste0(pathToYourDirectory, "/simulations")
dir.create(path, showWarnings = FALSE)
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
path <- paste0(path, "/", nameSim)
dir.create(path, showWarnings = FALSE)
path <- paste0(path, "/")


# -----------------------------
# GROUP LABELS
# -----------------------------
# Assign cells to discrete groups: divide cells in groups of the same dimensions
typeCell <- rep(NA, n_cells)      
start <- 0
step <- floor(n_cells/n_typeC)
for(i in 1:n_typeC){
    end <- start + step
    typeCell[(start+1):end] <- i
    start <- end 
}                              

# -----------------------------
# SUBGROUP LABELS
# -----------------------------
if(typeT == "T1"){
    n_clusters <- 1
    sdTimes <- 0
}else if(typeT == "T2"){
    n_clusters <- 3
    sdTimes <- sqrt(0.5)
}else if(typeT == "T3"){
    n_clusters <- 10
    sdTimes <- sqrt(0.5)
}

# Assign cells to subgroups
subtypeCell <- rep(NA, n_cells)
for(ty in 1:n_typeC){
    cellTy <- which(typeCell == ty)
    for(i in 1:n_clusters){
        subtypeCell[cellTy[((i-1)*floor(length(cellTy)/n_clusters) + 1):(floor(i*length(cellTy)/n_clusters))]] <- i + n_clusters*(ty - 1)
    }        
}   
n_subtypeC <- max(subtypeCell)

# -----------------------------
# SWITCHING CLUSTERS' LABELS
# -----------------------------
if(typeSW == "SW1"){
    t0_offClustReal <- FALSE                # all cells share the same switching times
    typeCellT0_off <- rep(1, n_cells)
}else if(typeSW == "SW2"){
    t0_offClustReal <- TRUE                 # Each cluster has its own switching times
    typeCellT0_off <- typeCell
}
n_typeT0_off <- max(typeCellT0_off)

# Grid of possible switching times, distributed across the on branch
possU0 <- c()
possU0[1] <- 0.3
possT0 <- c()
possT0[1] <- -log(1-possU0[1])
stp = possT0[1] 
for(i in 2:20){
    possT0[i] = possT0[i-1]+stp
    possU0[i] = 1-exp(-possT0[i])
}
possT0 <- possT0[which(possU0 < 0.995)] # exclude switching points that are too close to the upper steady state

# -----------------------------------------------------
# INITIALIZE PARAMETER VECTORS
# (alphas, betas, gammas, steady states, switching times, etc.)
# -----------------------------------------------------
# rate of the ODE
beta_real <- rep(NA, n_genes)       
alpha_real <- rep(NA, 2*n_genes)    # alpha_off and alpha_on
attr(alpha_real, "dim") <- c(n_genes, 2)
gamma_real <- rep(NA, n_genes)      # gamma reali per generare dati simulati

# lower steady state
u_SSoff_real <- rep(NA, n_genes)    # u_steady_state_OFF reale
s_SSoff_real <- rep(NA, n_genes)    # s_steady_state_OFF reale  
# difference between u-coordinates of the steady states
diffSS_off_real <- rep(NA, n_genes) # diff u_steady_state_ON e u_steady_state_OFF reale  

# off-switching time and associated coordinates
t0_off_real <- matrix(NA, n_typeT0_off, n_genes)     # tempo di switch tra il ramo on e il ramo off
u0_off_real <- matrix(NA, n_typeT0_off, n_genes)     # u(t0_off): punto iniziale del ramo off
s0_off_real <- matrix(NA, n_typeT0_off, n_genes)     # s(t0_off): punto iniziale del ramo off

# on-switching time (assumed to be 0 for identifiability) and associated coordinates 
t0_on_real <- rep(0, n_genes)      
u0_on_real  <- rep(NA, n_genes)
s0_on_real  <- rep(NA, n_genes)

# elapsed time in the dynamic
t_real <- rep(NA, n_cells*n_genes)  
attr(t_real, "dim") <- c(n_cells, n_genes)
# elapsed time since the last switching point
tau_real <- rep(NA, n_cells*n_genes)           # tempo medio di ogni tipo di cellula
attr(tau_real, "dim") <- c(n_cells, n_genes)
# associated coordinates in the dynamic
pos_u_real <- rep(NA, n_cells*n_genes)
attr(pos_u_real, "dim") <- c(n_cells, n_genes)
pos_s_real <- rep(NA, n_cells*n_genes)
attr(pos_s_real, "dim") <- c(n_cells, n_genes)
# branch in the dynamic (0: off, 2: on)
k_real <- rep(NA, n_cells*n_genes)       
attr(k_real, "dim") <- c(n_cells, n_genes)

# overdispersion
eta_real <- rep(NA, n_genes)      

# capture efficiency
catt_real <- rep(NA, n_cells)

# unspliced and spliced counts
Y_u <- rep(NA, n_cells*n_genes)
attr(Y_u, "dim") <- c(n_cells, n_genes)
Y_s <- rep(NA, n_cells*n_genes)
attr(Y_s, "dim") <- c(n_cells, n_genes)

# -----------------------------------------------------
# SIMULATE THE GENE-SPECIFIC DYNAMIC
# -----------------------------------------------------
# --- Capture efficiency
if(typeD == "D1"){
    catt_real <- rep(1, n_cells)
}else{
    catt_real <- runif(n_cells, 0.5, 1)
}
# Solve not identifiability of capture efficiencies
meanCatt_real <- mean(catt_real)
catt_real <- catt_real/meanCatt_real

# --- Overdispersion
if(typeD == "D1"){
    eta_real <- rep(0, n_genes)
}else{
    eta_real <- runif(n_genes,0.5, 1) # eccesso di varianza della distribuzione binomiale negativa dei dati 

}

# --- Rates and steady states parameters
set.seed(seed)
for(g in 1:n_genes){
    alpha_real[g,] <- rep(NA, 2)
    alpha_real[g,1] <- runif(1, 1, 5)
    alpha_real[g,2] <- runif(1, 6, 10)
    # update the parameters according to not identifiability of capture efficiencies
    alpha_real[g,] <- alpha_real[g,]*meanCatt_real

    beta_real[g] <- 1            # fixed beta for identifiability
    gamma_real[g] <- runif(1, 0.5, 1.5)

    # assume beta different from gamma (otherwise we have a different analytic solution of the ODE system --> To DO: implement solutions when beta = gamma)
    if(gamma_real[g] == beta_real[g]){
        gamma_real[g] <- gamma_real[g]-0.1
    }

    u_SSoff_real[g] <- alpha_real[g,1]/beta_real[g]
    s_SSoff_real[g] <- alpha_real[g,1]/gamma_real[g]
    u0_on_real[g] <- u_SSoff_real[g]
    s0_on_real[g] <- s_SSoff_real[g]

    diffSS_off_real[g] <- (alpha_real[g,2]- alpha_real[g,1])/beta_real[g]
}

#--- Switching times and coordinates
set.seed(seed)
tmpT0 <- matrix(NA, n_cells, n_genes)
tmpU0 <- matrix(NA, n_cells, n_genes)
tmpS0 <- matrix(NA, n_cells, n_genes)
for(g in 1:n_genes){
    # take into account the cluster labels
    for(tyT0 in 1:n_typeT0_off){
        cellTy <- which(typeCellT0_off == tyT0)
        t0_off_real[tyT0, g] <- sample(possT0, 1)
        u0_off_real[tyT0, g] <- u0(t0_off = t0_off_real[tyT0, g], t0_on = t0_on_real[g], u0_off = NA, u0_on = NA, k = 0, alpha = alpha_real[g,], beta_real[g])
        s0_off_real[tyT0, g] <- s0(t0_off = t0_off_real[tyT0, g], t0_on = t0_on_real[g], u0_off = u0_off_real[tyT0, g], u0_on = u0_on_real[g], s0_off = NA, s0_on = s0_on_real[g], k = 0, alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        tmpT0[cellTy,g] <- t0_off_real[tyT0, g]
        tmpU0[cellTy,g] <- u0_off_real[tyT0, g]
        tmpS0[cellTy,g] <- s0_off_real[tyT0, g]
    }
}
attr(t0_off_real, "dim") <- c(n_typeT0_off, n_genes)
attr(s0_off_real, "dim") <- c(n_typeT0_off, n_genes)
attr(u0_off_real, "dim") <- c(n_typeT0_off, n_genes)


#--- Elapsed time and associated parameters
set.seed(seed)

# Ensure that cells of the same group have a common underlying structure: generate a group-specific mean
muTypeCell <- matrix(NA, n_cells, n_genes)
for(g in 1:n_genes){
    for(ty in 1:n_typeC){
        cellTy <- which(typeCell == ty)
        theta <- runif(1, 0, 1) # decide in which branch the group-specific mean will be (on, off or lower steady state)
        if(theta < (1-0.01)/2){
            uProp <- u0_on_real[g] + (theta)/( (1-0.01)/2)*(u0_off_real[typeCellT0_off[cellTy][1], g] - u0_on_real[g])
            muTypeCell[cellTy, g] <- log((beta_real[g]*u0_on_real[g]- alpha_real[g, 2])/(beta_real[g]*uProp - alpha_real[g, 2]))

        }else if(theta < (1-0.01)){
            uProp <- u0_off_real[typeCellT0_off[cellTy][1], g] + (theta - ((1-0.01)/2))/((1-0.01)/2)*(u0_on_real[g] - u0_off_real[typeCellT0_off[cellTy][1], g])
            muTypeCell[cellTy, g] <- log((beta_real[g]*u0_off_real[typeCellT0_off[cellTy][1],g]- alpha_real[g, 1])/(beta_real[g]*uProp - alpha_real[g, 1])) + t0_off_real[typeCellT0_off[cellTy][1], g]
        }else{
            muTypeCell[cellTy, g] <- 0
        }
    }
}
# spread the different subgroups around the common group-specific mean
for(g in 1:n_genes){
    for(sty in 1:n_subtypeC){
        subcellTy <- which(subtypeCell == sty)
        t_real[subcellTy, g] <- rnorm(1, mean = muTypeCell[subcellTy[1], g], sd = sdTimes) 
    }
}

# put mass in the lower steady-state
t_real[t_real < 0] <- 0                    

# elapsed time since the last switching and branch of dynamic
tau_real <- t_real
for(g in 1:n_genes){
    for(c in 1:n_cells){
        if(t_real[c,g] <= t0_off_real[typeCellT0_off[c], g]){
            k_real[c,g] <- 2
        }else{
            k_real[c,g] <- 0
            tau_real[c,g] <- t_real[c,g] - t0_off_real[typeCellT0_off[c], g]
        }
    }
}



#--- Spliced and Unspliced counts
for(g in 1:n_genes){
        pos_u_real[,g] <- u(t = t_real[,g], t0_off = tmpT0[,g], t0_on = t0_on_real[g], u0_off = tmpU0[,g], u0_on = u0_on_real[g], k = k_real[,g], alpha = alpha_real[g,], beta = beta_real[g])
        pos_s_real[,g] <- s(t = t_real[,g], t0_off = tmpT0[,g], t0_on = t0_on_real[g], u0_off = tmpU0[,g], u0_on = u0_on_real[g], s0_off = tmpS0[,g], s0_on = s0_on_real[g], k = k_real[,g], alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        if(typeD == "D1"){ 
            # Poisson
            Y_u[,g] <- rpois(n_cells, lambda = pos_u_real[,g]*catt_real)
            Y_s[,g] <- rpois(n_cells, lambda = pos_s_real[,g]*catt_real)
        }else{
            # Negative Binomial
            Y_u[,g] <- rnbinom(n_cells, size= rep(1/eta_real[g], n_cells), mu = pos_u_real[,g]*catt_real) 
            Y_s[,g] <- rnbinom(n_cells, size= rep(1/eta_real[g], n_cells), mu = pos_s_real[,g]*catt_real) 
        }   
}

# -----------------------------------------------------
# SAVE RESULTS
# -----------------------------------------------------
save.image(file = paste0(path, "/", nameSim, "_", n_genes, ".RData"))
