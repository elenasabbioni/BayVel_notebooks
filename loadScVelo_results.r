# -----------------------------------------------------------
# Script to load scVelo results (following postprocessing of scVelo)
# And obtain the quantities necessary for the plots
# -----------------------------------------------------------

# --- load results
if(nameSim == "moments"){
    # real data
    path <- paste0(pathToResults, "/", typeSIM, "/", nameSim, "/output/res_", typeSIM)
    scVeloName <- paste0(pathToResults, "/", typeSIM, "/", nameSim, "/output/res_", typeSIM)
}else{
    # simulated data
    path <- paste0(pathToResults, "/", nameSim, "/", nameSim, "_", n_genes, "_", modelScVel)
    scVeloName <- paste0(pathToResults, "/", nameSim , "/output/res_", typeSIM, "_", nameSim, "_", n_genes, "genes_", modelScVel, "_scVelo")
}

# rates and steady states
scv$fit_alpha <- as.matrix(read.csv(paste(scVeloName, "_fit_alpha.csv", sep = ""), header = TRUE))
n_genes <- dim(scv$fit_alpha)[2]
colnames(scv$fit_alpha) <- as.character(seq(1,n_genes))
scv$fit_beta <- as.matrix(read.csv(paste(scVeloName, "_fit_beta.csv", sep = ""), header = TRUE))
colnames(scv$fit_beta) <- as.character(seq(1,n_genes))
scv$fit_gamma <- as.matrix(read.csv(paste(scVeloName, "_fit_gamma.csv", sep = ""), header = TRUE))
colnames(scv$fit_gamma) <- as.character(seq(1,n_genes))
scv$fit_steady_s <- as.matrix(read.csv(paste(scVeloName, "_fit_steady_s.csv", sep = ""), header = TRUE))
colnames(scv$fit_steady_s) <- as.character(seq(1,n_genes))
scv$fit_steady_u <- as.matrix(read.csv(paste(scVeloName, "_fit_steady_u.csv", sep = ""), header = TRUE))
colnames(scv$fit_steady_u) <- as.character(seq(1,n_genes))

# time and switching time
scv$fit_t <- as.matrix(read.csv(paste(scVeloName, "_fit_t.csv", sep = ""), header = TRUE))
colnames(scv$fit_t) <- as.character(seq(1,n_genes))
scv$fit_t_ <- as.matrix(read.csv(paste(scVeloName, "_fit_t_.csv", sep = ""), header = TRUE))
colnames(scv$fit_t_) <- as.character(seq(1,n_genes))
scv$fit_tau <- as.matrix(read.csv(paste(scVeloName, "_fit_tau.csv", sep = ""), header = TRUE))
colnames(scv$fit_tau) <- as.character(seq(1,n_genes))
scv$fit_tau_ <- as.matrix(read.csv(paste(scVeloName, "_fit_tau_.csv", sep = ""), header = TRUE))
colnames(scv$fit_tau_) <- as.character(seq(1,n_genes))
scv$fit_scaling <- as.matrix(read.csv(paste(scVeloName, "_fit_scaling.csv", sep = ""), header = TRUE))
colnames(scv$fit_scaling) <- as.character(seq(1,n_genes))
scv$fit_s0 <- as.matrix(read.csv(paste(scVeloName, "_fit_s0.csv", sep = ""), header = TRUE))
colnames(scv$fit_s0) <- as.character(seq(1,n_genes))
scv$fit_u0 <- as.matrix(read.csv(paste(scVeloName, "_fit_u0.csv", sep = ""), header = TRUE))
colnames(scv$fit_u0) <- as.character(seq(1,n_genes))
scv$top_genes <- as.matrix(read.csv(paste(scVeloName, "_top_genes.csv", sep = ""), header = TRUE))
# velocity
scv$velocity <- as.matrix(read.csv(paste(scVeloName, "_velocity.csv", sep = ""), header = TRUE))
colnames(scv$velocity) <- as.character(seq(1,n_genes))
scv$velocity_u <- as.matrix(read.csv(paste(scVeloName, "_velocity_u.csv", sep = ""), header = TRUE))
colnames(scv$velocity_u) <- as.character(seq(1,n_genes))
# pre-processed counts
scv$Mu <- as.matrix(read.csv(paste(path, "_Mu.csv", sep = ""), header = TRUE))
colnames(scv$Mu) <- as.character(seq(1,n_genes))
scv$Ms <- as.matrix(read.csv(paste(path, "_Ms.csv", sep = ""), header = TRUE))
colnames(scv$Ms) <- as.character(seq(1,n_genes))
scv$var <- as.matrix(read.csv(paste(scVeloName, "_/var.csv", sep = ""), header = TRUE))

# groups and subgroups labels
scv$subtypeCellReal <- as.matrix(read.csv(paste(scVeloName, "_/obs.csv", sep = ""), header = TRUE))[, "clusters"]
scv$subtypeCell <- as.numeric(as.factor(scv$subtypeCellReal))
scv$typeCellReal <- as.matrix(read.csv(paste(scVeloName, "_/obs.csv", sep = ""), header = TRUE))[, "clusters"]
scv$typeCell <- as.numeric(as.factor(scv$typeCellReal))

# number of cells
n_cells <- dim(scv$Mu)[1]

# switching clusters
if(grepl("SW1", nameSim)){
    scv$typeCellT0_off <- rep(1, n_cells)
}else{
    scv$typeCellT0_off <- scv$typeCell
}




# --- Process the parameters according to the steps done in scVelo plotting functions and obtain the parameters used in BayVel plotting functions
scv$beta <- t(scv$fit_beta * scv$fit_scaling ) # scale beta with scaling parameters
scv$alpha <- t(rbind(0, scv$fit_alpha)) # alpha_off assumed to be 0

# compute switching points accordingly to the processed rates
scv$u0_off <- u0_MCMC(t0_off = array(scv$fit_t_, dim = c(1, n_genes, 1)), t0_on = 0,k = 0, alpha = array(scv$alpha, dim = c(n_genes, 2, 1)), beta = array(scv$beta, dim = c(n_genes, 1)))
scv$s0_off <- s0_MCMC(t0_off = array(scv$fit_t_, dim = c(1, n_genes, 1)), t0_on = 0, k = 0, alpha = array(scv$alpha, dim = c(n_genes, 2, 1)), beta = array(scv$beta, dim = c(n_genes, 1)), gamma = array(scv$fit_gamma, dim = c(n_genes, 1)))

scv$u0_on <- array(0, dim = c(1, n_genes, 1))
scv$s0_on <- array(0, dim = c(1, n_genes, 1))

# define k and uniform the time notation of scVelo with the one of BayVel
scv$k <- matrix(0, nrow = n_cells, ncol = n_genes)
scv$fit_t0_matrix <- matrix(rep(scv$fit_t_[1,], n_cells), nrow = n_cells, byrow = TRUE)
scv$k[which(scv$fit_t < scv$fit_t0_matrix, arr.ind = TRUE)] <- 2
scv$tau <- scv$fit_t
scv$tau[which(scv$k == 0, arr.ind = TRUE)] <- scv$fit_t[which(scv$k == 0, arr.ind = TRUE)] - scv$fit_t0_matrix[which(scv$k == 0, arr.ind = TRUE)]
scv$k <- array(scv$k, dim = c(n_cells, n_genes, 1))
scv$tau <- array(scv$tau, dim = c(n_cells, n_genes, 1))

# compute the mean of Mu and Ms
scv$res <- u_and_s_withTauMCMC(tau = scv$tau, u0_off = scv$u0_off, u0_on = NA, s0_off = scv$s0_off, s0_on = NA, k = scv$k, alpha = array(scv$alpha, dim = c(n_genes, 2, 1)), beta = array(scv$beta, dim = c(n_genes, 1)), gamma = matrix(scv$fit_gamma, nrow = n_genes), subtypeCell = seq(1, n_cells), typeCellT0_off = rep(1, n_cells))
scv$pos_u <- scv$res$u
scv$pos_s <- scv$res$s

# scale the switching time accordingly to the plotting functions of scVelo
scv$u0_off <- scv$u0_off* array(scv$fit_scaling, dim = c(1, n_genes, 1))  +  array(scv$fit_u0, dim = c(1, n_genes, 1)) 
scv$s0_off <- scv$s0_off + array(scv$fit_s0, dim = c(1, n_genes, 1)) 
# scale the mean position of Mu and Ms accordingly to the plotting functions of scVelo
scv$pos_u[,,1] <- scv$pos_u[,,1] * matrix(scv$fit_scaling[1,], nrow = n_cells, ncol = n_genes, byrow = TRUE) + matrix(scv$fit_u0[1,], nrow = n_cells, ncol = n_genes, byrow = TRUE)
scv$pos_s[,,1] <- scv$pos_s[,,1] + matrix(scv$fit_s0[1,], nrow = n_cells, ncol = n_genes, byrow = TRUE)    
attr(scv$pos_u, "dim") <- c(n_cells, n_genes)
attr(scv$pos_s, "dim") <- c(n_cells, n_genes)        

           
        

                