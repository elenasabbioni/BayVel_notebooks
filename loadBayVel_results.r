# -----------------------------------------------------------
# Script to load BayVel results (posterior chains for the parameters)
# And obtain the quantities necessary for the plots#
# -----------------------------------------------------------

# Set if we will consider all the posterior sample of the MCMC, or just the sample that has the highest posterior likelihood
maxLog <- tryCatch({
    maxLog
    }, error = function(e) {maxLog = FALSE})


# load results
if(typeSIM != "sim"){
    # real data
    BayVelName <- paste0(pathToResults, "/", typeSIM, "/filter/")
    bv$typeCellReal <- loadObs(paste0(pathToYourDirectory, "/real data/",  typeSIM, "/filter"))[, "clusters"]
    bv$nameG <- nameGenes(paste0(pathToYourDirectory, "/real data/",  typeSIM, "/filter"))
}else{
    # simulated data
    BayVelName <- paste0(pathToResults, "/")
}

load(paste0(BayVelName, nameSim, "/output/res_", typeSIM, "_", nameSim, "_", n_genes, "genes_", mcmc, ".RData"), envir = bv)

# load total likelihood for the different iteration of mcmc
logLikFile <- paste0(BayVelName, nameSim, "/output/loglik_", typeSIM, "_", nameSim, "_", n_genes, "genes_", mcmc, ".csv")
if(file.exists(logLikFile)){
    bv$loglikTOT <- read.csv(logLikFile)
}

# set all the parameters that will be nedded
bv$n_typeC <- length(unique(bv$typeCell))
bv$n_subtypeC <- length(unique(bv$subtypeCell))
bv$n_typeT0_off <- dim(bv$T0_off_chain)[1]
bv$n_cells <- dim(bv$Y_u)[1]

bv$SampleToSave <- trunc((bv$mcmcIter - bv$mcmcBurnin)/bv$mcmcThin)


if(maxLog == TRUE){
    bv$iter <- which.max(bv$loglikTOT[,2])  # take the index of the iteration with the highest likelihood
}else{
    bv$iter <- seq(1, bv$SampleToSave)      # consider all the iterations
}

print("Everything loaded")

source(paste0(pathToYourDirectory, "/functions.R"))

# --- Resolve not-identifiability of capture efficiency
bv$Catt_chain <- invlogit(bv$LogitCatt_chain)
bv$meanCatt_Iter <- apply(bv$Catt_chain, FUN = mean, MARGIN = 2)
bv$Catt_chain <- sweep(bv$Catt_chain, 2, bv$meanCatt_Iter, "/")
bv$LogSS_chain[1:3, , ] <- sweep(bv$LogSS_chain[1:3, , ], 3, log(bv$meanCatt_Iter), "+")

# rescale appropriately all the real quantities (we should have already done this produce in the dataSimulation file, but we check, in case we missed it. If the procedure is already done, we are just rescling everything by 1) 
if(typeSIM == "sim"){
    bv$pos_u_real <- bv$pos_u_real*mean(bv$catt_real)
    bv$pos_s_real <- bv$pos_s_real*mean(bv$catt_real)
    bv$u0_off_real <- bv$u0_off_real*mean(bv$catt_real)
    bv$s0_off_real <- bv$s0_off_real*mean(bv$catt_real)
    bv$alpha_real <- bv$alpha_real*mean(bv$catt_real)   
    bv$catt_real <- bv$catt_real/mean(bv$catt_real)
    bv$s_SSoff_real <- bv$alpha_real[, 1]/bv$gamma_real
    bv$u_SSoff_real <- bv$alpha_real[, 1]/bv$beta_real
    bv$u_SSon_real <- bv$alpha_real[, 2]/bv$beta_real
    bv$s_SSon_real <- bv$alpha_real[, 2]/bv$gamma_real
}

bv$alpha_chain <- rep(NA, bv$n_genes*2*dim(bv$LogSS_chain)[3])
attr(bv$alpha_chain, "dim") <-  c(bv$n_genes, 2, dim(bv$LogSS_chain)[3])
bv$beta_chain <- rep(NA, bv$n_genes*dim(bv$LogSS_chain)[3])
attr(bv$beta_chain, "dim") <-  c(bv$n_genes, dim(bv$LogSS_chain)[3])
bv$gamma_chain <- rep(NA, bv$n_genes*dim(bv$LogSS_chain)[3])
attr(bv$gamma_chain, "dim") <-  c(bv$n_genes, dim(bv$LogSS_chain)[3])

bv$alpha_chain[,1,] <- exp(bv$LogSS_chain[1,,] + bv$LogSS_chain[4, ,])
bv$alpha_chain[,2,] <- bv$alpha_chain[,1,] + exp(bv$LogSS_chain[3,,] + bv$LogSS_chain[4, ,])
bv$gamma_chain <- bv$alpha_chain[,1,]/exp(bv$LogSS_chain[2,,])
bv$beta_chain <- exp(bv$LogSS_chain[4, ,])


# --- Extract just the iterations we are interested in 
# rates and steady state coordinates
bv$alpha_chain <- bv$alpha_chain[,, bv$iter, drop = FALSE]
bv$beta_chain <- bv$beta_chain[,bv$iter, drop = FALSE]
bv$gamma_chain <- bv$gamma_chain[,bv$iter, drop = FALSE ]

bv$uSS_off_chain = bv$alpha_chain[, 1, ]
bv$sSS_off_chain = bv$alpha_chain[, 1, ]/bv$gamma_chain
bv$uSS_on_chain = bv$alpha_chain[, 2, ]
bv$sSS_on_chain = bv$alpha_chain[, 2, ]/bv$gamma_chain
dim(bv$uSS_off_chain) <- c(bv$n_genes,  length(bv$iter))
dim(bv$sSS_off_chain) <- c(bv$n_genes,  length(bv$iter))
dim(bv$uSS_on_chain) <- c(bv$n_genes,  length(bv$iter))
dim(bv$sSS_on_chain) <- c(bv$n_genes,  length(bv$iter))


# switching time and associated coordinates
bv$T0_off_chain <- bv$T0_off_chain[, ,bv$iter, drop = FALSE]
bv$u0_off_chain <- u0_MCMC(t0_off = bv$T0_off_chain, t0_on = 0, k = 0, alpha = bv$alpha_chain, beta = bv$beta_chain)
bv$s0_off_chain <- s0_MCMC(t0_off = bv$T0_off_chain, t0_on = 0, k = 0, alpha = bv$alpha_chain, beta = bv$beta_chain, gamma = bv$gamma_chain)

# subgroup time and associated coordinates
if(max(bv$k_chain) == 1){
    bv$k_chain <- bv$k_chain[, ,bv$iter, drop = FALSE] *2
}

bv$Tau_chain <- bv$Tau_chain[, ,bv$iter, drop = FALSE]
bv$TStar_withM_chain <- bv$TStar_withM_chain[, , bv$iter, drop = FALSE]

bv$res <- u_and_s_withTauMCMC(tau = bv$Tau_chain, u0_off = bv$u0_off_chain, u0_on = NA, s0_off = bv$s0_off_chain, s0_on = NA, k = bv$k_chain, alpha = bv$alpha_chain, beta = bv$beta_chain, gamma = bv$gamma_chain, subtypeCell = bv$subtypeCell, typeCellT0_off = bv$typeCellT0_off)
bv$u_chain <- bv$res$u
bv$s_chain <- bv$res$s
dim(bv$u_chain) <- c(bv$n_subtypeC, bv$n_genes, length(bv$iter))
dim(bv$s_chain) <- c(bv$n_subtypeC, bv$n_genes, length(bv$iter))

# velocity 
bv$v <- array(NA, dim = c(bv$n_subtypeC, bv$n_genes, length(bv$iter)))
for(sty in 1:bv$n_subtypeC){
    bv$v[sty, , ] = bv$beta_chain*bv$u_chain[sty, ,] - bv$gamma_chain*bv$s_chain[sty, , ]
}

# overdispersion
bv$LogEta_chain <- bv$LogEta_chain[,bv$iter, drop = FALSE]
bv$Eta_chain <- exp(bv$LogEta_chain)

# capture efficiency
bv$Catt_chain <- bv$Catt_chain[, bv$iter, drop = FALSE]

# remove all the not necessary elements
bv$elEnv <- ls(envir = bv)
bv$elToRemove <- setdiff(bv$elEnv,  c("n_genes", "n_cells", "n_typeC", "n_subtypeC", "n_typeT0_off", "subtypeCell", "typeCell", "typeCellT0_off", "SampleToSave", "iter", "loglikTOT", "mcmcBurnin", "mcmcIter", "mcmcThin", "Catt_chain", "meanCatt_Iter", "LogSS_chain", "alpha_chain", "beta_chain", "gamma_chain", "uSS_off_chain", "sSS_off_chain", "uSS_on_chain", "sSS_on_chain", "T0_off_chain", "LogT0_off_chain", "u0_off_chain", "s0_off_chain", "k_chain", "Tau_chain",  "TStar_withM_chain", "u_chain", "s_chain", "v", "Eta_chain", "LogEta_chain", "Y_u", "Y_s", "typeCellReal", "nameG"))

for(el in bv$elToRemove){
    rm(list = el, envir = bv)
}

