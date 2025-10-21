# ------------- ANALYTIC FUNCTIONS ----------------

# -------------------------------------------------------------------------
# Function: u
# Purpose:
#   Computes the amount of *unspliced RNA* for multiple cells for a fixed gene,
#   based on RNA velocity ODE's model
#
# Parameters:
#   t        : Current time (scalar or vector).
#   t0_off   : Initial time when the gene switches to the "off" state.
#   t0_on    : Initial time when the gene switches to the "on" state (default = 0).
#   u0_off   : Initial unspliced RNA amount in the "off" state.
#   u0_on    : Initial unspliced RNA amount in the "on" state (optional).
#               If NA, it is set to the OFF steady-state value alpha_off/beta.
#   k        : scalar or vector of same length of t indicating the transcriptional state
#               for each cell (0 = off, 2 = on).
#   alpha    : Two-element vector of transcription rates:
#                 alpha[1] = rate when gene is off,
#                 alpha[2] = rate when gene is on.
#   beta     : Splicing rate constant (scalar). Equal to 1 by default.
#
# Returns:
#   A numeric vector of the same length of t containing the unspliced RNA values
#   for each cell at time t.
# -------------------------------------------------------------------------

u <- function(t, t0_off, t0_on = 0, u0_off, u0_on = NA, k, alpha, beta = 1) {
  
  # If not provided, initialize u0_on as the steady-state level for the "off" state
  if (any(is.na(u0_on))) {
    u0_on <- alpha[1] / beta
  }
  
  n_cells <- length(k)  # number of cells
  
  # Compute elapsed time (tau) depending on cell state
  tau <- ifelse(k == 2, t - t0_on, t - t0_off)
  
  # Precompute exponential decay terms for both states
  expBetaOFF <- exp(-beta * tau[k == 0])
  expBetaON  <- exp(-beta * tau[k == 2])
  
  # Initialize result vector
  res <- rep(NA, n_cells)
  
  # Assign appropriate initial unspliced RNA based on cell state
  u0 <- ifelse(k == 2, u0_on, u0_off)
  
  # Compute unspliced RNA 
  res[k == 0] <- u0[k == 0] * expBetaOFF + alpha[1]/beta * (1 - expBetaOFF)
  res[k == 2] <- u0[k == 2] * expBetaON  + alpha[2]/beta * (1 - expBetaON)
  
  return(res)
}

# -------------------------------------------------------------------------
# Function: s
# Purpose:
#   Computes the amount of *spliced RNA* for multiple cells for a fixed gene,
#   based on RNA velocity ODE's model
#
# Parameters:
#   t        : Current time (scalar or vector).
#   t0_off   : Initial time when the gene switches to the "off" state.
#   t0_on    : Initial time when the gene switches to the "on" state (default = 0).
#   u0_off   : Initial unspliced RNA amount in the "off" state.
#   u0_on    : Initial unspliced RNA amount in the "on" state (optional).
#               If NA, it is set to the OFF steady-state value alpha_off/beta.
#   s0_off   : Initial spliced RNA amount in the "off" state.
#   s0_on    : Initial spliced RNA amount in the "on" state (optional).
#               If NA, it is set to the OFF steady-state value alpha_off/gamma.
#   k        : scalar or vector of same length of t indicating the transcriptional state
#               for each cell (0 = off, 2 = on).
#   alpha    : Two-element vector of transcription rates:
#                 alpha[1] = rate when gene is off,
#                 alpha[2] = rate when gene is on.
#   beta     : Splicing rate constant (scalar). Equal to 1 by default.
#   gamma    : Degradation rate constant (scalar). 
#
# Returns:
#   A numeric vector of the same length of t containing the spliced RNA values
#   for each cell at time t.
# -------------------------------------------------------------------------
s <- function(t, t0_off, t0_on = 0, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta = 1, gamma){
  
  # If not provided, initialize u0_on and s0_on  as the steady-state level for the "off" state
  if(is.na(u0_on)){
    u0_on = alpha[1]/beta
  }
  if(is.na(s0_on)){
    s0_on = alpha[1]/gamma
  }

  n_cells <- length(k) # number of cells
  
  # Compute elapsed time (tau) depending on cell state
  tau <- ifelse(k == 2, t - t0_on, t - t0_off)
  # Precompute exponential decay terms for both states
  exp_gammaOFF <- exp(-gamma*tau[k == 0])
  exp_gammaON <- exp(-gamma*tau[k == 2])
  exp_betaOFF <- exp(-beta*tau[k == 0])
  exp_betaON <- exp(-beta*tau[k == 2])
  # Initialize result vector
  res <- rep(NA, n_cells)
  # Assign appropriate initial unspliced and spliced RNA based on cell state
  u0 <- ifelse(k == 2, u0_on,  u0_off)
  s0 <- ifelse(k == 2, s0_on,  s0_off)
  
  # Compute unspliced RNA 
  res[k == 0] <- s0[k == 0]*exp_gammaOFF + alpha[1]/gamma*(1-exp_gammaOFF)+ (alpha[1]-beta*u0[k == 0])/(gamma-beta)*(exp_gammaOFF - exp_betaOFF)
  res[k == 2] <- s0[k == 2]*exp_gammaON  + alpha[2]/gamma*(1-exp_gammaON) + (alpha[2]-beta*u0[k == 2])/(gamma-beta)*(exp_gammaON  - exp_betaON)
  return(res)
}

# -------------------------------------------------------------------------
# Function: s_di_u
# Purpose:
#   Computes the amount of *spliced RNA* for multiple cells for a fixed gene,
#   given the amount of the *unspliced RNA*. 
#   Note that this formula holds when beta is different from gamma.
#
# Parameters:
#   u        : Current value of unspliced RNA for he different cells (scalar or vector).
#   u0_off   : Initial unspliced RNA amount in the "off" state.
#   u0_on    : Initial unspliced RNA amount in the "on" state (optional).
#               If NA, it is set to the OFF steady-state value alpha_off/beta.
#   k        : scalar or vector of same length of t indicating the transcriptional state
#               for each cell (0 = off, 2 = on).
#   alpha    : Two-element vector of transcription rates:
#                 alpha[1] = rate when gene is off,
#                 alpha[2] = rate when gene is on.
#   beta     : Splicing rate constant (scalar). Equal to 1 by default.
#   gamma    : Degradation rate constant (scalar). 
#
# Returns:
#   A numeric vector of the same length of t containing the spliced RNA values
#   for each cell with unspliced RNA u.
# -------------------------------------------------------------------------
s_di_u <- function(u, u0_off, s0_off, u0_on = NA, s0_on = NA, k, alpha, beta = 1, gamma){
  # check that beta not equal to gamma
  if(beta == gamma){
    warning("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }

  # If not provided, initialize u0_on and s0_on  as the steady-state level for the "off" state
  if(is.na(u0_on)){
    u0_on = alpha[1]/beta
  }
  if(is.na(s0_on)){
    s0_on = alpha[1]/gamma
  }

  # set the parameters depending of the state
  if(k == 0){
    a <- alpha[1]
    u0 <- u0_off
    s0 <- s0_off
  }else{
    a <- alpha[2]
    u0 <- u0_on
    s0 <- s0_on
  }

  # compute the spliced RNA
  res <- (s0 - a/gamma + (a-beta*u0)/(gamma-beta))*((beta*u-a)/(beta*u0 - a))^(gamma/beta) -  (a*beta)/(gamma*(gamma-beta)) + beta/(gamma-beta)*u
  
  return(res)
}

# -------------------------------------------------------------------------
# Function: u0
# Purpose:
#   Computes the amount of u-coordinate of the initial point (for on or off phase) for a fixed gene.
#
# Parameters:
#   t0_off   : Initial time when the gene switches to the "off" state (optional).
#               If NA, it is set to Infinity, meaning that we reach the upper steady state and do not enter in the repressive phase.
#   t0_on    : Initial time when the gene switches to the "on" state (default = 0).
#   u0_off   : Initial unspliced RNA amount in the "off" state (optional).
#               If NA and k == 2, it is set to the ON steady-state value alpha_on/beta.
#   u0_on    : Initial unspliced RNA amount in the "on" state (optional).
#               If NA and k == 0, it is set to the OFF steady-state value alpha_off/beta.
#   k        : scalar indicating if we want to compute the initial point for the on (k = 2) or
#               off (k = 0) phase.
#   alpha    : Two-element vector of transcription rates:
#                 alpha[1] = rate when gene is off,
#                 alpha[2] = rate when gene is on.
#   beta     : Splicing rate constant (scalar). Equal to 1 by default.
#
# Returns:
#   A scalar with the computed u0.
# -------------------------------------------------------------------------
u0 <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, k, alpha, beta){

  # we want to compute the initial u-coordinate for the repressive phase
  if(k == 0){
    if(any(is.na(u0_on))){
      u0_on <- alpha[1]/beta      # the initial u-coordinate for the inductive phase coincides with the lower steady state 
    }

    if(any(is.na(t0_off))){
      t <- rep(NA, length(k))
      t[which(is.na(t0_off))] <- rep(Inf, length(which(is.na(t0_off)))) # the cells do not switch and remain in the inductive state, such that the upper steady state is reached
    }else{
      t <- t0_off # the upper steady state is not reached and the dynamic switch to the repressive phase
    }
    res <- u(t = t, t0_off = t0_off, t0_on = t0_on, u0_off = NA, u0_on = u0_on, k = rep(2, length(t)), alpha = alpha, beta = beta)
  }

  # we want to compute the initial u-coordinate for the inductive phase
  if(k == 2){
    if(any(is.na(u0_off))){
      u0_off <- alpha[2]/beta     # the initial u-coordinate for the inductive phase coincides with the lower steady state 
    }
    res <- u(t = Inf, t0_off = 0, t0_on = NA, u0_off = u0_off, u0_on = NA, k = rep(0, length(t)), alpha = alpha, beta = beta)
  }
  return(res)
}

# -------------------------------------------------------------------------
# Function: s0
# Purpose:
#   Computes the amount of s-coordinate of the initial point (for on or off phase) for a fixed gene.
#
# Parameters:
#   t0_off   : Initial time when the gene switches to the "off" state (optional).
#               If NA, it is set to Infinity, meaning that we reach the upper steady state and do not enter in the repressive phase.
#   t0_on    : Initial time when the gene switches to the "on" state (default = 0).
#   u0_off   : Initial unspliced RNA amount in the "off" state (optional).
#               If NA and k == 2, it is set to the ON steady-state value alpha_on/beta.
#               If NA and k == 0, the initial point is computed with the function u0()
#   u0_on    : Initial unspliced RNA amount in the "on" state (optional).
#               If NA and k == 0, it is set to the OFF steady-state value alpha_off/beta.
#   s0_off   : Initial spliced RNA amount in the "off" state (optional).
#               If NA and k == 2, it is set to the ON steady-state value alpha_on/gamma
#   s0_on    : Initial spliced RNA amount in the "on" state (optional).
#               If NA and k == 0, it is set to the OFF steady-state value alpha_off/gamma
#   k        : scalar indicating if we want to compute the initial point for the on (k = 2) or
#               off (k = 0) phase.
#   alpha    : Two-element vector of transcription rates:
#                 alpha[1] = rate when gene is off,
#                 alpha[2] = rate when gene is on.
#   beta     : Splicing rate constant (scalar). Equal to 1 by default.
#   gamma    : Degradation rate constant (scalar). 

# Returns:
#   A scalar with the computed s0.
# -------------------------------------------------------------------------
s0 <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k, alpha, beta, gamma){
  if(k == 0){
    if(any(is.na(s0_on))){
      s0_on <- alpha[1]/gamma
    }
    if(any(is.na(u0_on))){
      u0_on <- alpha[1]/beta
    }

    if(any(is.na(t0_off))){
      t <- rep(NA, length(k))
      t[which(is.na(t0_off))] <. rep(Inf, length(which(is.na(t0_off)))) # the cells do not switch and remain in the inductive state, such that the upper steady state is reached
    }else{
      t <- t0_off # the upper steady state is not reached and the dynamic switch to the repressive phase
    }
    if(any(is.na(u0_off))){
      # compute the initial u-coordinate
      u0_off <- u0(t0_off = t0_off, t0_on = t0_on, u0_off = NA, u0_on = u0_on, k = 0, alpha = alpha, beta = beta)
    }
    
    # compute s0
    res <- s(t = t, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(t)), alpha, beta, gamma)
  }

  if(k == 2){
    if(any(is.na(s0_off))){
      s0_off <- alpha[2]/gamma
    }
    if(any(is.na(u0_off))){
      u0_off <- alpha[2]/beta
    }

    res <- s(rep(Inf, length(t0_on)), t0_off = 0, t0_on = NA, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t0_on)), alpha = alpha, beta = beta, gamma = gamma)
    
  }
  
  return(res)
}

# -------------------------------------------------------------------------
# Function: u_and_s_withTauMCMC
# Purpose:
#   Compute unspliced (u) and spliced (s) RNA levels for multiple genes,
#   cells and MCMC samples (or other 3D structure) given per-cell elapsed
#   times (tau) and kinetic parameters. This implements the closed-form
#   solutions of linear ODEs for the 2-state transcription model
#   (states: off (k==0) and on (k==2)).
#
# Inputs (expected shapes):
#   tau             : array (dim: [subtype, n_genes, n_samples] or similar)
#                     elapsed time for each subgroup, gene and MCMC sample.
#   u0_off          : array with initial unspliced RNA in the "off" state
#   u0_on           : array with initial unspliced RNA in the "on" state.
#   s0_off          : array with initial spliced RNA in the "off" state.
#   s0_on           : array with initial spliced RNA in the "o" state.
#   k               : integer array with same dim as tau indicating
#                     transcriptional state per element (0 = off, 2 = on).
#   alpha           : numeric array of transcription rates. Expected to
#                     have dimensions [genes, state_index(1/2), samples]
#                     (uses alpha[,1,] for off and alpha[,2,] for on).
#   beta            : array for splicing rate ([n_genes, n_samples]).
#   gamma           : array for degradation rate ([n_genes, n_samples]).
#   subtypeCell     : vector indicating subtype index for each cell
#                     (used to assign off-state initial values by subtype).
#   typeCellT0_off  : vector indicating switching time labels 
#
# Returns:
#   A list with:
#     $u : array of same shape as tau with computed unspliced RNA values.
#     $s : array of same shape as tau with computed spliced   RNA values.
# -------------------------------------------------------------------------
u_and_s_withTauMCMC <- function(tau, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma, subtypeCell, typeCellT0_off){    
  
  # Set inductive initial values equal to off steady state coordinates, if not assigned.
  if(any(is.na(u0_on))){
    u0_on <- alpha[, 1, ]/beta
  }
  if(any(is.na(s0_on))){
    s0_on <- alpha[, 1, ]/gamma
  }
  # expand the elements, such that they will have the same dimension of tau
  u0_on = sweep(k == 2, MARGIN = c(2, 3), STATS = u0_on, FUN = "*")
  s0_on = sweep(k == 2, MARGIN = c(2, 3), STATS = s0_on, FUN = "*")

  # expand also the off initial coordinates, such that they will have the same dimension of tau
  u0_off_new <- array(NA, dim = dim(tau))
  s0_off_new <- array(NA, dim = dim(tau))
  for(sty in unique(subtypeCell)){
    # assign the corresponding off-state initial arrays to all elements
    # in each subtypes
      tyT0_off <- typeCellT0_off[which(subtypeCell == sty)[1]]
      u0_off_new[sty, , ] = u0_off[tyT0_off, , ]
      s0_off_new[sty, , ] = s0_off[tyT0_off, , ]
  }
  # Compose the full initial conditions: choose on vs off initial values
  # depending on k (state indicator).
  u0 = u0_on * (k == 2) + u0_off_new * (k == 0)
  s0 = s0_on * (k == 2) + s0_off_new * (k == 0)
  
  # upper-steady state u-coordinate
  uSS_on <- alpha[,2,]/beta
  
  # Exponential decay factors for splicing (beta) and degradation (gamma) rates.
  expBeta <- exp(-sweep(tau, MARGIN = c(2, 3), STATS = beta, FUN = "*"))    
  expGamma <- exp(-sweep(tau, MARGIN = c(2, 3), STATS = gamma, FUN = "*"))
  p2 <- sweep(k == 0, MARGIN = c(2, 3), STATS = alpha[, 1, ]/beta, FUN = "*")*(1-expBeta)

  # Transcription rate depending on the state
  alpha <- sweep(k == 2, MARGIN = c(2, 3), STATS = alpha[, 2, ], FUN = "*") + sweep(k == 0, MARGIN = c(2, 3), STATS = alpha[, 1, ], FUN = "*")

  # unspliced and spliced values
  resS <- s0*expGamma + sweep(alpha, MARGIN = c(2, 3), STATS = gamma, FUN = "/")*(1-expGamma) + sweep((alpha-sweep(u0, MARGIN = c(2, 3), STATS = beta, FUN = "*")), MARGIN = c(2, 3), STATS = (gamma-beta), FUN = "/")*(expGamma - expBeta)
  resU <- u0*expBeta + p2 + sweep(k==2, MARGIN = c(2, 3), STATS = uSS_on, FUN = "*")*(1-expBeta)

  return(list(u = resU, s = resS))
}

# -------------------------------------------------------------------------
# Function: u0_MCMC
# Purpose:
#   Compute initial unspliced RNA (u0) for MCMC-shaped arrays.
#
# Arguments:
#   t0_off : 3D array [n_switching_clusters, n_genes, n_samples] with time of OFF-switch.
#   t0_on  : scalar or array; default 0, for the starting time of the on phase.
#   k      : scalar or vector of states (0 or 2). 
#   alpha  : array [n_genes, 2, n_samples] for transcription rate (alpha[,1,] = off, alpha[,2,] = on).
#   beta   : scalar or array [n_genes, n_samples] for splicing rate. By default equal to 1.
#
# Returns:
#   Array of same shape as t0_off with u0 values for each element.
#   Note that we can compute just one between u0_off and u0_on (for different genes) for each call of this function. 
# -------------------------------------------------------------------------
u0_MCMC <- function(t0_off = NA, t0_on = 0, k, alpha, beta = 1){
  # If multiple k values provided, keep only unique values.
  if(length(k) > 1){
      k <- unique(k)
  }

  # initial point for repressive phase
  if(k == 0){
    # transform u0_on according to the dimension of t0_off
    u0_on <- alpha[, 1, ]/beta
    u0_on <- array(rep(as.vector(u0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    u0_on <- aperm(u0_on, c(3, 1, 2))
    
    t <- t0_off # we assume the switch from on to off phase occurs
    t0_on <- array(t0_on, dim = dim(u0_on))
    tau <- t - t0_on    # elapsed time in the dynamic    

    # coordinate of the on steady state
    uSS_ON <- alpha[, 2,]/beta
    uSS_ON <- array(rep(as.vector(uSS_ON), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    uSS_ON <- aperm(uSS_ON, c(3, 1, 2))

    # modify the splicing rate accordingly to the dimension of t0:off
    beta <- array(rep(as.vector(beta), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    beta <- aperm(beta, c(3, 1, 2))

    expBetaOFF <- exp(-beta*tau)
    # u0_off
    res <- u0_on*expBetaOFF + uSS_ON*(1-expBetaOFF)
  }

  # initial point for inductive phase
  if(any(k == 2)){
    res <- alpha[,1,]/beta
  }
  return(res)
}

# -------------------------------------------------------------------------
# Function: s0_MCMC
# Purpose:
#   Compute initial spliced RNA (s0) for MCMC-shaped arrays.
#
# Arguments:
#   t0_off : 3D array [n_switching_clusters, n_genes, n_samples] with time of OFF-switch.
#   t0_on  : scalar or array; default 0, for the starting time of the on phase.
#   k      : scalar or vector of states (0 or 2). 
#   alpha  : array [n_genes, 2, n_samples] for transcription rate (alpha[,1,] = off, alpha[,2,] = on).
#   beta   : scalar or array [n_genes, n_samples] for splicing rate. By default equal to 1.
#   gamma   : scalar or array [n_genes, n_samples] for degradation rate.
#
# Returns:
#   Array of same shape as t0_off with s0 values for each element.
#   Note that we can compute just one between s0_off and s0_on (for different genes) for each call of this function. 
# -------------------------------------------------------------------------
s0_MCMC <- function(t0_off = NA, t0_on = 0, k, alpha, beta, gamma){
  # If multiple k values provided, keep only unique values.
  if(length(k) > 1){
    k <- unique(k)
  }

  # initial point for repressive phase
  if(k == 0){
    # transform s0_on according to the dimension of t0_off
    s0_on = alpha[, 1, ]/gamma
    s0_on <- array(rep(as.vector(s0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    s0_on <- aperm(s0_on, c(3, 1, 2))

    # transform s0_on according to the dimension of t0_off
    u0_on = alpha[,1,]/beta
    u0_on <- array(rep(as.vector(u0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    u0_on <- aperm(u0_on, c(3, 1, 2))
      
    t <- t0_off # we assume the switch from on to off phase occurs
    t0_on <- array(t0_on, dim = dim(u0_on))
    tau <- t - t0_on # elapsed time in the dynamic    

    # coordinate of the on steady state
    sSS_ON <- alpha[,2,]/gamma
    sSS_ON <- array(rep(as.vector(sSS_ON), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    sSS_ON <- aperm(sSS_ON, c(3, 1, 2))

    # transform the rates according to the dimension of t0:off
    beta <- array(rep(as.vector(beta), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    beta <- aperm(beta, c(3, 1, 2))
    gamma <- array(rep(as.vector(gamma), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    gamma <- aperm(gamma, c(3, 1, 2))
    alpha2 <- array(rep(as.vector(alpha[,2,]), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
    alpha2 <- aperm(alpha2, c(3, 1, 2))

    exp_gammaOFF <- exp(-gamma*tau)
    exp_betaOFF <- exp(-beta*tau)
    
    res <- s0_on*exp_gammaOFF + sSS_ON*(1-exp_gammaOFF)+ (alpha2-beta*u0_on)/(gamma-beta)*(exp_gammaOFF - exp_betaOFF)    
  }

  # initial point for inductive phase
  if(any(k == 2)){
      res <-  alpha[,1,]/gamma
  }  
  return(res)
}


# ------------- LOAD THE DATA ----------------

# -------------------------------------------------------------------------
# Function: nameGenes
# Purpose:
#   Load the name of the considered genes for the real data.
#
# Parameters:
#   pathData: path where the file (named as "var.csv") with the name of the genes is stored.

# Returns:
#   A vector with the name of the genes.
# -------------------------------------------------------------------------
nameGenes <- function(pathData){
    var <- read.csv(paste(pathData, "/var.csv", sep = ""))
    nameG <- var$index
    return(nameG)
}

# -------------------------------------------------------------------------
# Function: nameCells
# Purpose:
#   Load the name of the considered cells for the real data.
#
# Parameters:
#   pathData: path where the file (named as "obs.csv") with the name of the cells is stored.

# Returns:
#   A vector with the name of the cells.
# -------------------------------------------------------------------------
nameCells <- function(pathData){
    obs <- read.csv(paste(pathData, "/obs.csv", sep = ""))
    nameC<- obs$index
    return(nameC)
}

# -------------------------------------------------------------------------
# Function: loadRealData
# Purpose:
#   Load the values of spliced and unspliced RNA for the real data.
#
# Parameters:
#  typeProcessing: type of pre-processing that have been previously applied to the data.
#                  It can be "filter", "filter_and_normalize", "filter_and_normalize_noLog", filter_and_normalize moments"
#  pathData: path where the data are stored.
#  names: boolean if we want to keep the name of cells and genes or not. By default it is FALSE.
#  sparse: boolean if we want to transform the loaded matrices in sparse objects or not. By default it is FALSE.

# Returns:
#   A list with 
#     -"unspliced": matrix unspliced RNA values.
#     -"spliced"  : matrix spliced   RNA values.
#     -"typeCell" : vector with the group label associated to each cell. 
# -------------------------------------------------------------------------
loadRealData <- function(typeProcessing, pathData, names = FALSE, sparse = FALSE){
  if(is.null(pathData)){
    error("Provide the directory where the real data are stored")
  }
  
  pathData <- paste0(pathData, "/", typeProcessing)

  if(typeProcessing == "moments"){ # import pre-processed continuous data
    unspliced <- loadMu(pathData, names, sparse)
    spliced   <- loadMs(pathData, names, sparse)
  }else{  # import raw counts
    unspliced <- loadUnspliced(pathData, names, sparse)
    spliced   <- loadSpliced(pathData, names, sparse)
  }

  # import the type of cells 
  typeCell <- loadObs(pathData)$clusters

  res <- list("unspliced" = unspliced, "spliced" = spliced, "typeCell" = typeCell)
  return(res)
}

# -------------------------------------------------------------------------
# Function: loadVar
# Purpose:
#   Load the var.csv file, containing the name of the genes for the real data.
#
# Parameters:
#   pathData: path where the file (named as "var.csv") is stored.
#
# Returns: a matrix with the content of the var file.
# -------------------------------------------------------------------------
loadVar <- function(pathData){
  var <- read.csv(file = paste(pathData, '/var.csv', sep =""))
  return(var)
}

# -------------------------------------------------------------------------
# Function: loadObs
# Purpose:
#   Load the obs.csv file, containing the name of the cells and the type of cells'labels for the real data.
#
# Parameters:
#   pathData: path where the file (named as "obs.csv") is stored.
#
# Returns: a matrix with the content of the obs file.
# -------------------------------------------------------------------------
loadObs <- function(pathData){
  obs <- read.csv(file = paste(pathData, '/obs.csv', sep = "")) 
  return(obs)
}

# -------------------------------------------------------------------------
# Function: loadUnspliced
# Purpose:
#   Load the matrix with unspliced RNA counts.
#
# Parameters:
#  pathData: path where the file with unspliced counts is stored.
#  names: boolean if wee want to store the name of cells and genes. By default it is FALSE.
#  sparse: boolean if we want to transform the loaded matrix in sparse object or not. By default it is FALSE.
   
# Returns: a matrix with the content of the unspliced data.
# -------------------------------------------------------------------------  
loadUnspliced <- function(pathData, names = FALSE, sparse = FALSE){
  # import unspliced raw data and transform them into a sparse matrix, if required
  unspliced <- read.csv(file = paste(pathData, '/unspliced.csv', sep = ""))

  if(names){
    nameG <- nameGenes(pathData)
    nameC <- nameCells(pathData) 
    colnames(unspliced) <- nameG
    rownames(unspliced) <- nameC
  }
 
  unspliced <- as(unspliced, "matrix") 

  if(sparse){
    unspliced <- as(unspliced, "dgCMatrix")
  }

  return(unspliced)
}

# -------------------------------------------------------------------------
# Function: loadSpliced
# Purpose:
#   Load the matrix with spliced RNA counts.
#
# Parameters:
#  pathData: path where the file with spliced counts is stored.
#  names: boolean if wee want to store the name of cells and genes. By default it is FALSE.
#  sparse: boolean if we want to transform the loaded matrix in sparse object or not. By default it is FALSE.
   
# Returns: a matrix with the content of the spliced data.
# -------------------------------------------------------------------------   
loadSpliced <- function(pathData, names = FALSE, sparse = FALSE){
  # import spliced raw data and transform them into a sparse matrix, if required
  spliced <- read.csv(file = paste(pathData, '/spliced.csv', sep = ""))

  if(names){
    nameG <- nameGenes(pathData)
    nameC <- nameCells(pathData) 
    colnames(unspliced) <- nameG
    rownames(unspliced) <- nameC
  }
 
  spliced <- as(spliced, "matrix") 

  if(sparse){
    spliced <- as(spliced, "dgCMatrix")
  }

  return(spliced)
}

# -------------------------------------------------------------------------
# Function: loadMu
# Purpose:
#   Load the matrix with unspliced RNA moments (previously computed by scVelo)
#
# Parameters:
#  pathData: path where the file with unspliced moments is stored.
#  names: boolean if wee want to store the name of cells and genes. By default it is FALSE.
#  sparse: boolean if we want to transform the loaded matrix in sparse object or not. By default it is FALSE.
   
# Returns: a matrix with the the unspliced moments.
# -------------------------------------------------------------------------   
loadMu <- function(pathData, names = FALSE, sparse = FALSE){
  # import unspliced pre-processed moments and transform them into a sparse matrix, if required
  Mu <- read.csv(file = paste(pathData, '/Mu.csv', sep = ""))

  if(names){
    nameG <- nameGenes(pathData)
    nameC <- nameCells(pathData) 
    colnames(Mu) <- nameG
    rownames(Mu) <- nameC
  }
 
  Mu <- as(Mu, "matrix") 

  if(sparse){
    Mu <- as(Mu, "dgCMatrix")
  }

  return(Mu)
}

# -------------------------------------------------------------------------
# Function: loadMs
# Purpose:
#   Load the matrix with spliced RNA moments (previously computed by scVelo)
#
# Parameters:
#  pathData: path where the file with spliced moments is stored.
#  names: boolean if wee want to store the name of cells and genes. By default it is FALSE.
#  sparse: boolean if we want to transform the loaded matrix in sparse object or not. By default it is FALSE.
   
# Returns: a matrix with the the spliced moments.
# -------------------------------------------------------------------------   
loadMs <- function(pathData, names = FALSE, sparse = FALSE){
  # import spliced pre-processed moments and transform them into a spare matrix, if required
  Ms <- read.csv(file = paste(pathData, '/Ms.csv', sep = ""))

  if(names){
    nameG <- nameGenes(pathData)
    nameC <- nameCells(pathData) 
    colnames(Ms) <- nameG
    rownames(Ms) <- nameC
  }
 
  Ms <- as(Ms, "matrix") 

  if(sparse){
    Ms <- as(Ms, "dgCMatrix")
  }

  return(Ms)
}

# -------------------------------------------------------------------------
# Function: nameReal_Pancreas
# Purpose:
#   Convert the name of thee simulations into the names used in Table 5 of bayVel's paper
#
# Parameters:
#  name: name of the considered simulation
   
# Returns: string with the corresponding name used in Table 5. 
# -------------------------------------------------------------------------   
nameReal_Pancreas <- function(name){
  label <- ""
  if(grepl("SW1", name)){
      label <- paste0(label,"K = 1,$")
  }else{
      label <- paste0(label,"K = 8,$")
  }

  if(grepl("T1", name)){
      label <- paste(label,"$R = 1$", sep = "-")
  }else if(grepl("T2", name)){
      label <- paste(label,"$R = 9$", sep = "-")
  }else if(grepl("T3", name)){
      label <- paste(label,"$R = 38$", sep = "-")
  }
  
  return(label)
}


# ------------- PLOT THE DATA ----------------

# -------------------------------------------------------------------------
# Function: plot_sVSu
# Purpose:
#   Plot the phase trajectory (s vs u) for a single gene using the
#   analytical solutions of the ODE model (functions `u()` and `s()`).
#   The function draws both the dynamic when the switch occurs and the 
#   potential dynamic when the switch does not occur. The function admit to 
#   overlay observed points for cell types.
#
# Arguments:
#   t0_off        : numeric scalar or Inf, time of the OFF switch (if Inf,
#                   the system stays in the ON steady-state before switching).
#   t0_on         : numeric scalar (default 0), time when ON state starts.
#   alpha         : transcription-rate vector.
#   beta          : splicing rate, equal to 1 by default.
#   gamma         : degradation rate.
#   pos_u         : vector of u positions (one per cell/subgroup). Omit if you want just the model dynamic.
#   pos_s         : vector of s positions (one per cell/subgroup). Omit if you want just the model dynamic.
#   subGrLabels   : vector of subgroup-labels for each cell.
#   g             : gene identifier (used in the plot title).
#   add           : boolean; if FALSE create a new ggplot, if TRUE add to `gg`. By default it is equal to FALSE.
#   gg            : an existing ggplot object to add layers to when add = TRUE.
#   colCell       : colors for subgroup points (vector or NA to auto-pick).
#   colDyn        : color for the model dynamic (optional, default "red").
#   xlim, ylim    : numeric scalars for axis limits (optional, if NA they are not set).
#   axisTitle.size, axisText.size, title.size : sizes for theme elements.
#   lineSize      : line width for model dynamic.
#   shapePoint    : shape of subgroup points (optional).
#   sizePoint     : size  of subgroup points (optional).
#
# Returns:
#   A ggplot object with the s vs u phase plot for the given gene.
# -------------------------------------------------------------------------
plot_sVSu <- function(
                      t0_off, 
                      t0_on = 0, 
                      alpha, 
                      beta = 1, 
                      gamma, 
                      pos_u, 
                      pos_s, 
                      subGrLabels, 
                      g, 
                      add = FALSE, 
                      gg = NA, 
                      colCell = NA, 
                      colDyn = "red", 
                      xlim = NA, 
                      ylim = NA, 
                      axisTitle.size = NA, 
                      axisText.size = NA, 
                      title.size = NA, 
                      lineSize = 1, 
                      shapePoint = 21, 
                      sizePoint = 2, 
                      ...
                      ){

  
  # Compute the initial points of the on and of the off phase
  u0_off <- u0(k = 0, alpha = alpha, beta = beta)
  u0_on <- u0(k = 2, alpha = alpha, beta = beta)  
  s0_off <- s0(k = 0, alpha = alpha, beta = beta, gamma = gamma)
  s0_on <- s0(k = 2, alpha = alpha, beta = beta, gamma = gamma)
 
  # Compute the points lying on the on branch of the dynamic. We assume here the switching point is infinite, in order to plot the potential behavior up to the upper steady state.
  t_seq <- seq(t0_on_real, 500, 0.1)
  u_plot <- u(t_seq, t0_off = Inf, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(t_seq)), alpha, beta)
  s_plot <- s(t_seq, t0_off = Inf, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(t_seq)), alpha, beta, gamma)
  
  df_on <- data.frame(s_plot = s_plot, u_plot = u_plot)

  # Title for the plot
  main <- paste("Gene", g)
  if(length(unique(subGrLabels)) == 1){
    main <- paste(main, ", typeCell", unique(subGrLabels)) # we are plotting just the position of one subgroup
  }

  # Initialize or extend ggplot with the ON branch (dashed line)
  if(!add){
    gg <- ggplot(data = df_on, aes(x = s_plot, y = u_plot)) + 
          geom_path(lty = 2, col = colDyn, linewidth = lineSize) + 
          xlab("s") + 
          ylab("u") + 
          labs(title = main)
  }else{
    gg <- gg + 
          geom_path(data = df_on, aes(x = s_plot, y = u_plot), lty = 2, col = colDyn, linewidth = lineSize)
  }

  # adjust graphical parameters
  if(is.numeric(axisText.size)){
    gg <- gg + theme(axis.text = element_text(size = axisText.size))
  }
  if(is.numeric(axisTitle.size)){
    gg <- gg + theme(axis.title = element_text(size = axisTitle.size))
  }
  if(is.numeric(title.size)){
    gg <- gg + theme(plot.title =  element_text(hjust = 0.5, size = title.size))
  }
  if(is.numeric(xlim)){
    gg <- gg + xlim(0, xlim) 
  }
  if(is.numeric(ylim)){
    gg <- gg + ylim(0, ylim) 
  }

  
  # Compute the points lying on the off branch of the dynamic, describing the potential behavior from to the upper steady state.
  u_plot <- u(t_seq, t0_off = 0, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(t_seq)), alpha, beta)
  s_plot <- s(t_seq, t0_off = 0, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t_seq)), alpha, beta, gamma)
  
  df_off <- data.frame(s_plot = s_plot, u_plot = u_plot)
  gg <- gg + geom_path(data = df_off, aes(x = s_plot, y = u_plot), lty = 2, col = colDyn, linewidth = lineSize) 
  
  # If a finite t0_off_real is provided, compute the transient branches:
  #  - ON branch from t0_on_real to t0_off_real (solid line)
  #  - OFF branch from t0_off_real to t0_off_real + 500 (solid line)
  # Else (t0_off_real == Inf) we reuse the previously computed branch(s).

  if(t0_off != Inf){
    # compute initial conditions at the switching times
    u0_off <- u0(t0_off = t0_off, k = 0, alpha = alpha, beta = beta)
    u0_on  <- u0(t0_on = t0_on, k = 2, alpha = alpha, beta = beta)
    s0_off <- s0(t0_off = t0_off, k = 0, alpha = alpha, beta = beta, gamma = gamma)
    s0_on  <- s0(t0_on = t0_on, k = 2, alpha = alpha, beta = beta, gamma = gamma)

 
    # ON transient: from t0_on_real up to t0_off_real (solid line)
    t_on_trans <- seq(t0_on, t0_off, 0.1)
    u_plot <- u(t_on_trans, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(t_on_trans)), alpha, beta)
    s_plot <- s(t_on_trans, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(t_on_trans)), alpha, beta, gamma)

    df_on_trans <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + 
          geom_path(data = df_on_trans, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 
   
    #  OFF transient: from t0_off_real onward (solid line)
    t_off_trans <- seq(t0_off, t0_off + 500, 0.1)
    u_plot <- u(t_off_trans, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(t_off_trans)), alpha, beta)
    s_plot <- s(t_off_trans, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t_off_trans)), alpha, beta, gamma)

    df_off_trans <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + 
          geom_path(data = df_off_trans, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 
  }else{
    # If infinite off-time: overlay the earlier computed branches.
    gg <- gg + 
          geom_path(data = df_on,  aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) + 
          geom_path(data = df_off, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 
  }
  
  # If pos_* vectors are provided (no NAs) add points (pos_s, pos_u) for each subgroup
  # Colors are auto-chosen if colCell is NA.
  subGrC <- unique(subGrLabels)

  if(sum(is.na(pos_s)) == 0){
    # If lengths already match unique types, use directly; otherwise pick
    # the first observed position for each subtype.
    if(length(pos_s) == length(subGrC)){
      df <- data.frame(s = pos_s, u = pos_u)
    }else{
      pos_s <- c()
      pos_u <- c()
      for(i in 1:length(subGrC)){
          x <- subGrC[i]
          pos_s <- c(pos_s, pos_s[which(subGrLabels == x)[1]])
          pos_u <- c(pos_u, pos_u[which(subGrLabels == x)[1]]) 
      }
      df <- data.frame(s = pos_s, u = pos_u)
    }

    # Choose palette if colCell not provided
    if(sum(is.na(colCell)) == 1){
      if(length(subGrC) >= 3){
        colCell <- brewer.pal(length(subGrC), "Spectral")
      }else if(length(subGrC) == 2){
        colCell <- c("darkgreen", "blue")
      }else if(length(subGrC) == 1){
        colCell <- c("darkgreen")
      }
    }
    df$colCell <- colCell
    # Add points to the plot
    gg <- gg + geom_point(data = df, aes(x = s, y = u, fill = colCell), fill = colCell, col = colDyn, size = sizePoint, shape = shapePoint)
  }

  return(gg)
}

# -------------------------------------------------------------------------
# Function: plot_sVSu_scVelo
# Purpose:
#   Plot the scVelo-style s vs u phase trajectory for a single gene using
#   analytical ODE solutions (functions u() and s()). It takes as input 
#   the output of scVelo's inference. The function draws:
#     - ON and OFF branches (dashed)
#     - transient branches (solid) if fit_t0_off is finite
#     - cell positions computed from fit_t (scaled/translated to model fit)
#
# Notes:
#   - scVelo uses beta_scaled = fit_beta * fit_scaling; alpha_off is assumed 0 and 
#     alpha_on = fit_alpha
#   - Trajectories are rescaled and translated with fit_scaling, fit_u0_offset,
#     fit_s0_offset to match scVelo plotting conventions.
#
# Arguments:
#   fit_t0_off      : numeric scalar or Inf, fitted switch time to OFF state.
#   fit_t0_on       : numeric scalar (default 0), fitted ON start time.
#   fit_alpha       : on transcription rate.
#   fit_beta        : splicing rate.
#   fit_gamma       : degradation rate-
#   fit_scaling     : scalar used to scale u (beta is multiplied by this).
#   fit_u0_offset   : scalar added to u after scaling (translation).
#   fit_s0_offset   : scalar added to s (translation).
#   fit_t           : numeric vector of per-cell times used to place cells on the curve.
#   subGrLabels     : vector of group labels for each element of fit_t (used for colors).
#   g               : gene identifier used in the plot title.
#   add             : boolean; if FALSE create new ggplot, if TRUE add to gg.
#   gg              : existing ggplot object to add layers to (when add = TRUE).
#   colCell         : color for points (vector or NA to auto-pick).
#   colDyn          : color for the model dynamic and points border (optional, default "red").
#   xlim, ylim      : numeric scalars for axis limits (optional, if NA they are not set).
#   axisTitle.size, axisText.size, title.size : sizes for theme elements.
#   lineSize        : line width for model dynamic.
#   shapePoint      : shape of points (optional).
#   sizePoint       : size  of points (optional). 
#
# Returns:
#   A ggplot object with the s vs u phase plot for the given gene, accordingg to scVelo plots.
# -------------------------------------------------------------------------
plot_sVSu_scVelo <- function(
                              fit_t0_off, 
                              fit_t0_on = 0, 
                              fit_alpha, 
                              fit_beta, 
                              fit_gamma, 
                              fit_scaling, 
                              fit_u0_offset, 
                              fit_s0_offset, 
                              fit_t, 
                              subGrLabels,
                              g, 
                              add = FALSE, 
                              gg = NA, 
                              colCell = NA, 
                              colDyn = "red", 
                              xlim = NA, 
                              ylim = NA, 
                              axisTitle.size = NA, 
                              axisText.size = NA, 
                              title.size = NA, 
                              lineSize = 1, 
                              shapePoint = 21, 
                              sizePoint = 1, 
                              ...
                            ){

  # Prepare kinetic parameters according to scVelo convention
  # scVelo uses beta_scaled = fit_beta * fit_scaling
  # alpha is c(0, fit_alpha) so that alpha[1] = 0 (off), alpha[2] = fit_alpha (on)  
  beta <- fit_beta * fit_scaling 
  alpha <- c(0, fit_alpha)

  # Compute the initial points of the on and of the off phase
  u0_off <- u0(k = 0, alpha = alpha, beta = beta)
  u0_on  <- u0(k = 2, alpha = alpha, beta = beta)
  s0_off <- s0(k = 0, alpha = alpha, beta = beta, gamma = fit_gamma)
  s0_on  <- s0(k = 2, alpha = alpha, beta = beta, gamma = fit_gamma)
 
  # Compute the points lying on the on branch of the dynamic. We assume here the switching point is infinite, in order to plot the potential behavior up to the upper steady state.
  t_on_seq <- seq(fit_t0_on, 500, 0.1)
  u_plot <- u(t_on_seq, t0_off = Inf, t0_on = 0, u0_off = u0_off, u0_on = 0, k = rep(2, length(t_on_seq)), alpha, beta)
  s_plot <- s(t_on_seq, t0_off = Inf, t0_on = 0, u0_off = u0_off, u0_on = 0, s0_off = s0_off, s0_on = 0, k = rep(2, length(t_on_seq)), alpha, beta, fit_gamma)

  # Apply scaling and translation to match scVelo plotting
  u_plot <- u_plot * fit_scaling + fit_u0_offset
  s_plot <- s_plot + fit_s0_offset
  
  df_on <- data.frame(s_plot = s_plot, u_plot = u_plot)

  # Title for the plot
  main <- paste("Gene", g)
  if(length(unique(subGrLabels))== 1){
    main <- paste(main, ", typeCell", unique(subGrLabels))
  }

  # Initialize or extend ggplot with the ON branch (dashed line)
  if(!add){
    gg <- ggplot(data = df_on, aes(x = s_plot, y = u_plot)) + 
          geom_path(lty = 2, col = colDyn, linewidth = lineSize) + 
          xlab("s") + 
          ylab("u") + 
          labs(title = main) 
          
    # adjust graphical parameters
    if(is.numeric(axisText.size)){
      gg <- gg + theme(axis.text = element_text(size = axisText.size))
    }
    if(is.numeric(axisTitle.size)){
      gg <- gg + theme(axis.title = element_text(size = axisTitle.size))
    }
    if(is.numeric(title.size)){
      gg <- gg + theme(plot.title =  element_text(hjust = 0.5, size = title.size))
    }
    if(is.numeric(xlim)){
      gg <- gg + xlim(0, xlim) 
    }
    if(is.numeric(ylim)){
      gg <- gg + ylim(0, ylim) 
    }
  }else{
    gg <- gg + 
          geom_path(data = df_on, aes(x = s_plot, y = u_plot), lty = 2, col = colDyn, linewidth = lineSize) 
  }

  # Compute the points lying on the off branch of the dynamic, describing the potential behavior from to the upper steady state.
  t_off_seq <- seq(0, 5000, 0.05)
  u_plot <- u(t_off_seq, t0_off = 0, t0_on = 0, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(t_off_seq)), alpha, beta)
  s_plot <- s(t_off_seq, t0_off = 0, t0_on = 0, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t_off_seq)), alpha, beta, fit_gamma)

  # Apply scaling and translation to match scVelo plotting
  u_plot <- u_plot * fit_scaling + fit_u0_offset
  s_plot <- s_plot + fit_s0_offset

  df_off <- data.frame(s_plot = s_plot, u_plot = u_plot)
  gg <- gg + geom_path(data = df_off, aes(x = s_plot, y = u_plot), lty = 2, col = colDyn, linewidth = lineSize) 
  
  # If a finite t0_off_real is provided, compute the transient branches:
  #  - ON branch from fit_t0_on to fit_t0_off (solid line)
  #  - OFF branch from fit_t0_off to fit_t0_off + 5000 (solid line)
  # Else (t0_off_real == Inf) we reuse the previously computed branch(s).
  if(fit_t0_off != Inf){
    # compute initial conditions at the switching times
    u0_off <- u0(t0_off = fit_t0_off, k = 0, alpha = alpha, beta = beta)
    u0_on <- u0(t0_on = fit_t0_on, k = 2, alpha = alpha, beta = beta)
    s0_off <- s0(t0_off = fit_t0_off, k = 0, alpha = alpha, beta = beta, gamma = fit_gamma)
    s0_on <- s0(t0_on = fit_t0_on, k = 2, alpha = alpha, beta = beta, gamma = fit_gamma)

    # ON transient: from fit_t0_on up to fit_t0_off (solid line)
    t_on_trans <- seq(fit_t0_on, fit_t0_off, 0.01)
    u_plot <- u(t_on_trans, t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(t_on_trans)), alpha, beta)
    s_plot <- s(t_on_trans, t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(t_on_trans)), alpha, beta, fit_gamma)
    # Apply scaling and translation to match scVelo plotting
    u_plot <- u_plot * fit_scaling + fit_u0_offset
    s_plot <- s_plot + fit_s0_offset

    df_on_trans <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + geom_path(data = df_on_trans, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 

    #  OFF transient: from fit_t0_off onward (solid line)
    t_off_trans <- seq(fit_t0_off, fit_t0_off + 5000, 0.01)
    u_plot <- u(t_off_trans, t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(t_off_trans)), alpha, beta)
    s_plot <- s(t_off_trans, t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t_off_trans)), alpha, beta, fit_gamma)
    # Apply scaling and translation to match scVelo plotting
    u_plot = u_plot * fit_scaling + fit_u0_offset
    s_plot = s_plot + fit_s0_offset

    df_off_trans <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + 
          geom_path(data = df_off_trans, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 
  }else{
    # If infinite off-time: overlay the earlier computed branches.
    gg <- gg + 
          geom_path(data = df_on, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) + 
          geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colDyn, linewidth = lineSize) 
  }

  # If fit_t is provided (no NAs) add points (pos_s, pos_u) 
  # Colors are auto-chosen if colCell is NA.
  n_cells <- length(fit_t)
  pos_u <- rep(NA, n_cells)
  pos_s <- rep(NA, n_cells)
  k <- rep(0, n_cells)
  k[which(fit_t < fit_t0_off)] <- 2
  pos_u <- u(fit_t, t0_off = fit_t0_off, t0_on = 0, u0_off = u0_off, u0_on = u0_on, k = k, alpha = alpha, beta = beta)
  pos_s <- s(fit_t, t0_off = fit_t0_off, t0_on = 0, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = k, alpha = alpha, beta = beta, gamma = fit_gamma)
  # Apply scaling and translation to match scVelo plotting
  pos_u <- pos_u * fit_scaling + fit_u0_offset
  pos_s <- pos_s + fit_s0_offset

  df <- data.frame(s = pos_s, u = pos_u)

  # Choose palette if colCell not provided
  subGrC <- unique(subGrLabels)
  if(sum(is.na(colCell)) == 1){
    colCell <- brewer.pal(length(subGrC), "Spectral")
  }
  color <- c()
  for(i in 1:length(subGrC)){
    x <- subGrC[i]
    color <- c(color, rep(colCell[i], length(which(subGrLabels == x))))
  }
  if(nrow(df) == length(color)){
    df$colCell <- color
  }else{
    df$colCell <- color[1:nrow(df)]
  }
  # Add points to the plot
  gg <- gg + geom_point(data = df, aes(x = s, y = u, fill = colCell), fill = df$colCell, col = colDyn, size = sizePoint, shape = shapePoint)

  return(gg)
}

# -------------------------------------------------------------------------
# Function: plot_GeneDynamic_withNotes
# Purpose:
#   Draw the s-vs-u phase diagram for a single gene with annotation for Figure 1
#
# Arguments:
#   t0_off : numeric scalar; time of OFF switch (use Inf if no switch before SS)
#   u0_off : numeric scalar; u-coordinate of the switching point (u^omega)
#   s0_off : numeric scalar; s-coordinate of the switching point (s^omega)
#   t0_on  : numeric scalar; time when ON phase starts (default 0)
#   alpha  : numeric vector length 2 with transcription rates: c(alpha_off, alpha_on)
#   beta   : numeric scalar for splicing rate
#   gamma  : numeric scalar fore degradation rate
#   r      : optional label used in title (default NA)
#   g      : gene identifier (used in title; may be NA)
#
# Returns:
#   A ggplot object (s vs u) with annotations.
# -------------------------------------------------------------------------
plot_GeneDynamic_withNotes <- function(
                                        t0_off, 
                                        u0_off, 
                                        s0_off, 
                                        t0_on = 0, 
                                        alpha, 
                                        beta, 
                                        gamma, 
                                        r = NA, 
                                        g
                                      ){
  
  # compute the coordinate of the steady states
  u_SS_off <- alpha[1]/beta
  u_SS_on  <- alpha[2]/beta
  s_SS_off <- alpha[1]/gamma        
  s_SS_on  <- alpha[2]/gamma        

  # Compute the points lying on the on branch of the dynamic. We assume here the switching point is infinite, in order to plot the potential behavior up to the upper steady state.
  u_plotON <- seq(u_SS_off, u_SS_on, 0.01)
  s_plotON <- s_di_u(u_plotON, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 2, alpha, beta, gamma)

  # Compute the points lying on the off branch of the dynamic, describing the potential behavior from to the upper steady state.
  u_plotOFF <- seq(u_SS_off, u_SS_on, 0.01)
  s_plotOFF <- s_di_u(u_plotOFF, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 0, alpha, beta, gamma)

  df <- data.frame(uON_SS = as.vector(u_plotON), sON_SS = as.vector(s_plotON), uOFF_SS = as.vector(u_plotOFF), sOFF_SS = as.vector(s_plotOFF), gene = g)
  # Plot on and off branches
  pl <- ggplot(df)  +
    geom_path(aes(x = sON_SS, y = uON_SS, color = "red"), linetype = "dotted", linewidth = 3) + 
    geom_path(aes(x = sOFF_SS, y = uOFF_SS, color = "blue"), linetype = "dotted", linewidth = 3)

  # If a finite t0_off_real is provided, compute the transient branches.
  if(t0_off != Inf){
    # ON transient: from the lower steady state to the switching point.
    u_plotON_switch <- seq(u_SS_off, u0_off, 0.001)
    s_plotON_switch <- s_di_u(u_plotON_switch, u0_off, s0_off, u_SS_off, s_SS_off, k = 2, alpha, beta, gamma)

    # OFF transient: from the switching point to the lower steady state
    u_plotOFF_switch <- seq(u_SS_off, u0_off, 0.001)
    s_plotOFF_switch <- s_di_u(u_plotOFF_switch, u0_off, s0_off, u_SS_off, s_SS_off, k = 0, alpha, beta, gamma)

    dfSwitchON <- data.frame(uON_switch = u_plotON_switch, sON_switch = s_plotON_switch)
    dfSwitchOFF <- data.frame(uOFF_switch = u_plotOFF_switch, sOFF_switch = s_plotOFF_switch)

    maxSoff <- which.max(s_plotOFF_switch)

    # add the two switching branches
    pl <- pl +
      geom_path(data = dfSwitchON, aes(x = sON_switch, y = uON_switch, color = "red"), linewidth = 3) +     
      geom_path(data = dfSwitchOFF, aes(x = sOFF_switch, y = uOFF_switch, color = "blue"), linewidth = 3)   

  }else{
    # If infinite off-time: overlay the earlier computed branches.
    pl <- pl +
      geom_path(aes(x = sON_SS, y = uON_SS, color = "red"), linewidth = 3) +
      geom_path(aes(x = sOFF_SS, y = uOFF_SS, color = "blue"), linewidth = 3)
  }

  # Title for the plot
  if(is.na(g)){
    g <- "g"
  }
  main <- ""
  if(!is.na(r)){
    main <- paste(main, " for group r")
  }

 # Annotations: steady-state and switching point labels
  pl <- pl + 
        theme(legend.position = "none") + 
        annotate(geom = "text", y = alpha[1]*0.88, x = (alpha[1]/gamma)*0.9,  label = TeX("$SS^{off}$", output = "character"), size = 35, parse = TRUE,  family = "serif") +
        annotate(geom = "text", y = alpha[2]*1.05, x = (alpha[2]/gamma)*1.02, label = TeX("$SS^{on}$",  output = "character"), size = 35, parse = TRUE,  family = "serif") +
        annotate(geom = "text", y = u0_off*1.05, x = s0_off*0.92,
        label =TeX(r"($(s^{omega}, u^{omega})$)", output = "character"), 
        size = 35, parse = TRUE, family = "serif") +
        annotate(geom = "text", x = min(df$sON_SS) + min(df$sON_SS)*0.45, y = min(df$uON_SS) + (max(df$uON_SS) - min(df$uON_SS))/2, label = "Induction", size = 35, angle = 62, family = "serif")  +
        annotate(geom = "text", x = max(df$sON_SS)*0.82, y = min(df$uON_SS) + (max(df$uON_SS) - min(df$uON_SS))/2 -0.1, label = "Repression", size = 35, angle = 56, family = "serif") + 
        theme(plot.title = element_text(family = "serif", size=70,  hjust = 0.5)) + 
        labs(x = "", y = "") 

  # add arrows for direction of transient induction branch
  ind_on_trans <- which.min(abs(dfSwitchON$uON_switch - (min(dfSwitchON$uON_switch) + (max(dfSwitchON$uON_switch) -  min(dfSwitchON$uON_switch))/2)))  
  arr_on_u1 <-  dfSwitchON$uON_switch[ind_on_trans]
  arr_on_s1 <- dfSwitchON$sON_switch[which(dfSwitchON$uON_switch == arr_on_u1)[1]]
  arr_on_u2 <- dfSwitchON$uON_switch[ind_on_trans + 1]
  arr_on_s2 <- dfSwitchON$sON_switch[which(dfSwitchON$uON_switch == arr_on_u2)[1]]

  pl <- pl  + 
        geom_segment(aes(color = "red"), x = arr_on_s1, y = arr_on_u1, xend = arr_on_s2, yend = arr_on_u2, arrow = arrow( length = unit(0.3, "inches")), size = 3)

  # add arrows for direction of transient repressive branch
  ind_off_trans <- which.min(abs(dfSwitchOFF$uOFF_switch - (min(dfSwitchOFF$uOFF_switch) + (max(dfSwitchOFF$uOFF_switch) -  min(dfSwitchOFF$uOFF_switch))/2)))
  arr_off_u1 <-  dfSwitchOFF$uOFF_switch[ind_off_trans]
  arr_off_s1 <- dfSwitchOFF$sOFF_switch[which(dfSwitchOFF$uOFF_switch == arr_off_u1)[1]]
  arr_off_u2 <- dfSwitchOFF$uOFF_switch[ind_off_trans + 1]
  arr_off_s2 <- dfSwitchOFF$sOFF_switch[which(dfSwitchOFF$uOFF_switch == arr_off_u2)[1]]

  pl <- pl + geom_segment(aes(color = "blue"), x = arr_off_s2, y = arr_off_u2, xend = arr_off_s1, yend = arr_off_u1, arrow = arrow( length = unit(0.3, "inches")), size = 3)       

  # add arrows for direction of not-transient repressive branch
  ind_off <- which.min(abs(df$uOFF_SS - (min(df$uOFF_SS) + (max(df$uOFF_SS) -  min(df$uOFF_SS))/2)))
  arr_off_u1 <- df$uOFF_SS[ind_off]     
  arr_off_s1 <- df$sOFF_SS[which(df$uOFF_SS == arr_off_u1)[1]]
  arr_off_u2 <- df$uOFF_SS[ind_off + 5]  
  arr_off_s2 <- df$sOFF[which(df$uOFF_SS == arr_off_u2)[1]]

  pl <- pl + geom_segment(aes(colour = "blue"), x = arr_off_s2, y = arr_off_u2, xend = arr_off_s1, yend = arr_off_u1,  arrow = arrow(length = unit(0.3, "inches")), size = 3)

  pl <- pl +
    coord_cartesian(xlim =  c(min(rbind(df$sON_SS, df$sOFF_SS)) - 0.5, max(rbind(df$sON_SS, df$sOFF_SS)) + 0.5), ylim = c(min(rbind(df$uON_SS, df$uOFF_SS)) - 0.2, max(rbind(df$uON_SS, df$uOFF_SS)) + 0.2)) + theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

  # draw cartesian axes on the side
  pl <- pl + geom_segment(aes(x=3.7, y=0.85, xend=4.5, yend=0.85), arrow = arrow(length=unit(.5, 'cm')), color='black', linewidth=2) + 
  geom_segment(aes(x=3.7, y=0.846, xend=3.7, yend=1.35), arrow = arrow(length=unit(.5, 'cm')), color='black', linewidth=2) + 
  annotate(geom = "text", y = 0.78, x = 4.1,  label = TeX("$s$", output = "character"), size = 30, parse = TRUE,  family = "serif") +
  annotate(geom = "text", y = 1.046, x = 3.55,  label = TeX("$u$", output = "character"), size = 30, parse = TRUE,  family = "serif") 

  return(pl)
}


