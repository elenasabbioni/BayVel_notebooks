# Input: the result of rmultinom(n_cells, 1, pk): it has 4 rows and n_cells columns
# We want to transform it in a 1 x n_cells array, with each element denoting the state
# k (from 0 to 3)
assign_state <- function(v){
  n_cells <- dim(v)[2]
  n_states <- dim(v)[1]
  res <- rep(NA, n_cells)
  for(i in 1:n_cells){
    res[i] <- which(v[,i]==1) - 1
  }
  return(res)
}

# unspliced RNA (soluzione ODE) per n_cells (1 gene fixed)
# tau: vettore n_cells dimensionale
# k: vettore n_cells dimensionale
# alpha: vettore bidimensionale (2 stati)
# beta: scalare
u <- function(t, t0_off, t0_on = 0, u0_off, u0_on = NA, k, alpha, beta){
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1 (off), alpha[2] = alpha_2/3 (on)
  if(any(is.na(u0_on))){
    u0_on = alpha[1]/beta
  }
  n_cells <- length(k)
  
  tau <- ifelse(k == 2, t - t0_on, t - t0_off)
  expBetaOFF <- exp(-beta*tau[k == 0])
  expBetaON <- exp(-beta*tau[k == 2])
  res <- rep(NA, n_cells)
  
  u0 <- ifelse(k == 2, u0_on, u0_off)

  res[k == 0] <- u0[k == 0]*expBetaOFF + alpha[1]/beta*(1-expBetaOFF)
  res[k == 2] <- u0[k == 2]*expBetaON  + alpha[2]/beta*(1-expBetaON)
  
  return(res)
}


# spliced RNA (soluzione ODE) per n_cells (1 gene fixed)
# tau: vettore n_cells dimensionale
# k: vettore n_cells dimensionale
# alpha: vettore bidimensionale (2 stati)
# beta: scalare
# gamma: scalare
s <- function(t, t0_off, t0_on = 0, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma){
  
  if(is.na(u0_on)){
    u0_on = alpha[1]/beta
  }
  if(is.na(s0_on)){
    s0_on = alpha[1]/gamma
  }
  n_cells <- length(k)
  
  tau <- ifelse(k == 2, t - t0_on, t - t0_off)
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
  exp_gammaOFF <- exp(-gamma*tau[k == 0])
  exp_gammaON <- exp(-gamma*tau[k == 2])
  exp_betaOFF <- exp(-beta*tau[k == 0])
  exp_betaON <- exp(-beta*tau[k == 2])
  
  u0 <- ifelse(k == 2, u0_on,  u0_off)
  s0 <- ifelse(k == 2, s0_on,  s0_off)
  res <- rep(NA, n_cells)
  
  res[k == 0] <- s0[k == 0]*exp_gammaOFF + alpha[1]/gamma*(1-exp_gammaOFF)+ (alpha[1]-beta*u0[k == 0])/(gamma-beta)*(exp_gammaOFF - exp_betaOFF)
  res[k == 2] <- s0[k == 2]*exp_gammaON  + alpha[2]/gamma*(1-exp_gammaON) + (alpha[2]-beta*u0[k == 2])/(gamma-beta)*(exp_gammaON  - exp_betaON)
  return(res)
}


u_withTau <- function(tau, u0_off, u0_on = NA, k, alpha, beta){
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1 (off), alpha[2] = alpha_2/3 (on)
  if(is.na(u0_on)){
    u0_on = alpha[1]/beta
  }
  n_cells <- length(k)
  
  
  expBetaOFF <- exp(-beta*tau[k == 0])
  expBetaON <- exp(-beta*tau[k == 2])
  res <- rep(NA, n_cells)
  
  u0 <- ifelse(k == 2, u0_on, u0_off)

  res[k == 0] <- u0[k == 0]*expBetaOFF + alpha[1]/beta*(1-expBetaOFF)
  res[k == 2] <- u0[k == 2]*expBetaON  + alpha[2]/beta*(1-expBetaON)
  
  return(res)
}

# calcola u per un gene e per una cellula fissati, per tutte le iterazioni dell'MCMC
# passiamo, per ogni parametro una catena con il valore di questo parametro nelle varie iterazioni
u_withTauMCMC <- function(tau, u0_off, u0_on = NA, k, alpha, beta){
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1 (off), alpha[2] = alpha_2/3 (on)
  if(any(is.na(u0_on))){
    u0_on = alpha[1,]/beta
  }
  sampleToSave <- length(tau)
  
  
  expBeta <- exp(-beta*tau)
 
  res <- rep(NA, sampleToSave)
  
  u0 <- ifelse(k == 2, u0_on, u0_off)
  res <- u0*expBeta + alpha[1,]/beta*(1-expBeta)*c(k==0) + alpha[2,]/beta*(1-expBeta)*c(k==2) 
  
  return(res)
}


# u_withTauMCMC_new <- function(tau, u0_off, u0_on = NA, k, alpha, beta, subtypeCell, typeCellT0_off){
#   u0_off_new <- array(NA, dim = dim(tau))

#   for(sty in unique(subtypeCell)){
#     tyT0_off <- typeCellT0_off[which(subtypeCell == sty)[1]]
#     u0_off_new[sty, , ] = u0_off[tyT0_off, , ]
#   }

#   uSS_on <- alpha[,2,]/beta
#   uSS_on <- array(rep(as.vector(uSS_on),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   uSS_on <- aperm(uSS_on, c(3, 1, 2))

#   u0_on <- alpha[, 1, ]/beta
#   u0_on <- array(rep(as.vector(u0_on),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   u0_on <- aperm(u0_on, c(3, 1, 2))

#   u0 <- ifelse(k == 2, u0_on, u0_off_new)

#   beta <- array(rep(as.vector(beta),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   beta <- aperm(beta, c(3, 1, 2))
#   expBeta <- exp(-beta*tau)
 
#   res <- u0*expBeta + u0_on*(k==0)*(1-expBeta) + uSS_on*(k==2)*(1-expBeta)

#   return(res)
# }

# spliced RNA (soluzione ODE) per n_cells (1 gene fixed)
# tau: vettore n_cells dimensionale
# k: vettore n_cells dimensionale
# alpha: vettore bidimensionale (2 stati)
# beta: scalare
# gamma: scalare
s_withTau <- function(tau, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma){
  
  if(is.na(u0_on)){
    u0_on = alpha[1]/beta
  }
  if(is.na(s0_on)){
    s0_on = alpha[1]/gamma
  }
  n_cells <- length(k)
  
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
  exp_gammaOFF <- exp(-gamma*tau[k == 0])
  exp_gammaON <- exp(-gamma*tau[k == 2])
  exp_betaOFF <- exp(-beta*tau[k == 0])
  exp_betaON <- exp(-beta*tau[k == 2])
  
  u0 <- ifelse(k == 2, u0_on,  u0_off)
  s0 <- ifelse(k == 2, s0_on,  s0_off)
  res <- rep(NA, n_cells)
  
  res[k == 0] <- s0[k == 0]*exp_gammaOFF + alpha[1]/gamma*(1-exp_gammaOFF)+ (alpha[1]-beta*u0[k == 0])/(gamma-beta)*(exp_gammaOFF - exp_betaOFF)
  res[k == 2] <- s0[k == 2]*exp_gammaON  + alpha[2]/gamma*(1-exp_gammaON) + (alpha[2]-beta*u0[k == 2])/(gamma-beta)*(exp_gammaON  - exp_betaON)
  return(res)
}


# calcola s per un gene e per una cellula fissati, per tutte le iterazioni dell'MCMC
# passiamo, per ogni parametro una catena con il valore di questo parametro nelle varie iterazioni
s_withTauMCMC <- function(tau, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma){
  
  if(any(is.na(u0_on))){
    u0_on = alpha[1,]/beta
  }
  if(any(is.na(s0_on))){
    s0_on = alpha[1,]/gamma
  }
  sampleToSave <- length(tau)
  
  # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
  exp_gamma <- exp(-gamma*tau)
  exp_beta <- exp(-beta*tau)
  
  u0 <- ifelse(k == 2, u0_on,  u0_off)
  s0 <- ifelse(k == 2, s0_on,  s0_off)
  res <- rep(NA, sampleToSave)
  
  res <- s0*exp_gamma + (alpha[1,]*c(k == 0) + alpha[2,]*c(k==2))/gamma*(1-exp_gamma)+ (alpha[1,]*c(k == 0) + alpha[2,]*c(k == 2)-beta*u0)/(gamma-beta)*(exp_gamma - exp_beta)
  
  return(res)
}

# s_withTauMCMC_new <- function(tau, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma, subtypeCell, typeCellT0_off){
#   u0_off_new <- array(NA, dim = dim(tau))
#   s0_off_new <- array(NA, dim = dim(tau))

#   for(sty in unique(subtypeCell)){
#     tyT0_off <- typeCellT0_off[which(subtypeCell == sty)[1]]
#     u0_off_new[sty, , ] = u0_off[tyT0_off, , ]
#     s0_off_new[sty, , ] = s0_off[tyT0_off, , ]
#   }

  
#   sSS_on <- alpha[,2,]/gamma
#   sSS_on <- array(rep(as.vector(sSS_on),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   sSS_on <- aperm(sSS_on, c(3, 1, 2))

#   s0_on <- alpha[, 1, ]/gamma
#   s0_on <- array(rep(as.vector(s0_on),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   s0_on <- aperm(s0_on, c(3, 1, 2))

#   u0_on <- alpha[, 1, ]/beta
#   u0_on <- array(rep(as.vector(u0_on),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   u0_on <- aperm(u0_on, c(3, 1, 2))
 

#   u0 <- ifelse(k == 2, u0_on, u0_off_new)
#   s0 <- ifelse(k == 2, s0_on, s0_off_new)
#   alpha <- ifelse(k == 2, alpha[, 2,], alpha[, 1,])

#   beta <- array(rep(as.vector(beta),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   beta <- aperm(beta, c(3, 1, 2))
#   gamma <- array(rep(as.vector(gamma),dim(tau)[1]), dim = c(dim(tau)[2], dim(tau)[3], dim(tau)[1]))
#   gamma <- aperm(gamma, c(3, 1, 2))
#   expBeta <- exp(-beta*tau)
#   expGamma <- exp(-gamma*tau)

  
#   res <- s0*expGamma + (s0_on*(k == 0) + sSS_on*(k ==2))*(1-expGamma)+ (alpha-beta*u0)/(gamma-beta)*(expGamma - expBeta)
  
#   return(res)
# }



u_and_s_withTauMCMC_new <- function(tau, u0_off, u0_on = NA, s0_off, s0_on = NA, k, alpha, beta, gamma, subtypeCell, typeCellT0_off){    
    u0_on <- alpha[, 1, ]/beta
    s0_on <- alpha[, 1, ]/gamma
        
    u0_on = sweep(k == 2, MARGIN = c(2, 3), STATS = u0_on, FUN = "*")
    s0_on = sweep(k == 2, MARGIN = c(2, 3), STATS = s0_on, FUN = "*")

    u0_off_new <- array(NA, dim = dim(tau))
    s0_off_new <- array(NA, dim = dim(tau))

    for(sty in unique(subtypeCell)){
        tyT0_off <- typeCellT0_off[which(subtypeCell == sty)[1]]
        u0_off_new[sty, , ] = u0_off[tyT0_off, , ]
        s0_off_new[sty, , ] = s0_off[tyT0_off, , ]
    }

    u0 = u0_on * (k == 2) + u0_off_new * (k == 0)
    s0 = s0_on * (k == 2) + s0_off_new * (k == 0)
    
    uSS_on <- alpha[,2,]/beta
    
    expBeta <- exp(-sweep(tau, MARGIN = c(2, 3), STATS = beta, FUN = "*"))
    p2 <- sweep(k == 0, MARGIN = c(2, 3), STATS = alpha[, 1, ]/beta, FUN = "*")*(1-expBeta)

    expGamma <- exp(-sweep(tau, MARGIN = c(2, 3), STATS = gamma, FUN = "*"))
    alpha = sweep(k == 2, MARGIN = c(2, 3), STATS = alpha[, 2, ], FUN = "*") + sweep(k == 0, MARGIN = c(2, 3), STATS = alpha[, 1, ], FUN = "*")

    resS <- s0*expGamma + sweep(alpha, MARGIN = c(2, 3), STATS = gamma, FUN = "/")*(1-expGamma) + sweep((alpha-sweep(u0, MARGIN = c(2, 3), STATS = beta, FUN = "*")), MARGIN = c(2, 3), STATS = (gamma-beta), FUN = "/")*(expGamma - expBeta)
    resU <- u0*expBeta + p2 + sweep(k==2, MARGIN = c(2, 3), STATS = uSS_on, FUN = "*")*(1-expBeta)

  return(list(u = resU, s = resS))
}




s_di_uMCMC <- function(u, u0_off, s0_off, k, alpha, beta, gamma){
  if(any(beta == gamma)){
    print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }
  
  u0_on = alpha[1,]/beta
  s0_on = alpha[1,]/gamma

  sampleToSave <- length(u)
  
  u0 <- ifelse(k == 2, u0_on,  u0_off)
  s0 <- ifelse(k == 2, s0_on,  s0_off)
  res <- rep(NA, sampleToSave)
   
  
  res <- (s0- (alpha[1,]*c(k == 0) + alpha[2,]*c(k == 2))/gamma + (alpha[1,]*c(k == 0) + alpha[2,]*c(k == 2) - beta*u0)/(gamma-beta))*((beta*u-alpha[1,]*c(k == 0) - alpha[2,]*c(k == 2))/(beta*u0- alpha[1,]*c(k == 0) - alpha[2,]*c(k == 2)))^(gamma/beta) -  ((alpha[1,]*c(k == 0) + alpha[2,]*c(k == 2))*beta)/(gamma*(gamma-beta)) + beta/(gamma-beta)*u
  
  return(res)
}






u_withTau_matrice <- function(tau, u0_off, u0_on = NA, k, alpha, beta, typeCell, typeCellT0_off){
  print("CONTROLLARE")
  if(is.null(dim(k))){
    u <- u_withTau(tau, u0_off, u0_on = NA, k, alpha, beta)
  }else if(length(dim(k)) == 2){ # n_cells x n_genes
    if(sum(which(is.na(u0_on)))>0){
      u0_on <- alpha[,1]/beta
      u0_on <- matrix(rep(u0_on, max(typeCell)), nrow = max(typeCell), ncol = n_genes, byrow = TRUE)
    }
    
    betaMat <- matrix(rep(beta, max(typeCell)), nrow = max(typeCell), ncol = n_genes, byrow = TRUE)
    expBetaOFF <- exp(-tau*betaMat)
    expBetaON <- expBetaOFF
    expBetaOFF[k == 0] <- 0
    expBetaON[k == 2] <- 0
    res <- matrix(NA, nrow = dim(tau)[1], ncol = dim(tau)[2])
    u0 <- res

    res1 <- expBetaOFF*u0_off + (1-expBetaOFF)*dim(matrix(rep(alpha[,2]/beta, max(typeCell)), nrow = max(typeCell), ncol = n_genes, byrow = TRUE))
    res1[k == 2] <- 0

    res2 <- expBetaON*u0_on  + (1-expBetaON)%*%diag(alpha[,2]/beta)
    res2[k == 0] <- 0
    res <- res1 + res2
    return(res)
  }else if(length(dim(k)) == 3){ # n_cells x n_genes x n_iter
     if(sum(which(is.na(u0_on)))>0){
      u0_on <- alpha[,1,]/beta
      u0_on <- array(rep(u0_on, max(typeCell)), dim(max(typeCell), n_genes, dim(tau)[3]), byrow = TRUE)
    }
    
    betaMat <- matrix(rep(beta, max(typeCell)), dim(max(typeCell), n_genes, dim(tau)[3]), byrow = TRUE)
    expBetaOFF <- exp(-tau*betaMat)
    expBetaON <- expBetaOFF
    expBetaOFF[k == 0] <- 0
    expBetaON[k == 2] <- 0
    res <- matrix(NA, nrow = dim(tau)[1], ncol = dim(tau)[2])
    u0 <- res

    res1 <- expBetaOFF*u0_off + (1-expBetaOFF)*dim(matrix(rep(alpha[,2]/beta, max(typeCell)), nrow = max(typeCell), ncol = n_genes, byrow = TRUE))
    res1[k == 2] <- 0

    res2 <- expBetaON*u0_on  + (1-expBetaON)%*%diag(alpha[,2]/beta)
    res2[k == 0] <- 0
    res <- res1 + res2
  }
}


s_matrice <- function(tau, k, alpha1, alpha2, beta, gamma){
  print("La formula va riiplementata usando u0, s0 e t in modo diverso")
  # if(is.null(dim(k))){
  #   n_typeC <- 1
  #   n_genes <- length(k)
  # }else{
  #   n_typeC <- dim(k)[2]
  #   n_genes <- dim(k)[1]
  # }
  
  # tau0 <- ifelse(k==0, tau, 0)
  # tau2 <- ifelse(k==2, tau, 0)
  
  # exp_betaOFF <- exp(-beta*tau0)
  # exp_betaON <- exp(-beta*tau2)
  # exp_gammaOFF <- exp(-gamma*tau0)
  # exp_gammaON <- exp(-gamma*tau2)
  
  # # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
  # res <- rep(NA, n_typeC*n_genes)
  # attr(res, "dim") <- c(n_genes, n_typeC)
  
  # res <- ifelse(k == 0, alpha2/gamma*exp_gammaOFF + alpha1/gamma*(1-exp_gammaOFF)+ (alpha1-beta*alpha2/beta)/(gamma-beta)*(exp_gammaOFF - exp_betaOFF),
  #           alpha1/gamma*exp_gammaON + alpha2/gamma*(1-exp_gammaON)+ (alpha2-beta*alpha1/beta)/(gamma-beta)*(exp_gammaON - exp_betaON))
  
  # return(res)
}


#------------------------------------------
#    S DI U 
#------------------------------------------
# s di u  (for beta /not = gamma)

# s_di_u <- function(u, u0_off, s0_off, k, alpha, beta, gamma){
#   if(beta == gamma){
#     print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
#   }
#   n_cells <- length(k)
  
#   u0 <- rep(NA, n_cells)              # n_cells dimensionale perche per ogni cellula possiamo partire o dallo steady state on o da quello off, a seconda del suo stato 
#   # ma potre avere solo due possibili valori 
#   s0 <- rep(NA, n_cells)
#   res <- rep(NA, n_cells)
#   a <- rep(NA, n_cells)
  
#   repr <- which(k == 0)
#   induc <- which(k == 2)
  
#   # fase repressiva/ss_off
#   u0[repr] <- u0_off
#   s0[repr] <- s0_off
#   a[repr]  <- alpha[1]
  
#   # fase induttiva/ss_on
#   u0[induc] <- alpha[1]/beta
#   s0[induc] <- alpha[1]/gamma
#   a[induc]  <- alpha[2] 
  
#   res <- (s0- a/gamma + (a-beta*u0)/(gamma-beta))*((beta*u-a)/(beta*u0- a))^(gamma/beta) -  (a*beta)/(gamma*(gamma-beta)) + beta/(gamma-beta)*u
  
#   return(res)
# }


s_di_u <- function(u, u0_off, s0_off, u0_on, s0_on, k, alpha, beta, gamma){
  if(beta == gamma){
    print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }

  if(k == 0){
    a <- alpha[1]
    u0 <- u0_off
    s0 <- s0_off
  }else{
    a <- alpha[2]
    u0 <- u0_on
    s0 <- s0_on
  }
  res <- (s0 - a/gamma + (a-beta*u0)/(gamma-beta))*((beta*u-a)/(beta*u0 - a))^(gamma/beta) -  (a*beta)/(gamma*(gamma-beta)) + beta/(gamma-beta)*u
  
  return(res)
}




# s_tilde ON di u_tilde (for beta /not = gamma), where
#           - s_tilde /in [0,1] � la posizione di s nell'intervallo riscalato
#           - u_tilde /in [0,1] � la posizione di u nell'intervallo riscalato
# restituisce sempre il valore di s_tilde t.c. (u_tilde, s_tilde) � sul ramo ON
s_tilde_di_u_tilde <- function(u_tilde, alpha, beta, gamma){
  if(beta == gamma){
    print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }
  print("La funzione deve essere rimplementata con u0, s0, t0")
  # n_cells <- length(u_tilde)
  # res <- rep(NA, n_cells)
  # res <- beta/(gamma - beta)*((1-u_tilde)^(gamma/beta) - (1-u_tilde*gamma/beta))
  # return(res)
}


# -----------------------------
#     U DI S
# -----------------------------
u_di_s <- function(s, k, alpha, beta, gamma){
  if(beta == gamma){
    print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }
  print("La funzione deve essere rimplementata con u0, s0, t0")
  # n_cells <- length(k)
  # u <- rep(NA, n_cells)
  # for(c in 1:n_cells){
  #   if(k[c] != 1 & k[c]!= 3){
  #     r = uniroot((function (y) s_di_u(y,k = k[c], alpha = alpha, beta = beta, gamma = gamma ) - s[c]), lower =(alpha[1])/beta, upper = (alpha[2])/beta )[1]
  #     u[c] = r[[1]]
  #   }
  # }
  # return(u)
}

# u_tilde ON di s_tilde (for beta /not = gamma), where
#           - s_tilde /in [0,1] � la posizione di s nell'intervallo riscalato
#           - u_tilde /in [0,1] � la posizione di u nell'intervallo riscalato
# restituisce sempre il valore di u_tilde t.c. (u_tilde, s_tilde) � sul ramo ON
u_tilde_di_s_tilde <- function(s_tilde, alpha, beta, gamma){
  if(beta == gamma){
    print("WARNING: this formula holds for beta not equal gamma, so in this case you can not use it!")
  }
  print("La funzione deve essere rimplementata con u0, s0, t0")
  # n_cells <- length(s_tilde)
  # u_tilde <- rep(NA, n_cells)
  # for(c in 1:n_cells){
  #   r = uniroot((function (y) s_tilde_di_u_tilde(y,alpha = alpha, beta = beta, gamma = gamma) - s_tilde[c]), lower =0, upper = 1)
  #   u_tilde[c] = r[[1]]
  # }
  
  # return(u_tilde)
}

# likelihood function: prodotto delle likelihood delle osservazioni delle singole cellule
#                      per ogni cellula abbiamo il prodotto di due binomiali negative (indipendenti)
# Y_u e Y_s sono vettori n_cells dimensionali con le osservazioni per lo spliced
# e per l'unspliced per quel gene (sono colonne delle matrici che avevamo in Python)
# alpha vettore bidimensionale
# beta, gamma scalari
# pos_u vettore n_cells dimensionale
# pos_s vettore n_cells dimensionale
# k vettore n_cells dimensionale
# eta scalare



likelihood <- function(Y_u, Y_s, alpha, beta, gamma, t, t0_off, t0_on, u0_off, u0_on, s0_off, s0_on, k, eta, catt, log = TRUE){
  n_cells <- length(k)
  
  pos_u <- u(t, t0_off, t0_on, u0_off, u0_on, k, alpha, beta)
  pos_s <- s(t, t0_off, t0_on, u0_off, u0_on, s0_off, s0_on, k, alpha, beta, gamma)
  # media delle due binomiali negative
  mu_u <- catt*pos_u          
  mu_s <- catt*pos_s
  
  loglik <- rep(NA, n_cells)
   
  loglik <- dnbinom(Y_u, size= 1/eta, mu = mu_u, log = TRUE) + dnbinom(Y_s, size= 1/eta, mu = mu_s, log = TRUE)
  loglik <- sum(loglik)
  
  if(log == FALSE){
    loglik <- exp(loglik)
  }
  return(loglik)
}

## It performs the accept-reject procedure of the MCMC
# input:
# - logNum: the logarithm of the numerator of the acceptance ratio
# - logDen: the logarithm of the denominator of the acceptance ratio
# - prop: the new proposed values for the parameters
# - prev: the previous values of the parameters
acceptMH <- function(logNum, logDen, prop, prev){
  Alpha <- min(1, exp(logNum - logDen))
  # vediamo se accettare
  acc <- runif(1,0.0,1.0)
  moved <- FALSE
  if(acc<Alpha){
    # accetto
    res <- prop
    moved <- TRUE
  }else{
    # mantengo quelli precedenti
    res <- prev
  }
  return(list("Alpha"= Alpha, "update"= res, "moved" = moved))
}


# X ~ gamma(shape, scale)
# W:=log(X)--> F_W(w) = P(W<= w)= P(log(X)<=w) = P(X <= e^w) = F_X(e^w)
#          --> f_W(w) = f_X(e^w)*e^w
dgamma_log<- function(w, shape, scale, log = TRUE){
  # Abbiamo tolto direttamente la costante di normalizzazione
  res <- w*shape-exp(w)/scale
  # res <- dgamma(exp(w), shape = shape, scale = scale, log = FALSE)*exp(w)
  # se exp(w) = Inf la densit� calcolata verrebbe 0*inf e quindi ritorna NaN--> Facendo i conti viene 0
  # if(exp(w) == Inf){
  #   res <- 0
  # }
  # # se exp(w) = 0 d� NaN-> facendo i conti viene 0
  # if(exp(w) == 0){
  #   res <- 0
  # }
  
  # if log = TRUE, we directl return the value we have just computed, 
  # that is already the logarithm 
  # otherwise we compute its exponential
  if(!log){
    res <- exp(res) 
  }
  return(res)
}

 

#------------------------
# !!!ANCORA DA CONTROLLARE!!!
# G FISSATO
# beta, gamma sono scalari (nell'MCMC saranno i parametri del passo precedente PER QUEL GENE)
# alpha � un vettore bidimensionale con alpha_off e alpha_on
# Y_u e Y_s sono vettori n_cells - dimensionali, con l'espressione per il gene fissato

approxTau <- function(Y_u, Y_s, alpha, beta, gamma, k){
  # n_cells <- length(Y_u)
  # tau <- rep(NA, n_cells)
  

  # beta_tilde <- beta/(gamma-beta)
  # sinf <- rep(NA, n_cells)
  # uinf <- rep(NA, n_cells)
  # s0 <- rep(NA, n_cells)
  # u0 <- rep(NA, n_cells)
  
  # sinf[which(k == 0 | k == 1 )] <- alpha[1]/gamma
  # sinf[which(k == 2 | k == 3 )] <- alpha[2]/gamma
  # uinf[which(k == 0 | k == 1)] <- alpha[1]/beta 
  # uinf[which(k == 2 | k == 3)] <- alpha[2]/beta 
  
  # s0[which(k == 0 | k == 1)] <- alpha[2]/gamma
  # s0[which(k == 2 | k == 3)] <- alpha[1]/gamma
  # u0[which(k == 0 | k == 1)] <- alpha[2]/beta
  # u0[which(k == 2 | k == 3)] <- alpha[1]/beta
 
  # s_tilde <- Y_s - beta_tilde*Y_u     # vettore n_cells - dimensionale
  # sinf_tilde <- sinf - beta_tilde*uinf 
  # s0_tilde <- s0 - beta_tilde*u0      
  
  
  # if(beta>gamma){
  #   tau <- (s_tilde - sinf_tilde)/(s0_tilde - sinf_tilde)
  #   # alcuNi Tau cos� sono negativi...
  #   # loro li limitano tra 0 e 1 per fare il logaritmo
  #   tau <- -1/gamma*sapply(tau, clipped_log) 
    
  # }else{
  #   tau <- (Y_u - uinf)/(u0 - uinf)
  #   tau <-  -1/beta*sapply(tau,clipped_log)
  # }
  
  # # tau[which(k == 1 | k == 3)] <- +Inf 
  # return(tau)
  print("Dunzione da implementare usando anche u0, s0, t0")
}

clipped_log <- function(x, lb = 0, up = 1, epsilon = 1e-6){
  if(x <= lb){
    x <- lb +epsilon
  }else if(x>= up){
    x <- up - epsilon
  }
  return(log(x))
}


tau_inv <- function(u, u0_off, u0_on, k, alpha, beta, gamma){
  tau <- rep(NA, length(k))
  uinf <- rep(NA, length(k))
  u0 <- rep(NA, length(k))
  
  
  uinf[which(k == 0)] <- alpha[1]/beta
  u0[which(k == 0)] <- u0_off
  uinf[which(k == 2)] <- alpha[2]/beta
  u0[which(k == 2)] <- u0_on
  
  tau <- (u - uinf)/(u0 - uinf)
  # tau <-  -1/beta*sapply(tau,clipped_log)
  tau <-  -1/beta*log(tau)
  return(tau)
}

# If y ~ truncatedNormal, then log(y) has this distribution. Here we use only the kernel of the density 
dTrunc_norm_log <- function(x, mean, sd, a = 0, b = +Inf, log = TRUE){
  # res <- log(dtruncnorm(exp(x), a = a, b = b, mean = mean , sd = sd)) + x
  

  # logaritmo del kernel di una normale troncata tra 0 e +Inf
  xi <- (exp(x)-mean)/sd
  res <- -1/2*xi^2 + x

  if(log == FALSE){
    res <- exp(res)
  }
  return(res)
}

# if    W := logit(X) is distributed as a normal with mean and standard deviation given 
# then  U has this distribution
d_invLogit <- function(x, mean, sd, log = TRUE){
  d <- dnorm(logit(x), mean = mean, sd = sd, log = log)
  if(log){
    res <- d -log(x*(1-x))
  }else{
    res <- d/(x*(1-x))
  }
  return(res)
}

u0 <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, k, alpha, beta){
  if(k == 0){
    if(any(is.na(u0_on))){
      u0_on = alpha[1]/beta
    }
    if(any(is.na(t0_off))){
      t <- rep(NA, length(k))
      t[which(is.na(t0_off))] = rep(Inf, length(which(is.na(t0_off)))) # raggiungiamo effettivamente lo osteady state alto
    }else{
      t = t0_off # c'è lo switch prima di raggiungere lo SS alto
    }
    res <- u(t, t0_off, t0_on, u0_off = NA, u0_on = u0_on, k = rep(2, length(t)), alpha, beta)
  }
  if(k == 2){
    if(any(is.na(u0_off))){
      u0_off = alpha[2]/beta
    }
    res <- u(Inf, t0_off = 0, t0_on = NA, u0_off = u0_off, u0_on = NA, k = rep(0, length(t)), alpha, beta)
  }
  return(res)
}



s0 <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k, alpha, beta, gamma){
  if(k == 0){
    if(any(is.na(s0_on))){
      s0_on = alpha[1]/gamma
    }
    if(any(is.na(u0_on))){
      u0_on = alpha[1]/beta
    }
    if(any(is.na(t0_off))){
      t <- rep(NA, length(k))
      t[which(is.na(t0_off))] = rep(Inf, length(which(is.na(t0_off)))) # raggiungiamo effettivamente lo osteady state alto
    }else{
      t = t0_off # c'è lo switch prima di raggiungere lo SS alto
    }
    if(any(is.na(u0_off))){
      u0_off <- u0(t0_off = t0_off, t0_on = t0_on, u0_off = NA, u0_on = u0_on, k = 0, alpha = alpha, beta = beta)
    }
   
    res <- s(t = t, t0_off = t0_off, t0_on = t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(t)), alpha, beta, gamma)
    
  }
  if(k == 2){
    if(any(is.na(s0_off))){
      s0_off = alpha[2]/gamma
    }
    if(any(is.na(u0_off))){
      u0_off = alpha[2]/beta
    }
    res <- s(rep(Inf, length(t0_on)), t0_off = 0, t0_on = NA, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(t0_on)), alpha = alpha, beta = beta, gamma = gamma)
    
  }
  
  return(res)
}



u0_MCMC_g <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, k, alpha, beta, g, tyT0_off){
    if(length(k) > 1){
        k <- unique(k)
    }

    if(k == 0){
        u0_on = alpha[g, 1, ]/beta[g,]
        t = t0_off[tyT0_off, g,] # c'è lo switch prima di raggiungere lo SS alto

        if(any(is.na(t0_off[tyT0_off, g, ]))){
            t <- rep(NA, length(u0_on))
            t[which(is.na(t0_off[tyT0_off, g,]))] = rep(Inf, length(which(is.na(t0_off[tyT0_off, g,])))) # raggiungiamo effettivamente lo osteady state alto
        }

        
        t0_on <- rep(0, length(u0_on))

        tau <- t - t0_on
        expBetaOFF <- exp(-beta[g,]*tau)

        res <- u0_on*expBetaOFF + alpha[g, 2,]/beta[g,]*(1-expBetaOFF)
    }
    if(any(k == 2)){
        res <-  alpha[g,1,]/beta[g,]
    }
    return(res)
}


u0_MCMC <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, k, alpha, beta){
    if(length(k) > 1){
        k <- unique(k)
    }

    if(k == 0){
        u0_on = alpha[, 1, ]/beta
        u0_on <- array(rep(as.vector(u0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
        u0_on <- aperm(u0_on, c(3, 1, 2))


        t <- t0_off # c'è lo switch prima di raggiungere lo SS alto

        # if(any(is.na(t0_off))){
        #     t <- matrix(NA, dim(u0_on))
        #     t[which(is.na(t0_off[tyT0_off, ,]))] = rep(Inf, length(which(is.na(t0_off[tyT0_off, ,])))) # raggiungiamo effettivamente lo steady state alto
        # }
        
        t0_on <- array(0, dim = dim(u0_on))

        tau <- t - t0_on

        uSS_ON <- alpha[, 2,]/beta
        uSS_ON <- array(rep(as.vector(uSS_ON), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
        uSS_ON <- aperm(uSS_ON, c(3, 1, 2))


        beta <- array(rep(as.vector(beta), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
        beta <- aperm(beta, c(3, 1, 2))

        expBetaOFF <- exp(-beta*tau)

        res <- u0_on*expBetaOFF + uSS_ON*(1-expBetaOFF)
    }


    if(any(k == 2)){
        res <-  alpha[,1,]/beta
    }
    return(res)
}



s0_MCMC_g <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k, alpha, beta, gamma, g, tyT0_off){
    if(length(k) == 1){
        k <- rep(k, length(t0_off[tyT0_off,g,]))
    }

    if(any(k == 0)){
        if(any(is.na(s0_on))){
            s0_on = alpha[g,1,]/gamma[g,]
        }
        if(any(is.na(u0_on))){
            u0_on = alpha[g,1,]/beta[g,]
        }
        t = t0_off[tyT0_off, g,] # c'è lo switch prima di raggiungere lo SS alto
        if(any(is.na(t0_off[tyT0_off, g,]))){
            t <- rep(NA, length(u0_on))
            t[which(is.na(t0_off[tyT0_off, g,]))] = rep(Inf, length(which(is.na(t0_off[tyT0_off, g,])))) # raggiungiamo effettivamente lo osteady state alto
        }

   
        t0_on <- rep(0, length(u0_on))
        tau <- t - t0_on
        # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
        exp_gammaOFF <- exp(-gamma[g,]*tau)
        exp_betaOFF <- exp(-beta[g,]*tau)
        
  
        res <- s0_on*exp_gammaOFF + alpha[g,2,]/gamma[g,]*(1-exp_gammaOFF)+ (alpha[g,2,]-beta[g,]*u0_on)/(gamma[g,]-beta[g,])*(exp_gammaOFF - exp_betaOFF)    
    }
    if(any(k == 2)){
        res <-  alpha[g,1,]/gamma[g,]
    }  
    return(res)
}


s0_MCMC <- function(t0_off = NA, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k, alpha, beta, gamma){
     if(length(k) > 1){
        k <- unique(k)
    }

    if(k == 0){
      s0_on = alpha[,1,]/gamma
      s0_on <- array(rep(as.vector(s0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
      s0_on <- aperm(s0_on, c(3, 1, 2))


      u0_on = alpha[,1,]/beta
      u0_on <- array(rep(as.vector(u0_on),dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
      u0_on <- aperm(u0_on, c(3, 1, 2))
        
      t = t0_off # c'è lo switch prima di raggiungere lo SS alto
          
      t0_on <- array(0, dim = dim(u0_on))

      tau <- t - t0_on

      sSS_ON <- alpha[,2,]/gamma
      sSS_ON <- array(rep(as.vector(sSS_ON), dim(t0_off)[1]), dim = c(dim(t0_off)[2], dim(t0_off)[3], dim(t0_off)[1]))
      sSS_ON <- aperm(sSS_ON, c(3, 1, 2))

        # alpha vettore bidimensionale: alpha[1] = alpha_0/1, alpha[2] = diff t.c. alpha_2/3 = alpha_0/1 + diff
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
    if(any(k == 2)){
        res <-  alpha[,1,]/gamma
    }  
    return(res)
}

# ------------------------
# PLOT 

# ------------------------


OLDplot_sVSu <- function(t0_off_real, t0_on_real = 0, alpha_real, beta_real, gamma_real, pos_u_real, pos_s_real, k_real, g, point = NA, Y_u, Y_s, tipoCellula, catt_real, legendType = TRUE, legendNumber = TRUE,...){
  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  u0_off <- u0(k = 0, alpha = alpha_real, beta = beta_real)
  
  u0_on <- u0(k = 2, alpha = alpha_real, beta = beta_real)
  
  s0_off <- s0(k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
  s0_on <- s0(k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
  
  
  # ramo on
  if(any(!is.na(Y_u))){
    xlimY_u_up <- max(Y_u/catt_real)
    xlimY_u_down <- min(Y_u/catt_real)
  }else{
    xlimY_u_up <- alpha_real[2]/beta_real
    xlimY_u_down <- alpha_real[1]/beta_real
  }
  if(any(!is.na(Y_s))){
    xlimY_s_up <- max(Y_s/catt_real)
    xlimY_s_down <- min(Y_s/catt_real)
  }else{
    xlimY_s_up <- alpha_real[2]/gamma_real
    xlimY_s_down <- alpha_real[1]/gamma_real
  }


  u_plot <- u(seq(t0_on_real, 5000, 0.05), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, 5000, 0.05))), alpha_real, beta_real)
  s_plot <- s(seq(t0_on_real, 5000, 0.05), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, 5000, 0.05))), alpha_real, beta_real, gamma_real)

  res <- rbind(u_plot, s_plot)
  main <- paste("Gene", g)
  
  if(length(unique(tipoCellula))== 1){
    main <- paste(main, ", typeCell", unique(tipoCellula))
  }
  plot(s_plot, u_plot, type = "l", xlab = "s", ylab = "u", main = main, lty = 2,  ylim = c(min(xlimY_u_down, alpha_real[1]/beta_real)-0.3, max(xlimY_u_up, alpha_real[2]/beta_real)+0.3), xlim = c(min(xlimY_s_down, alpha_real[1]/gamma_real)-0.3, max(xlimY_s_up, alpha_real[2]/gamma_real)+0.3), col = "blue", cex.axis = 1.7, cex.lab = 1.5, ...)

  # ramo off
  u_plot <- u(seq(0, 5000, 0.05), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha_real, beta_real)
  s_plot <- s(seq(0, 5000, 0.05), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha_real, beta_real, gamma_real)
  

  lines(s_plot, u_plot, type = "l", lty = 2, col = "red",...)
  

  if(t0_off_real != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
    u0_off <- u0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real)
    u0_on <- u0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real)
    s0_off <- s0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
    s0_on <- s0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)

  # ramo on 
    u_plot <- u(seq(t0_on_real, t0_off_real, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.01))), alpha_real, beta_real)
    s_plot <- s(seq(t0_on_real, t0_off_real, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.01))), alpha_real, beta_real, gamma_real)
    lines(s_plot, u_plot, type = "l", lty = 1, col = "blue", ...)
  # ramo off

    u_plot <- u(seq(t0_off_real, t0_off_real + 5000, 0.05), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 5000, 0.05))), alpha_real, beta_real)
    s_plot <- s(seq(t0_off_real, t0_off_real + 5000, 0.05), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 5000, 0.05))), alpha_real, beta_real, gamma_real)
    lines(s_plot, u_plot, type = "l", lty = 1, col = "red", ...)
  }else{
    lines(res[2,], res[1,], type = "l", lty = 1, col = "blue", ...)
    lines(s_plot, u_plot, type = "l", lty = 1, col = "red", ...)
  }
  
 
 res <- list("upper" = res, "lower" = rbind(u_plot, s_plot))

typeC <- unique(tipoCellula)
colors <- rainbow(length(typeC))
  for(i in 1:length(typeC)){
      x <- typeC[i]
      # points(pos_s_real[tipoCellula == x],  pos_u_real[tipoCellula == x], col = tipoCellula[tipoCellula == x], pch = 16, ...)
      points(pos_s_real[tipoCellula == x],  pos_u_real[tipoCellula == x], col = colors[i], pch = i, ...)
    }
  
  if(any(!is.na(Y_u)) & any(!is.na(Y_s))){
    for(i in 1:length(typeC)){
          x <- typeC[i]
          points(Y_s[tipoCellula == x]/catt_real[tipoCellula == x], Y_u[tipoCellula == x]/catt_real[tipoCellula == x], col = tipoCellula[tipoCellula == x], pch = 8,...)
    }
  }
  tab <- data.frame(table(tipoCellula))
  
  if(is.na(point[1])){
    if(length(unique(tipoCellula))>1){
      legend <- paste(tab$tipoCellula, tab$Freq)
    }else{
      legend <- paste(unique(tipoCellula), length(tipoCellula))
    }
    if(legendType) legend("topright", legend = legend, col = typeC, pch = 16, cex = 0.8 ) 
  }
  totOn <- 0
  totOff <- 0
  if(legendNumber){
    if(length(k_real)>1){
        for(i in 1:n_typeC){
          if(k_real[i] == 0){
            totOff <- totOff + length(which(tipoCellula == typeC[i]))
          }else{
            totOn <- totOn + length(which(tipoCellula == typeC[i]))
          }
        } 
      }else{
        if(k_real == 0){
          totOff <- totOff + length(tipoCellula)
        }else{
          totOn <- totOn + length(tipoCellula)
        }
      }
      legend("bottomright", legend = c(paste("Inductive",totOn), paste("Repressive", totOff)), cex = 0.8)
  }
  return(res)
}




OLDplot_sVSu_scVelo <- function(t0_off_real, t0_on_real = 0, alpha_real, beta_real, gamma_real, scaling, u0_offset, s0_offset, pos_u_real, pos_s_real, k_real, g, point = NA, Y_u, Y_s, tipoCellula, catt_real, legendType = TRUE, legendNumber = TRUE, add = FALSE, col, ...){

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  u0_off <- u0(k = 0, alpha = alpha_real, beta = beta_real)
  
  u0_on <- u0(k = 2, alpha = alpha_real, beta = beta_real)
  
  s0_off <- s0(k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
  s0_on <- s0(k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
 
  # ramo on

  u_plot <- u(seq(t0_on_real, 5000, 0.05), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, 5000, 0.05))), alpha_real, beta_real)
  s_plot <- s(seq(t0_on_real, 5000, 0.05), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, 5000, 0.05))), alpha_real, beta_real, gamma_real)

  u_plot= u_plot * scaling + u0_offset
  s_plot = s_plot + s0_offset
  
  main <- paste("Gene", g)
  
  if(length(unique(tipoCellula))== 1){
    main <- paste(main, ", typeCell", unique(tipoCellula))
  }

  if(!add){
      plot(s_plot, u_plot, type = "l", xlab = "s", ylab = "u", main = main, lty = 2,  ylim = c(0,60), xlim = c(0,60), cex.axis = 1.7, cex.lab = 1.5, col = col,...)
  }else{
    lines(s_plot, u_plot, type = "l", xlab = "s", ylab = "u", main = paste("Gene", g), lty = 2, col = col,...)
  }
  res <- rbind(u_plot, s_plot)

  # ramo off
  u_plot <- u(seq(0, 5000, 0.05), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha_real, beta_real)
  s_plot <- s(seq(0, 5000, 0.05), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha_real, beta_real, gamma_real)
  
  u_plot= u_plot * scaling + u0_offset
  s_plot = s_plot + s0_offset
  lines(s_plot, u_plot, type = "l", lty = 2, col = col,...)
  
  if(t0_off_real != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
    u0_off <- u0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real)
    u0_on <- u0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real)
    s0_off <- s0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
    s0_on <- s0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)

  # ramo on 
    u_plot <- u(seq(t0_on_real, t0_off_real, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.01))), alpha_real, beta_real)
    s_plot <- s(seq(t0_on_real, t0_off_real, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.01))), alpha_real, beta_real, gamma_real)

    u_plot= u_plot * scaling + u0_offset
    s_plot = s_plot + s0_offset

    lines(s_plot, u_plot, type = "l", lty = 1, col = col,...)
  # ramo off

    u_plot <- u(seq(t0_off_real, t0_off_real + 5000, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 5000, 0.01))), alpha_real, beta_real)
    s_plot <- s(seq(t0_off_real, t0_off_real + 5000, 0.01), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 5000, 0.01))), alpha_real, beta_real, gamma_real)

    u_plot= u_plot * scaling + u0_offset
    s_plot = s_plot + s0_offset

    lines(s_plot, u_plot, type = "l", lty = 1,col = col,...)
  }else{
    lines(res[2,], res[1,], type = "l", lty = 1, col = col, ...)
    lines(s_plot, u_plot, type = "l", lty = 1, col = col, ...)
  }

####################################################################################
####################################################################################
####################################################################################
####################################################################################

typeC <- unique(tipoCellula)
  for(i in 1:length(typeC)){
      x <- typeC[i]
      points(pos_s_real[tipoCellula == x],  pos_u_real[tipoCellula == x], col = tipoCellula[tipoCellula == x], ...)
    }
  
}





plot_sVSu <- function(t0_off_real, t0_on_real = 0, alpha_real, beta_real, gamma_real, pos_u_real, pos_s_real, g, tipoCellula, catt_real, add = FALSE, colCell = NA, xlim = 15, ylim = 15, axisTitle.size = 5, axisText.size, title.size = 20, colReal = "red", lineSize = 1, shapePoint = 21, sizePoint = 2, gg = NA, ...){

  ####################################################################################
  ####################################################################################
  ####################################################################################
  ####################################################################################
  ####################################################################################

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  u0_off <- u0(k = 0, alpha = alpha_real, beta = beta_real)
  u0_on <- u0(k = 2, alpha = alpha_real, beta = beta_real)  
  s0_off <- s0(k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
  s0_on <- s0(k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
 
  # ramo on
  u_plot <- u(seq(t0_on_real, 500, 0.1), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, 500, 0.1))), alpha_real, beta_real)
  s_plot <- s(seq(t0_on_real, 500, 0.1), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, 500, 0.1))), alpha_real, beta_real, gamma_real)
  
  main <- paste("Gene", g)
  
  if(length(unique(tipoCellula))== 1){
    main <- paste(main, ", typeCell", unique(tipoCellula))
  }

  df <- data.frame(s_plot = s_plot, u_plot = u_plot)

  if(!add){
    gg <- ggplot(data = df, aes(x = s_plot, y = u_plot)) + geom_path(lty = 2, col = colReal, linewidth = lineSize) + xlab("s") + ylab("u") + labs(title = main) + theme(plot.title = element_text(hjust = 0.5, size = title.size), axis.title.x = element_text(size = axisTitle.size), axis.title.y = element_text(size = axisTitle.size), axis.text.x = element_text(size = axisText.size), axis.text.y = element_text(size = axisText.size)) 

    if(is.numeric(xlim)){
      gg <- gg + xlim(0, xlim) 
    }
    if(is.numeric(ylim)){
      gg <- gg + ylim(0, ylim) 
    }

  }else{
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 2, col = colReal, linewidth = lineSize)
  }



  res <- df

  # ramo off
  u_plot <- u(seq(0, 500, 0.1), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(0, 500, 0.1))), alpha_real, beta_real)
  s_plot <- s(seq(0, 500, 0.1), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(0, 500, 0.1))), alpha_real, beta_real, gamma_real)
  
  df <- data.frame(s_plot = s_plot, u_plot = u_plot)
  gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 2, col = colReal, linewidth = lineSize) 

  
  if(t0_off_real != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
    u0_off <- u0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real)
    u0_on <- u0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real)
    s0_off <- s0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
    s0_on <- s0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)

    # ramo on 
    u_plot <- u(seq(t0_on_real, t0_off_real, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.1))), alpha_real, beta_real)
    s_plot <- s(seq(t0_on_real, t0_off_real, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.1))), alpha_real, beta_real, gamma_real)

    df <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 
   
    # ramo off
    u_plot <- u(seq(t0_off_real, t0_off_real + 500, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 500, 0.1))), alpha_real, beta_real)
    s_plot <- s(seq(t0_off_real, t0_off_real + 500, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 500, 0.1))), alpha_real, beta_real, gamma_real)

    df <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 
  }else{
    gg <- gg + geom_path(data = res, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 
  }

  ####################################################################################
  ####################################################################################
  ####################################################################################
  ####################################################################################
  

  typeC <- unique(tipoCellula)

  if(sum(is.na(pos_s_real)) == 0){
    if(length(pos_s_real) == length(typeC)){
      df <- data.frame(s = pos_s_real, u = pos_u_real)
    }else{
      pos_s <- c()
      pos_u <- c()
      for(i in 1:length(typeC)){
          x <- typeC[i]
          pos_s <- c(pos_s, pos_s_real[which(tipoCellula == x)[1]])
          pos_u <- c(pos_u, pos_u_real[which(tipoCellula == x)[1]]) 
      }
      df <- data.frame(s = pos_s, u = pos_u)
    }
    if(sum(is.na(colCell)) == 1){
      if(length(typeC) >= 3){
        colCell <- brewer.pal(length(typeC), "Spectral")
      }else if(length(typeC) == 2){
        colCell <- c("darkgreen", "blue")
      }else if(length(typeC) == 1){
        colCell <- c("darkgreen")
      }
    }
    df$colCell <- colCell
    gg <- gg + geom_point(data = df, aes(x = s, y = u, fill = colCell), fill = colCell, col = colReal, size = sizePoint, shape = shapePoint)
  }

  return(gg)
}




plot_sVSu_scVelo <- function(fit_t0_off, fit_t0_on = 0, fit_alpha, fit_beta, fit_gamma, fit_scaling, fit_u0_offset, fit_s0_offset, fit_t, g, tipoCellula, add = FALSE, colCell = NA, xlim = 15, ylim = 15, axisTitle.size = 3, axisText.size = 5, title.size = 20, colReal = "red", lineSize = 1, shapePoint = 21, sizePoint = 1, gg = NA, ...){

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

# il beta che viene usato in scVelo è fit_beta*fit_scaling 
  beta <- fit_beta * fit_scaling 
  alpha <- c(0, fit_alpha)

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  u0_off <- u0(k = 0, alpha = alpha, beta = beta)
  u0_on  <- u0(k = 2, alpha = alpha, beta = beta)

  s0_off <- s0(k = 0, alpha = alpha, beta = beta, gamma = fit_gamma)
  s0_on  <- s0(k = 2, alpha = alpha, beta = beta, gamma = fit_gamma)
 
  # ramo on
  u_plot <- u(seq(fit_t0_on, 5000, 0.05), t0_off = Inf, t0_on = 0, u0_off = u0_off, u0_on = 0, k = rep(2, length(seq(fit_t0_on, 5000, 0.05))), alpha, beta)
  s_plot <- s(seq(fit_t0_on, 5000, 0.05), t0_off = Inf, t0_on = 0, u0_off = u0_off, u0_on = 0, s0_off = s0_off, s0_on = 0, k = rep(2, length(seq(fit_t0_on, 5000, 0.05))), alpha, beta, fit_gamma)

  # I RAMI VENGONO RISCALATI E TRASLATI 
  u_plot <- u_plot * fit_scaling + fit_u0_offset
  s_plot <- s_plot + fit_s0_offset
  
  main <- paste("Gene", g)
  
  if(length(unique(tipoCellula))== 1){
    main <- paste(main, ", typeCell", unique(tipoCellula))
  }

  df <- data.frame(s_plot = s_plot, u_plot = u_plot)

  if(!add){
    gg <- ggplot(data = df, aes(x = s_plot, y = u_plot)) + geom_path(lty = 2, col = colReal, linewidth = lineSize) + xlab("s") + ylab("u") + labs(title = main) + theme(plot.title = element_text(hjust = 0.5, size = title.size), axis.title.x = element_text(size = axisTitle.size), axis.title.y = element_text(size = axisTitle.size), axis.text.x = element_text(size = axisText.size), axis.text.y = element_text(size = axisText.size)) 
    
    if(is.numeric(xlim)){
      gg <- gg + xlim(0, xlim) 
    }
    if(is.numeric(ylim)){
      gg <- gg + ylim(0, ylim) 
    }
  }else{
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 2, col = colReal, linewidth = lineSize) 
  }

  
  res <- df


  # ramo off
  u_plot <- u(seq(0, 5000, 0.05), t0_off = 0, t0_on = 0, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha, beta)
  s_plot <- s(seq(0, 5000, 0.05), t0_off = 0, t0_on = 0, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(0, 5000, 0.05))), alpha, beta, fit_gamma)
  
  u_plot <- u_plot * fit_scaling + fit_u0_offset
  s_plot <- s_plot + fit_s0_offset

  df <- data.frame(s_plot = s_plot, u_plot = u_plot)
  gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 2, col = colReal, linewidth = lineSize) 
  
  if(fit_t0_off != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
    u0_off <- u0(t0_off = fit_t0_off, k = 0, alpha = alpha, beta = beta)
    u0_on <- u0(t0_on = fit_t0_on, k = 2, alpha = alpha, beta = beta)
    s0_off <- s0(t0_off = fit_t0_off, k = 0, alpha = alpha, beta = beta, gamma = fit_gamma)
    s0_on <- s0(t0_on = fit_t0_on, k = 2, alpha = alpha, beta = beta, gamma = fit_gamma)

    # ramo on 
    u_plot <- u(seq(fit_t0_on, fit_t0_off, 0.01), t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(fit_t0_on, fit_t0_off, 0.01))), alpha, beta)
    s_plot <- s(seq(fit_t0_on, fit_t0_off, 0.01), t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(fit_t0_on, fit_t0_off, 0.01))), alpha, beta, fit_gamma)

    u_plot = u_plot * fit_scaling + fit_u0_offset
    s_plot = s_plot + fit_s0_offset

    df <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 

    # ramo off
    u_plot <- u(seq(fit_t0_off, fit_t0_off + 5000, 0.01), t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(fit_t0_off, fit_t0_off + 5000, 0.01))), alpha, beta)
    s_plot <- s(seq(fit_t0_off, fit_t0_off + 5000, 0.01), t0_off = fit_t0_off, t0_on = fit_t0_on, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(fit_t0_off, fit_t0_off + 5000, 0.01))), alpha, beta, fit_gamma)

    u_plot = u_plot * fit_scaling + fit_u0_offset
    s_plot = s_plot + fit_s0_offset

    df <- data.frame(s_plot = s_plot, u_plot = u_plot)
    gg <- gg + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 
  }else{
    gg <- gg + geom_path(data = res, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) + geom_path(data = df, aes(x = s_plot, y = u_plot), lty = 1, col = colReal, linewidth = lineSize) 
  }

  ####################################################################################
  ####################################################################################
  ####################################################################################
  ####################################################################################
  # pos_s e pos_u
  n_cells <- length(fit_t)
  pos_u <- rep(NA, n_cells)
  pos_s <- rep(NA, n_cells)
  k <- rep(0, n_cells)
  k[which(fit_t < fit_t0_off)] <- 2
  pos_u <- u(fit_t, t0_off = fit_t0_off, t0_on = 0, u0_off = u0_off, u0_on = u0_on, k = k, alpha, beta)
  pos_s <- s(fit_t, t0_off = fit_t0_off, t0_on = 0, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = k, alpha, beta, fit_gamma)

  # I RAMI VENGONO RISCALATI E TRASLATI 
  pos_u <- pos_u * fit_scaling + fit_u0_offset
  pos_s <- pos_s + fit_s0_offset
  df <- data.frame(s = pos_s, u = pos_u)

  typeC <- unique(tipoCellula)
  if(sum(is.na(colCell)) == 1){
    colCell <- brewer.pal(length(typeC), "Spectral")
  }
  colori <- c()
  for(i in 1:length(typeC)){
    x <- typeC[i]
    colori <- c(colori, rep(colCell[i], length(which(tipoCellula == x))))
  }
  if(nrow(df) == length(colori)){
    df$colCell <- colori
  }else{
    df$colCell <- colori[1:nrow(df)]
  }
  gg <- gg + geom_point(data = df, aes(x = s, y = u, fill = colCell), fill = df$colCell, col = colReal, size = sizePoint, shape = shapePoint)

  return(gg)
}




plot_sVSu_MCMC <- function(t0_off_real, t0_on_real = 0, alpha_real, beta_real, gamma_real, pos_u_real, pos_s_real, k_real, g, point = NA, Y_u, Y_s, tipoCellula, catt_real, legendType = TRUE, legendNumber = TRUE, add = FALSE, colCell = NA, col, xlim = 15, ylim = 15, ...){

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  u0_off <- u0(k = 0, alpha = alpha_real, beta = beta_real)
  
  u0_on <- u0(k = 2, alpha = alpha_real, beta = beta_real)
  
  s0_off <- s0(k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
  s0_on <- s0(k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
 
  # ramo on

  u_plot <- u(seq(t0_on_real, 500, 0.1), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, 500, 0.1))), alpha_real, beta_real)
  s_plot <- s(seq(t0_on_real, 500, 0.1), t0_off = Inf, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, 500, 0.1))), alpha_real, beta_real, gamma_real)


  
  main <- paste("Gene", g)
  
  if(length(unique(tipoCellula))== 1){
    main <- paste(main, ", typeCell", unique(tipoCellula))
  }

  if(!add){
      plot(s_plot, u_plot, type = "l", xlab = "s", ylab = "u", main = main, lty = 2,  ylim = c(0,ylim), xlim = c(0,xlim), cex.axis = 1.7, cex.lab = 1.5, col = col,...)
  }else{
    lines(s_plot, u_plot, type = "l", xlab = "s", ylab = "u", main = paste("Gene", g), lty = 2, col = col,...)
  }
  res <- rbind(u_plot, s_plot)

  # ramo off
  u_plot <- u(seq(0, 500, 0.1), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(0, 500, 0.1))), alpha_real, beta_real)
  s_plot <- s(seq(0, 500, 0.1), t0_off = 0, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(0, 500, 0.1))), alpha_real, beta_real, gamma_real)
  

  lines(s_plot, u_plot, type = "l", lty = 2, col = col,...)
  
  if(t0_off_real != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
    u0_off <- u0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real)
    u0_on <- u0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real)
    s0_off <- s0(t0_off = t0_off_real, k = 0, alpha = alpha_real, beta = beta_real, gamma = gamma_real)
    s0_on <- s0(t0_on = t0_on_real, k = 2, alpha = alpha_real, beta = beta_real, gamma = gamma_real)

  # ramo on 
    u_plot <- u(seq(t0_on_real, t0_off_real, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.1))), alpha_real, beta_real)
    s_plot <- s(seq(t0_on_real, t0_off_real, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(2, length(seq(t0_on_real, t0_off_real, 0.1))), alpha_real, beta_real, gamma_real)
    lines(s_plot, u_plot, type = "l", lty = 1, col = col,...)
 
  # ramo off

    u_plot <- u(seq(t0_off_real, t0_off_real + 500, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 500, 0.1))), alpha_real, beta_real)
    s_plot <- s(seq(t0_off_real, t0_off_real + 500, 0.1), t0_off = t0_off_real, t0_on = t0_on_real, u0_off = u0_off, u0_on = u0_on, s0_off = s0_off, s0_on = s0_on, k = rep(0, length(seq(t0_off_real, t0_off_real + 500, 0.1))), alpha_real, beta_real, gamma_real)
    lines(s_plot, u_plot, type = "l", lty = 1,col = col,...)
  }else{
    lines(res[2,], res[1,], type = "l", lty = 1, col = col, ...)
    lines(s_plot, u_plot, type = "l", lty = 1, col = col, ...)
  }

####################################################################################
####################################################################################
####################################################################################
####################################################################################

  
typeC <- unique(tipoCellula)
if(sum(is.na(pos_s_real)) == 0){
  if(length(pos_s_real) <= 100){
      for(i in 1:length(typeC)){
        points(pos_s_real[i],  pos_u_real[i], bg = colCell[i], col = colCell[i],...)
      }
}else{
for(i in 1:length(typeC)){
    x <- typeC[i]
    points(pos_s_real[tipoCellula == x],  pos_u_real[tipoCellula == x], col = alpha(colCell[tipoCellula[tipoCellula == x]], 0.3), bg = alpha(colCell[tipoCellula[tipoCellula == x]], 0.3),...)
  }
}
}

}
  






plot_GSdynamics <- function(iterToPlot, t0_off_real, u0_off_real, s0_off_real, t0_on_real = 0, pos_u_real, pos_s_real, alpha_real, beta_real, gamma_real, u0_off_chain, s0_off_chain, u_chain, s_chain,alpha_chain, beta_chain, gamma_chain, subtypeCell, typeCell, g, alphaDef = 1, maxLOG = TRUE, sizePoints = 1,...){

  #################################
  ######### REAL STRUCUTRE ########
  #################################

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto  
  # ramo on
  u_SS_off_real <- alpha_real[1]/beta_real
  u_SS_on_real  <- alpha_real[2]/beta_real
  s_SS_off_real <- alpha_real[1]/gamma_real
  s_SS_on_real  <- alpha_real[2]/gamma_real
  
  # ramo on tutto   
  u_plotON_real <- seq(u_SS_off_real, u_SS_on_real, 0.01)
  s_plotON_real <- s_di_u(u_plotON_real, u_SS_on_real, s_SS_on_real, u_SS_off_real, s_SS_off_real, k = 2, alpha_real, beta_real, gamma_real)


  # ramo off
  u_plotOFF_real <- seq(u_SS_off_real, u_SS_on_real, 0.01)
  s_plotOFF_real <- s_di_u(u_plotOFF_real, u_SS_on_real, s_SS_on_real, u_SS_off_real, s_SS_off_real, k = 0, alpha_real, beta_real, gamma_real)

  df_real <- data.frame(uON_SS = as.vector(u_plotON_real), sON_SS = as.vector(s_plotON_real), uOFF_SS = as.vector(u_plotOFF_real), sOFF_SS = as.vector(s_plotOFF_real), gene = g) 

  pl <- ggplot(df_real)  + geom_path(aes(x = sON_SS, y = uON_SS), col = alpha("black", 1), linetype = "dotted", ...)+ geom_path(aes(x = sOFF_SS, y = uOFF_SS), col = alpha("black", 1), linetype = "dotted", ...) 
  

  #################################
  ######### RESULTS ###############
  #################################
  
  dfSS_chain <- data.frame(uON_SS = alpha_chain[2, iterToPlot]/beta_chain[iterToPlot], sON_SS = alpha_chain[2, iterToPlot]/gamma_chain[iterToPlot], uOFF_SS = alpha_chain[1, iterToPlot]/beta_chain[iterToPlot], sOFF_SS = alpha_chain[1, iterToPlot]/gamma_chain[iterToPlot])


  pl <- pl + geom_point(data = dfSS_chain, aes(x = sON_SS, y = uON_SS), col = alpha("blue", alphaDef), shape = 8, size = sizePoints) + geom_point(data = dfSS_chain, aes(x = sOFF_SS, y = uOFF_SS, color = "green", size = sizePoints), shape = 8, size = sizePoints)






  ##########################################
  if(length(t0_off_real) == 1){
    attr(u0_off_chain, "dim") <- c(1, length(u0_off_chain))
    attr(s0_off_chain, "dim") <- c(1, length(s0_off_chain))
  }
  for(sty in 1:length(t0_off_real)){

    #################################
    ######### REAL STRUCUTRE ########
    #################################
      if(t0_off_real[sty] != Inf){
        #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
        # ramo on 
        u_plotON_switch_real <- seq(u_SS_off_real, u0_off_real[sty], 0.01)

        
        s_plotON_switch_real <- s_di_u(u_plotON_switch_real, u0_off_real[sty], s0_off_real[sty], u_SS_off_real, s_SS_off_real, k = 2, alpha_real, beta_real, gamma_real)
        
              
        # ramo off 
        u_plotOFF_switch_real <- seq(u_SS_off_real, u0_off_real[sty], 0.01)
        s_plotOFF_switch_real <- s_di_u(u_plotOFF_switch_real, u0_off_real[sty], s0_off_real[sty], u_SS_off_real, s_SS_off_real, k = 0, alpha_real, beta_real, gamma_real)
        
        dfSwitchON_real <- data.frame(uON_switch = u_plotON_switch_real, sON_switch = s_plotON_switch_real)
        dfSwitchOFF_real <- data.frame(uOFF_switch = u_plotOFF_switch_real, sOFF_switch = s_plotOFF_switch_real)


        # dfSwitchON <- dfSwitchON[order(dfSwitchON$sON_switch), order(dfSwitchON$sON_switch)]

        pl <- pl + geom_path(data = dfSwitchON_real, aes(x = sON_switch, y = uON_switch), color = alpha("black", alphaDef), ...) + geom_path(data = dfSwitchOFF_real, aes(x = sOFF_switch, y = uOFF_switch), color = alpha("black", alphaDef), ...) 

      }else{
        pl <- pl + geom_path(aes(x = sON_SS, y = uON_SS, color = alpha("black", alphaDef)), ...) + geom_path(aes(x = sOFF_SS, y = uOFF_SS), color = alpha("black", alphaDef), ...)
        
      }
      
      #################################
      ############# RESULTS ###########
      #################################
      dfSwitch_chain <- data.frame(uOFF_switch = u0_off_chain[sty,iterToPlot], sOFF_switch = s0_off_chain[sty,iterToPlot])

      pl <- pl + geom_point(data = dfSwitch_chain, aes(x = sOFF_switch, y = uOFF_switch), col = alpha("green", alphaDef), shape = 8, size = sizePoints) 

  }


  u_chain = as.vector(u_chain[, iterToPlot])
  s_chain = as.vector(s_chain[, iterToPlot])

  dfPos_chain <- data.frame(u = u_chain, s = s_chain, sty = rep(seq(1, max(subtypeCell)), length(iterToPlot)))

  
  pl <- pl + geom_point(data = dfPos_chain, aes(x = s, y = u, color = sty), colour = alpha(rep(seq(1, max(subtypeCell)), length(iterToPlot)), alphaDef), shape = 20, size = sizePoints) 

  dfPos_real <- data.frame(u = pos_u_real, s = pos_s_real, sty = seq(1, max(subtypeCell)))

  if(maxLOG){
    sizePointsREAL <- sizePoints/2
  }else{
    sizePointsREAL <- 2
  }
  pl <- pl + geom_point(data = dfPos_real, aes(x = s, y = u, color = sty), fill = seq(1, max(subtypeCell)),  colour = "black", shape = 21, size = sizePointsREAL) 



  main <- paste("Gene", g)
#   if(!is.na(r)){
#     main <- paste(main, " for group r")
#   }

  pl <- pl + theme(legend.position = "none") + labs(title = paste(main), x = expression(s[g]), y = expression(u[g]), ...)
 
# + xlim(min(rbind(df_real$sON_SS, df_real$sOFF_SS)) - 1, max(rbind(df_real$sON_SS, df_real$sOFF_SS)) + 1) + ylim(min(rbind(df_real$uON_SS, df_real$uOFF_SS)) - 1, max(rbind(df_real$uON_SS, df_real$uOFF_SS)) + 1)

  return(pl)
}


# SIMILE A plot_sVSu_GGPLOT()
plot_GSdynamics2024 <- function(iterToPlot, t0_off_real, u0_off_real, s0_off_real, t0_on_real = 0, pos_u_real, pos_s_real, alpha_real, beta_real, gamma_real, u0_off_chain, s0_off_chain, u_chain, s_chain,alpha_chain, beta_chain, gamma_chain, subtypeCell, typeCell, g, alphaDef = 1, maxLOG = TRUE, sizePoints = 1,...){

  #################################
  ######### REAL STRUCUTRE ########
  #################################

  #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto  
  # ramo on
  u_SS_off_real <- alpha_real[1]/beta_real
  u_SS_on_real  <- alpha_real[2]/beta_real
  s_SS_off_real <- alpha_real[1]/gamma_real
  s_SS_on_real  <- alpha_real[2]/gamma_real
  
  # ramo on tutto   
  u_plotON_real <- seq(u_SS_off_real, u_SS_on_real, 0.01)
  s_plotON_real <- s_di_u(u_plotON_real, u_SS_on_real, s_SS_on_real, u_SS_off_real, s_SS_off_real, k = 2, alpha_real, beta_real, gamma_real)


  # ramo off
  u_plotOFF_real <- seq(u_SS_off_real, u_SS_on_real, 0.01)
  s_plotOFF_real <- s_di_u(u_plotOFF_real, u_SS_on_real, s_SS_on_real, u_SS_off_real, s_SS_off_real, k = 0, alpha_real, beta_real, gamma_real)

  df_real <- data.frame(uON_SS = as.vector(u_plotON_real), sON_SS = as.vector(s_plotON_real), uOFF_SS = as.vector(u_plotOFF_real), sOFF_SS = as.vector(s_plotOFF_real), gene = g) 

  pl <- ggplot(df_real)  + geom_path(aes(x = sON_SS, y = uON_SS), col = alpha("black", 1), linetype = "dotted")+ geom_path(aes(x = sOFF_SS, y = uOFF_SS), col = alpha("black", 1), linetype = "dotted") 

  #################################
  ######### RESULTS ########
  #################################
  u0_on_chain = alpha_chain[1]/beta_chain
  s0_on_chain = alpha_chain[1]/gamma_chain

  # ramo on tutto   
  u_plotON_chain <- seq(u0_on_chain, u0_off_chain, 0.01)
  s_plotON_chain <- s_di_u(u_plotON_chain, u0_off_chain, s0_off_chain, u0_off_chain, s0_on_chain, k = 2, alpha_chain, beta_chain, gamma_chain)

  # ramo off
  u_plotOFF_chain <- seq(u0_on_chain, u0_off_chain, 0.01)
  s_plotOFF_chain <- s_di_u(u_plotOFF_chain, u0_off_chain, s0_off_chain, u0_off_chain, s0_on_chain, k = 0, alpha_chain, beta_chain, gamma_chain)

  df_chain <- data.frame(uON_SS = as.vector(u_plotON_chain), sON_SS = as.vector(s_plotON_chain), uOFF_SS = as.vector(u_plotOFF_chain), sOFF_SS = as.vector(s_plotOFF_chain), gene = g) 


  pl <- pl +  geom_path(data = df_chain, aes(x = sON_SS, y = uON_SS), col = alpha("red", 1), linetype = "dotted")+ geom_path(aes(x = sOFF_SS, y = uOFF_SS), col = alpha("red", 1), linetype = "dotted") 



  # ##########################################
  # if(length(t0_off_real) == 1){
  #   attr(u0_off_chain, "dim") <- c(1, length(u0_off_chain))
  #   attr(s0_off_chain, "dim") <- c(1, length(s0_off_chain))
  # }
  # for(sty in 1:length(t0_off_real)){

  #   #################################
  #   ######### REAL STRUCUTRE ########
  #   #################################
  #     if(t0_off_real[sty] != Inf){
  #       #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto 
  #       # ramo on 
  #       u_plotON_switch_real <- seq(u_SS_off_real, u0_off_real[sty], 0.01)

        
  #       s_plotON_switch_real <- s_di_u(u_plotON_switch_real, u0_off_real[sty], s0_off_real[sty], u_SS_off_real, s_SS_off_real, k = 2, alpha_real, beta_real, gamma_real)
        
              
  #       # ramo off 
  #       u_plotOFF_switch_real <- seq(u_SS_off_real, u0_off_real[sty], 0.01)
  #       s_plotOFF_switch_real <- s_di_u(u_plotOFF_switch_real, u0_off_real[sty], s0_off_real[sty], u_SS_off_real, s_SS_off_real, k = 0, alpha_real, beta_real, gamma_real)
        
  #       dfSwitchON_real <- data.frame(uON_switch = u_plotON_switch_real, sON_switch = s_plotON_switch_real)
  #       dfSwitchOFF_real <- data.frame(uOFF_switch = u_plotOFF_switch_real, sOFF_switch = s_plotOFF_switch_real)


  #       # dfSwitchON <- dfSwitchON[order(dfSwitchON$sON_switch), order(dfSwitchON$sON_switch)]

  #       pl <- pl + geom_path(data = dfSwitchON_real, aes(x = sON_switch, y = uON_switch), color = alpha("black", alphaDef), ...) + geom_path(data = dfSwitchOFF_real, aes(x = sOFF_switch, y = uOFF_switch), color = alpha("black", alphaDef), ...) 

  #     }else{
  #       pl <- pl + geom_path(aes(x = sON_SS, y = uON_SS, color = alpha("black", alphaDef)), ...) + geom_path(aes(x = sOFF_SS, y = uOFF_SS), color = alpha("black", alphaDef), ...)
        
  #     }
      
  #     #################################
  #     ############# RESULTS ###########
  #     #################################
  #     if(!is.na(iterToPlot)){
  #       dfSwitch_chain <- data.frame(uOFF_switch = u0_off_chain[sty,iterToPlot], sOFF_switch = s0_off_chain[sty,iterToPlot])
  #     }else{
  #       dfSwitch_chain <- data.frame(uOFF_switch = u0_off_chain[sty], sOFF_switch = s0_off_chain[sty])
  #     }

  #     pl <- pl + geom_point(data = dfSwitch_chain, aes(x = sOFF_switch, y = uOFF_switch), col = alpha("green", alphaDef), shape = 8, size = sizePoints) 

  # }

  # if(!is.na(iterToPlot)){
  #   u_chain = as.vector(u_chain[, iterToPlot])
  #   s_chain = as.vector(s_chain[, iterToPlot])
  # }
  # dfPos_chain <- data.frame(u = u_chain, s = s_chain, sty = rep(seq(1, max(subtypeCell)), length(iterToPlot)))

  
  # pl <- pl + geom_point(data = dfPos_chain, aes(x = s, y = u, color = sty), colour = alpha(rep(seq(1, max(subtypeCell)), length(iterToPlot)), alphaDef), shape = 20, size = sizePoints) 

  # dfPos_real <- data.frame(u = pos_u_real, s = pos_s_real, sty = seq(1, max(subtypeCell)))

  # if(maxLOG){
  #   sizePointsREAL <- sizePoints/2
  # }else{
  #   sizePointsREAL <- 2
  # }
  # pl <- pl + geom_point(data = dfPos_real, aes(x = s, y = u, color = sty), fill = seq(1, max(subtypeCell)),  colour = "black", shape = 21, size = sizePointsREAL) 



  main <- paste("Gene", g)
#   if(!is.na(r)){
#     main <- paste(main, " for group r")
#   }

  pl <- pl + theme(legend.position = "none") + labs(title = paste(main), x = expression(s[g]), y = expression(u[g]), ...)
 
# + xlim(min(rbind(df_real$sON_SS, df_real$sOFF_SS)) - 1, max(rbind(df_real$sON_SS, df_real$sOFF_SS)) + 1) + ylim(min(rbind(df_real$uON_SS, df_real$uOFF_SS)) - 1, max(rbind(df_real$uON_SS, df_real$uOFF_SS)) + 1)

  return(pl)
}


# simile a plot_sVSu_scVelo_GGPLOT()
plot_GSdynamicsDATI_REALI <- function(loglik, u0_off_chain, s0_off_chain, u_chain, s_chain,alpha_chain, beta_chain, gamma_chain, subtypeCell, typeCell, g, alphaDef = 1, maxLOG = TRUE, sizePoints = 1,...){
 

  #################################
  ######### RESULTS ###############
  #################################
  if(maxLOG){
    iterToPlot <- which.max(apply(loglik, MARGIN = 2, FUN = sum))
  }else{
    iterToPlot <- seq(1, dim(loglik)[2])
  }

  if(is.null(dim(u0_off_chain))){
    dim(u0_off_chain) <- c(1, length(u0_off_chain))
    dim(s0_off_chain) <- c(1, length(s0_off_chain))
  }
  dfSS_chain <- data.frame(uON_SS = alpha_chain[2, iterToPlot]/beta_chain[iterToPlot], sON_SS = alpha_chain[2, iterToPlot]/gamma_chain[iterToPlot], uOFF_SS = alpha_chain[1, iterToPlot]/beta_chain[iterToPlot], sOFF_SS = alpha_chain[1, iterToPlot]/gamma_chain[iterToPlot])

  pl <- ggplot(dfSS_chain) + geom_point(data = dfSS_chain, aes(x = sON_SS, y = uON_SS), col = alpha("blue", alphaDef), shape = 8, size = sizePoints) + geom_point(data = dfSS_chain, aes(x = sOFF_SS, y = uOFF_SS, color = "green", size = sizePoints), shape = 8, size = sizePoints)


  for(sty in 1:dim(u0_off_chain)[1]){

      #################################
      ############# RESULTS ###########
      #################################
      dfSwitch_chain <- data.frame(uOFF_switch = u0_off_chain[sty,iterToPlot], sOFF_switch = s0_off_chain[sty,iterToPlot])

      pl <- pl + geom_point(data = dfSwitch_chain, aes(x = sOFF_switch, y = uOFF_switch), col = alpha("green", alphaDef), shape = 8, size = sizePoints) 

  }


  u_chain = as.vector(u_chain[, iterToPlot])
  s_chain = as.vector(s_chain[, iterToPlot])

  dfPos_chain <- data.frame(u = u_chain, s = s_chain, sty = rep(seq(1, max(subtypeCell)), length(iterToPlot)))

  
  pl <- pl + geom_point(data = dfPos_chain, aes(x = s, y = u, color = sty), colour = alpha(rep(seq(1, max(subtypeCell)), length(iterToPlot)), alphaDef), shape = 20, size = sizePoints) 


  main <- paste("Gene", g)
#   if(!is.na(r)){
#     main <- paste(main, " for group r")
#   }

  pl <- pl + theme(legend.position = "none") + labs(title = paste(main), x = expression(s[g]), y = expression(u[g]), ...)
 
# + xlim(min(rbind(df_real$sON_SS, df_real$sOFF_SS)) - 1, max(rbind(df_real$sON_SS, df_real$sOFF_SS)) + 1) + ylim(min(rbind(df_real$uON_SS, df_real$uOFF_SS)) - 1, max(rbind(df_real$uON_SS, df_real$uOFF_SS)) + 1)

  return(pl)
}



plot_GeneDynamic_withNotes <- function(t0_off, u0_off, s0_off, t0_on = 0, alpha, beta, gamma, r = NA, g){
  
   #----- (u_SS_basso, s_SS_basso) fino a (u_SS_alto, s_SS_alto) ----> partiamo dallo steady state basso e arriviamo a quello alto
  # ramo on
  u_SS_off <- alpha[1]/beta
  u_SS_on  <- alpha[2]/beta
  s_SS_off <- alpha[1]/gamma        
  s_SS_on  <- alpha[2]/gamma        

  # ramo on tutto
  u_plotON <- seq(u_SS_off, u_SS_on, 0.01)
  s_plotON <- s_di_u(u_plotON, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 2, alpha, beta, gamma)


  # ramo off
  u_plotOFF <- seq(u_SS_off, u_SS_on, 0.01)
  s_plotOFF <- s_di_u(u_plotOFF, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 0, alpha, beta, gamma)

  df <- data.frame(uON_SS = as.vector(u_plotON), sON_SS = as.vector(s_plotON), uOFF_SS = as.vector(u_plotOFF), sOFF_SS = as.vector(s_plotOFF), gene = g)

  pl <- ggplot(df)  +
    geom_path(aes(x = sON_SS, y = uON_SS, color = "red"), linetype = "dotted", linewidth = 3) + 
    geom_path(aes(x = sOFF_SS, y = uOFF_SS, color = "blue"), linetype = "dotted", linewidth = 3)

  if(t0_off != Inf){
    #----- (u0_on, s0_on) fino a (u0_off, s0_off)   ----> c'è lo switch prima di raggiungere lo steady state alto
    # ramo on
    u_plotON_switch <- seq(u_SS_off, u0_off, 0.001)
    s_plotON_switch <- s_di_u(u_plotON_switch, u0_off, s0_off, u_SS_off, s_SS_off, k = 2, alpha, beta, gamma)

    # ramo off
    u_plotOFF_switch <- seq(u_SS_off, u0_off, 0.001)
    s_plotOFF_switch <- s_di_u(u_plotOFF_switch, u0_off, s0_off, u_SS_off, s_SS_off, k = 0, alpha, beta, gamma)

    dfSwitchON <- data.frame(uON_switch = u_plotON_switch, sON_switch = s_plotON_switch)
    dfSwitchOFF <- data.frame(uOFF_switch = u_plotOFF_switch, sOFF_switch = s_plotOFF_switch)

    maxSoff <- which.max(s_plotOFF_switch)

    # dfSwitchON <- dfSwitchON[order(dfSwitchON$sON_switch), order(dfSwitchON$sON_switch)]

    pl <- pl +
      geom_path(data = dfSwitchON, aes(x = sON_switch, y = uON_switch, color = "red"), linewidth = 3) +     
      geom_path(data = dfSwitchOFF, aes(x = sOFF_switch, y = uOFF_switch, color = "blue"), linewidth = 3)   

  }else{
    pl <- pl +
      geom_path(aes(x = sON_SS, y = uON_SS, color = "red"), linewidth = 3) +
      geom_path(aes(x = sOFF_SS, y = uOFF_SS, color = "blue"), linewidth = 3)
  }

  if(is.na(g)){
    g <- "g"
  }
  main <- ""
  if(!is.na(r)){
    main <- paste(main, " for group r")
  }


  pl <- pl + 
        theme(legend.position = "none") + 
        annotate(geom = "text", y = alpha[1]*0.88, x = (alpha[1]/gamma)*0.9,  label = TeX("$SS^{off}$", output = "character"), size = 35, parse = TRUE,  family = "serif") +
        annotate(geom = "text", y = alpha[2]*1.05, x = (alpha[2]/gamma)*1.02, label = TeX("$SS^{on}$",  output = "character"), size = 35, parse = TRUE,  family = "serif") +
        labs(x = "", y = "") + theme(plot.title = element_text(family = "serif", size=70,  hjust = 0.5)) +
        annotate(geom = "text", y = u0_off*1.05, x = s0_off*0.92,
        # label =TeX(r"($(s(t_0^{on} + omega), u(t_0^{on} + omega))$)", output = "character"), 
        label =TeX(r"($(s^{omega}, u^{omega})$)", output = "character"), 
        size = 35, parse = TRUE, family = "serif") +
        annotate(geom = "text", x = min(df$sON_SS) + min(df$sON_SS)*0.45, y = min(df$uON_SS) + (max(df$uON_SS) - min(df$uON_SS))/2, label = "Induction", size = 35, angle = 62, family = "serif")  +
        annotate(geom = "text", x = max(df$sON_SS)*0.82, y = min(df$uON_SS) + (max(df$uON_SS) - min(df$uON_SS))/2 -0.1, label = "Repression", size = 35, angle = 56, family = "serif")

indexU1 <- which.min(abs(dfSwitchON$uON_switch - (min(dfSwitchON$uON_switch) + (max(dfSwitchON$uON_switch) -  min(dfSwitchON$uON_switch))/2)))  
provaU1 <-  dfSwitchON$uON_switch[indexU1]
provaS1 <- dfSwitchON$sON_switch[which(dfSwitchON$uON_switch == provaU1)[1]]
provaU2 <- dfSwitchON$uON_switch[indexU1 + 1]
provaS2 <- dfSwitchON$sON_switch[which(dfSwitchON$uON_switch == provaU2)[1]]

pl <- pl  + 
      geom_segment(aes(color = "red"), x = provaS1, y = provaU1,  xend =   provaS2, yend = provaU2, arrow = arrow( length = unit(0.3, "inches")), size = 3)


indexU1 <- which.min(abs(dfSwitchOFF$uOFF_switch - (min(dfSwitchOFF$uOFF_switch) + (max(dfSwitchOFF$uOFF_switch) -  min(dfSwitchOFF$uOFF_switch))/2)))
provaU1 <-  dfSwitchOFF$uOFF_switch[indexU1]
provaS1 <- dfSwitchOFF$sOFF_switch[which(dfSwitchOFF$uOFF_switch == provaU1)[1]]
provaU2 <- dfSwitchOFF$uOFF_switch[indexU1 + 1]
provaS2 <- dfSwitchOFF$sOFF_switch[which(dfSwitchOFF$uOFF_switch == provaU2)[1]]


pl <- pl + geom_segment(aes(color = "blue"), x = provaS2, y = provaU2, xend = provaS1, yend = provaU1, arrow = arrow( length = unit(0.3, "inches")), size = 3)       


indexU1 <- which.min(abs(df$uOFF_SS 
- (min(df$uOFF_SS) + (max(df$uOFF_SS) -  min(df$uOFF_SS))/2)))
provaU1 <-  df$uOFF_SS[indexU1]     
provaS1 <- df$sOFF_SS[which(df$uOFF_SS == provaU1)[1]]
provaU2 <- df$uOFF_SS[indexU1 + 5]  
provaS2 <- df$sOFF[which(df$uOFF_SS 
== provaU2)[1]]


pl <- pl + geom_segment(aes(colour = "blue"), x = provaS2, y = provaU2, 
 xend = provaS1, yend = provaU1,  arrow = arrow(length = unit(0.3, "inches")), size = 3)

pl <- pl +
  coord_cartesian(xlim =  c(min(rbind(df$sON_SS, df$sOFF_SS)) - 0.5, max(rbind(df$sON_SS, df$sOFF_SS)) + 0.5), ylim = c(min(rbind(df$uON_SS, df$uOFF_SS)) - 0.2, max(rbind(df$uON_SS, df$uOFF_SS)) + 0.2)) + theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())



pl <- pl + geom_segment(aes(x=3.7, y=0.85, xend=4.5, yend=0.85), arrow = arrow(length=unit(.5, 'cm')), color='black', linewidth=2) + 
geom_segment(aes(x=3.7, y=0.846, xend=3.7, yend=1.35), arrow = arrow(length=unit(.5, 'cm')), color='black', linewidth=2) + 
annotate(geom = "text", y = 0.78, x = 4.1,  label = TeX("$s$", output = "character"), size = 30, parse = TRUE,  family = "serif") +
annotate(geom = "text", y = 1.046, x = 3.55,  label = TeX("$u$", output = "character"), size = 30, parse = TRUE,  family = "serif") 

  return(pl)
}











ss_to_rates <- function(beta, u_off, s_off, diffU){
  u_on <- u_off + diffU
  alpha_off <- beta*u_off
  alpha_on  <- beta*u_on
  gamma     <- u_off/s_off * beta
  return(c(alpha_off, alpha_on, gamma))
  
}



# generiamo da p t.c.
# p(x) = p0*delta(x = 0) + p1*delta(x = 1) + (1-p0-p1)*U(0,1)
# Quindi assegnamo con probabilit� p0 il valore 0 (steady state OFF)
#                  con probabilit� p1 il valore 1 (steady state ON)
#                  con probabilit� , prev1-p0-p1 un valore proveniente da un unifore in (0,1)
rProporzioni <- function(p0, p1, p2, prev){
  type <- which(rmultinom(1, 1,  c(p0, p1,p2, 1-(p0+p1+p2)))==1)
  if(type == 1){
    res <- 0
  }else if(type == 2){
    res <- 1
  }else if(type == 3){
    res <- runif(1, 0, 1)
  }else{
    res <- runif(1,max(prev-0.1,0), min(prev+0.1, 1))
  }
  return(res)
}

# densit� con massa nei due steady- state
dProporzioni <- function(x, p0, p1, p2, prev){
  x = round(x, 5) 
  res <- p0*I(x == 0) + p1*I(x ==1) + p2*dunif(x, 0, 1)*I(x>0 & x<1) + (1-p0-p1-p2)*dunif(x, max(prev-0.1, 0), min(prev+0.1, 1))*I(x>0 & x<1)
  return(res)
}

# dalla proporzione otteniamo la coordinata tra uOFF e uON
proporz_to_coord <- function(prop, alpha, beta, gamma, typeCoord){
  print("Da sistemare con u0, to, s0")
  # if(typeCoord == "u"){
  #   den <- beta
  # }else if(typeCoord == "s"){
  #   den <- gamma
  # }
  # res <- prop*(alpha[2]-alpha[1])/den + alpha[1]/den
  # return(res)
}

coord_to_proporz <- function(coord, alpha, beta, gamma, typeCoord){
  print("Da risistemare con u0, s0, t0")
  # if(typeCoord == "u"){
  #   den <- beta
  # }else if(typeCoord == "s"){
  #   den <- gamma
  # }
  # res <- (coord- alpha[1]/den)/((alpha[2]-alpha[1])/den)
  # return(res)
}

# s/u to its symmetric (in the original coordinates)
simmetria_coord <- function(coord, alpha, beta, gamma, typeCoord){
  if(typeCoord == "u"){
    den <- beta
  }else if(typeCoord == "s"){
    den <- gamma
  }
  res <- sum(alpha)/den - coord
  return(res)
}

# assumiamo di essere in un modello con sigma2 = IG(shape_prior, rate_prior)
#                                       mu|sigma2 = N(mu_prior, sigma2/n0)
#                                       dato_i | mu, sigma2 = N(mu, sigma2) con i = 1, ..., n
# Allora otteniamo che le distribuzioni a posteriori sono 
#                                       sigma2|dati = IG(v1/2, v1*phi1/2)
#                                       mu|sigma2, dati = N(mu1, sigma2/n1)
# con specifici parametri 
r_invgamma_posterior <- function(nSample, shape_prior, rate_prior, mu_prior, n0, dati){
  # v0 <- 2*shape_prior
  # phi0 <- 2*rate_prior/v0
  # n <- length(dati)
  # n1 <- n0 + n
  # v1 <- v0 + n
  # s2 <- sum((dati - mean(dati))^2)/(n-1)
  # phi1 <- 1/v1*((n-1)*s2 + v0*phi0 + n*n0/n1*(mean(dati) - mu_prior)^2)
  # res <- invgamma:::rinvgamma(nSample, shape = v1/2, rate = v1*phi1/2)
  n <- length(dati)
  n_posterior <- n0 + n
  mu_posterior <- (n0*mu_prior + n*mean(dati))/n_posterior
  shape_posterior <- shape_prior + n/2
  rate_posterior <- rate_prior + 1/2*(mu_prior^2*n0 + sum(dati^2) - mu_posterior^2*n_posterior)
  res <- invgamma:::rinvgamma(x, shape = shape_posterior, rate = rate_posterior)
  
  return(res)
}
d_invgamma_posterior <- function(x, shape_prior, rate_prior, mu_prior, n0, dati, log = FALSE){
  n <- length(dati)
  n_posterior <- n0 + n
  mu_posterior <- (n0*mu_prior + n*mean(dati))/n_posterior
  shape_posterior <- shape_prior + n/2
  rate_posterior <- rate_prior + 1/2*(mu_prior^2*n0 + sum(dati^2) - mu_posterior^2*n_posterior)
  # s2 <- sum((dati-mean(dati))^2)
  # rate_posterior <- rate_prior + s2/2 + n0*n*(mean(dati)- mu_prior)^2/(2*(n0 + n))
  res <- invgamma:::dinvgamma(x, shape = shape_posterior, rate = rate_posterior, log = log)
  
  return(res)
}

r_norm_posterior <- function(nSample, mu_prior, n0, sigma2, dati){
  n <- length(dati)
  n_posterior <- n0 + n
  mu_posterior <- (n0*mu_prior + n*mean(dati))/n_posterior
  res <- rnorm(nSample, mean = mu_posterior, sd = sqrt(sigma2/n_posterior))
  return(res)
}
d_norm_posterior <- function(x, mu_prior, n0, sigma2, dati, log = FALSE){
  n <- length(dati) 
  n_posterior <- n0 + n
  mu_posterior <- (n0*mu_prior + n*mean(dati))/n_posterior
  res <- dnorm(x, mean = mu_posterior, sd = sqrt(sigma2/n_posterior), log = log)
  return(res)
}

# Assumiamo che p = Beta(alpha, beta)
#               dato_i = bern(p) con i = 1, ..., n
# Allora la distribusione a posteriori di p � 
#               p | dati = Beta(alpha + sum(dati), beta + n - sum(dati))
r_beta_posterior <- function(nSample, alpha_prior, beta_prior, dati){
  n <- length(dati)
  successi <- sum(dati)
  res <- rbeta(nSample, shape1 = alpha_prior + successi, shape2 = beta_prior + n - successi)
  return(res)
}
d_beta_posterior <- function(x, alpha_prior, beta_prior, dati, log = FALSE){
  n <- length(dati)
  successi <- sum(dati)
  res <- dbeta(x, shape1 = alpha_prior + successi, shape2 = beta_prior + n - successi)
  if(log == TRUE){
    res <- log(res)
  }
  return(res)
}


# SPIEGARE CHE FUNZIONE E' 
g_tau <- function(x, logINPUT = TRUE, logOUTPUT = TRUE, k, alpha, beta, gamma){
  if(logINPUT == TRUE){
    x <- exp(x)
  }
  
  res <- -log1p(-3/4*exp(-x))
  if(logOUTPUT == TRUE){
    res <- log(res)
  }
  
  res[res < log(1e-5)] <- log(1e-5)
  res[res > log(1e5)]  <- log(1e5)
  return(res)
}

dg_tau <- function(x){
 res <- log(3/4) + x - exp(x) - log1p(-3/4*exp(-exp(x))) - log(abs(log1p(-3/4*exp(-exp(x)))))
  return(res)
}



plot_variabilita <- function(){
  tau_Q <- apply(exp(tau_log[g,,]), MARGIN = 2, quantile, c(0.25, 0.75))
  uOFF_Q <- quantile(exp(uOFF_log), c(0.25, 0.75))
  sOFF_Q <- quantile(exp(sOFF_log), c(0.25, 0.75))
  diffUOFF_Q <- quantile(exp(diffU_log), c(0.25, 0.75))
  gamma_Q <- quantile(gamma[g], c(0.25, 0.75))
  
  alphaOFF_Q <- quantile(alpha[g,,1], c(0.25, 0.75))
  alphaON_Q <- quantile(alpha[g,,2], c(0.25, 0.75))
  
  uON_Q <- quantile(alpha[g,,2]/beta[g], c(0.25, 0.75))
  sON_Q <- quantile(alpha[g,,2]/gamma[g], c(0.25, 0.75))
    
  library(plotrix)
  
  plotCI(meanSS_u[g], ui=uOFF_Q[2], li=uOFF_Q[1])
  plotCI(meanSS_s[g], ui=sOFF_Q[2], li=sOFF_Q[1], err = "x")
  plotCI(meanSS_u[g] + meandDiff[g], ui = uON_Q[2], li = uON_Q[1])
  plotCI(meanSS_u[g] + meandDiff[g], ui = uON_Q[2], li = uON_Q[1], err = "x")
  polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
  polygon(
    c(u(seq(0, 30, 0.01), k = rep(2, length(seq(0, 30, 0.01))), c(alphaOFF_Q[1], alphaON_Q[1]), 1), u(seq(0, 30, 0.01), k = rep(2, length(seq(0, 30, 0.01))), c(alphaOFF_Q[2], alphaON_Q[2]), 1)),
    c(s(seq(0, 30, 0.01), k = rep(2, length(seq(0, 30, 0.01))), c(alphaOFF_Q[1], alphaON_Q[1]), 1, gamma_Q[1]), s(seq(0, 30, 0.01), k = rep(2, length(seq(0, 30, 0.01))), c(alphaOFF_Q[2], alphaON_Q[2]), 1, gamma_Q[2]))
    , col = "grey", border = TRUE)
  lines(u(seq(0, 30, 0.01), k = rep(0, length(seq(0, 30, 0.01))), meanAlpha[g,], 1), s(seq(0, 30, 0.01), k = rep(0, length(seq(0, 30, 0.01))), meanAlpha[g,], 1, meanGamma[g]), type = "l", col = "red", lwd = 1)
}




#################################### ANALYSIS RESULTS functions #################################
moda <- function(v) {
   value <- c(0, 2)
   return(value[which.max(tabulate(match(v, value)))])
}

nameGenes <- function(path){
    var <- read.csv(paste(path, "/var.csv", sep = ""))
    nameGenes <- var$index
    return(nameGenes)
}

nameCells <- function(path){
    obs <- read.csv(paste(path, "/obs.csv", sep = ""))
    nameCells <- obs$index
    return(nameCells)
}

nameTypeCells <- function(path){
    obs <- read.csv(paste(path, "/obs.csv", sep = ""))
    typeC_string <- unique(obs$clusters)
    typeC_numeric <- unique(as.numeric(as.factor(obs$clusters)))

    res <- list()
    for(i in 1:length(typeC_numeric)){
        res[[as.character(typeC_numeric[i])]] <- typeC_string[i]
    }
    return(res)
}

driverGenes <- function(driverScVelo, nameGenes){
    index <- rep(NA, length(driverScVelo))
    for(i in 1:length(index)){
        index[i] <- which(nameGenes == driverScVelo[i])
    }
    return(index)
}

riscalChain <- function(LogSS_chain, LogitCatt_chain){
    maxCatt_perIter = rep(NA, dim(LogitCatt_chain)[2])
    LogSS_chainRiscal = LogSS_chain
    Catt_chainRiscal = invlogit(LogitCatt_chain)
    for(iter in 1:length(maxCatt_perIter)){
        maxCatt_perIter[iter] = max(invlogit(LogitCatt_chain[,iter]))
        Catt_chainRiscal[,iter] = Catt_chainRiscal[,iter]/maxCatt_perIter[iter]
        LogSS_chainRiscal[1:3, ,iter] = LogSS_chainRiscal[1:3, ,iter] + log(maxCatt_perIter[iter])
    }
    return(list("LogSS_chainRiscal" = LogSS_chainRiscal, "Catt_chainRiscal" = Catt_chainRiscal))
}


################################ PLOT ################################
histConteggi <- function(Y_u, Y_s, nameGenes, geni){
    for(g in geni){
        histConteggi_gene(Y_u, Y_s, nameGenes, g)
    }
}

histConteggi_gene <- function(Y_u, Y_s, nameGenes, g){
    hist(Y_u[,g], xlab = "unspliced", breaks = seq(min(Y_u[,g])-0.5, max(Y_u[,g])+0.5), xlim = c(0, quantile(Y_u[,g], 0.90)), main = paste("Unspliced gene", g, nameGenes[g]))
    hist(Y_s[,g], xlab = "spliced"  , breaks = seq(min(Y_s[,g])-0.5, max(Y_s[,g])+0.5), xlim = c(0, quantile(Y_s[,g], 0.90)), main = paste("Spliced gene", g, nameGenes[g]))
}

identifySS <- function(dictTau, dictk, dictT0_off, typeCell){
    for(i in 1:length(dictTau)){
        SS1 <- which(dictk[[i]] == 0 & dictTau[[i]] >= 1000, arr.ind = TRUE)
        dictk[[i]][SS1] <- 2
        dictTau[[i]][SS1] <- 0

        # T0_off_rep <- rep(NA, prod(dim(dictTime[[i]])))
        # attr(T0_off_rep, "dim") <- c(dim(dictTime[[i]]))
        # n_typeC <- max(typeCell)

        # for(ty in 1:n_typeC){
        #   cellTy <- which(typeCell == ty)
        #   print(cellTy)
        #   print(ty)
        #   print(dim(T0_off_rep))
        #   print(dim(dictT0_off[[i]]))
        #   T0_off_rep[cellTy, ,] <- dictT0_off[[i]][ty, , ]          
        # }

        SS2 <- which(dictk[[i]] == 2 & dictTau[[i]] >= 1000, arr.ind = TRUE)
        dictk[[i]][SS2] <- 0
        dictTau[[i]][SS2] <- 0
        
    }
    res <- list("Tau" = dictTau, "k" = dictk)
    return(res)
}

plotLogSS_chain <- function(dictLogSS, u_SSoff_real, s_SSoff_real, diffSS_off_real, nameGenes, geni, plotBeta){
    for(g in geni){
        plotLogSS_chain_gene(dictLogSS, u_SSoff_real, s_SSoff_real, diffSS_off_real, nameGenes, g, plotBeta)
    }
}

plotLogSS_chain_gene <- function(LogSS_chain, u_SSoff_real, s_SSoff_real, diffSS_off_real, g, plotBeta){
    namesSS <- c("Log-uOFF", "Log_sOFF", "Log-diffU", "Log-Beta")
    real <- log(rbind(u_SSoff_real, s_SSoff_real, diffSS_off_real, rep(1, length(diffSS_off_real))))

    n_comp <- dim(LogSS_chain)[1]
    if(!plotBeta){
      n_comp <- n_comp - 1
    }

    for(j in 1:n_comp){
        plot(LogSS_chain[j,g,], type="l", ylab = namesSS[j], main = paste(namesSS[j], round(real[j,g], 3), ", g", g), ylim = range(LogSS_chain[j, g,], real[j, g]))
        lines(x = seq(1,mcmcIter), y = rep(mean(LogSS_chain[j,g,]), mcmcIter), type = "l", col = "green")
        lines(x = seq(1,mcmcIter), y = rep(median(LogSS_chain[j,g,]), mcmcIter), type = "l", col = "purple")
        abline(h = round(real[j,g], 3), col = "red", lwd = 2.5)

    }    
}


plotSS_chain <- function(SS_chain, u_SSoff_real, s_SSoff_real, diffSS_off_real, geni){
    for(g in geni){
        plotSS_chain_gene(SS_chain, u_SSoff_real, s_SSoff_real, diffSS_off_real, g)
    }
}

plotSS_chain_gene <- function(SS_chain, u_SSoff_real, s_SSoff_real, diffSS_off_real, g){
    namesSS <- c("uOFF", "sOFF", "diffU")
    real <- rbind(u_SSoff_real, s_SSoff_real, diffSS_off_real, rep(1, length(diffSS_off_real)))

    n_comp <- dim(SS_chain)[1] - 1
   
    for(j in 1:n_comp){
        plot(SS_chain[j,g,], type="l", ylab = namesSS[j], main = paste(namesSS[j], round(real[j,g], 3), ", g", g), ylim = range(SS_chain[j, g,], real[j, g]))
        abline(h = round(real[j,g], 3), col = "red", lwd = 2.5)
    }    
}

plotRate_chain <- function(alpha_chain, gamma_chain, alpha_real, gamma_real){
    for(g in geni){
        plotRate_chain_gene(alpha_chain, gamma_chain, alpha_real, gamma_real, g)
    }
}

plotRate_chain_gene <- function(alpha_chain, gamma_chain, alpha_real, gamma_real, g){
    nameRates <- c("Alpha_off", "Alpha_on", "Gamma")


    real <- rbind(alpha_real[,1], alpha_real[,2], gamma_real, rep(1, length(gamma_real)))
   
    plot(alpha_chain[g,1,], type="l", ylab = nameRates[1], main = paste(nameRates[1], round(real[1,g], 3), ", g", g), ylim = range(alpha_chain[g,1,], real[1, g]))
    abline(h = round(real[1,g], 3), col = "red", lwd = 2.5)

    plot(alpha_chain[g,2,], type="l", ylab = nameRates[2], main = paste(nameRates[2], round(real[2,g], 3), ", g", g), ylim = range(alpha_chain[g,2,], real[2, g]))
    abline(h = round(real[2,g], 3), col = "red", lwd = 2.5)

    plot(gamma_chain[g,], type="l", ylab = nameRates[3], main = paste(nameRates[3], round(real[3,g], 3), ", g", g), ylim = range(gamma_chain[g,], real[3, g]))
    abline(h = round(real[3,g], 3), col = "red", lwd = 2.5)
    
}

plotLogTime_chain <- function(dictLogTau, LogTau_real, nameGenes, geni, lb, ub, nameTypeCells, typeCell, typeCellReal){
    for(g in geni){
        plotLogTime_chain_gene(dictLogTau, LogTau_real, nameGenes, g, lb, ub, nameTypeCells, typeCell, typeCellReal)
    }
}

plotLogTime_chain_gene <- function(dictLogTau, LogTau_real, nameGenes, g, lb, ub, nameTypeCells, typeCell, typeCellReal){
    for(ty in 1:dim(dictLogTau[[1]])[1]){
      plotLogTime_chain_gene_ty(dictLogTau, LogTau_real, nameGenes, g, lb, ub, nameTypeCells, typeCell, typeCellReal, ty)
    }  
}


plotLogTime_chain_gene_ty <- function(dictLogTau, LogTau_real, nameGenes, g, lb, ub, nameTypeCells, typeCell,typeCellReal, ty){
  
      ty_real = unique(typeCellReal[which(typeCell == ty)])
      cellTyReal = which(typeCellReal == ty_real)
        for(i in 1:length(dictLogTau)){
            plot(dictLogTau[[i]][ty,g,], type="l", ylab = "LogTau_star",main = "")
            title(main = paste( "LogTau_star", round(LogTau_real[cellTyReal[1], g],3), ",g", g, ",ty", ty, ",tyReal", ty_real) ,cex.main = 1)
            lines(x = seq(1,mcmcIter), y = rep((mean(dictLogTau[[i]][ty, g,])), mcmcIter), type = "l", col = "green") 
            lines(x = seq(1,mcmcIter), y = rep((median(dictLogTau[[i]][ty, g,])), mcmcIter), type = "l", col = "purple") 
            abline(h = round(LogTau_real[cellTyReal[1],g], 3), col = "red", lwd = 1.5)
            lines(x = seq(1,mcmcIter), y = rep(lb, mcmcIter), type = "l", lty = 2, col = "blue") 
            lines(x = seq(1,mcmcIter), y = rep(ub, mcmcIter), type = "l", lty = 2, col = "blue") 
        }
}

plotLogT0_off_chain <- function(dictLogT0_off, LogT0_off_real, nameGenes, geni, lb, ub, nameTypeClusters, typeCellT0_off, typeCellT0_offReal){
    for(g in geni){
        plotLogT0_off_chain_gene(dictLogT0_off, LogT0_off_real, nameGenes, g, lb, ub, nameTypeClusters, typeCellT0_off, typeCellT0_offReal)
    }
}

plotLogT0_off_chain_gene <- function(dictLogT0_off, LogT0_off_real,  nameGenes, g, lb, ub, nameTypeClusters, typeCellT0_off, typeCellT0_offReal){
    for(tyT0_off in 1:dim(dictLogT0_off[[1]])[1]){
       plotLogT0_off_chain_gene_tyT0(dictLogT0_off, LogT0_off_real,  nameGenes, g, lb, ub, nameTypeClusters, tyT0_off, typeCellT0_off, typeCellT0_offReal)
    }  
}

plotLogT0_off_chain_gene_tyT0 <- function(dictLogT0_off, LogT0_off_real,  nameGenes, g, lb, ub, nameTypeClusters, tyT0_off, typeCellT0_off, typeCellT0_offReal){
  for(i in 1:length(dictLogT0_off)){
            tyT0_real = unique(typeCellT0_offReal[which(typeCellT0_off == tyT0_off)])
            plot(dictLogT0_off[[i]][tyT0_off,g,], type="l", ylab = "LogT0_off",main = "")
            title(main = paste("LogT0", round(LogT0_off_real[tyT0_real, g], 3), ", g", g, ", tyT0 sim", tyT0_off, ", tyT0_Real" ,tyT0_real), cex.main = 1)
            lines(x = seq(1,mcmcIter), y = rep((mean(dictLogT0_off[[i]][tyT0_off, g,])), mcmcIter), type = "l", col = "green") 
            lines(x = seq(1,mcmcIter), y = rep((median(dictLogT0_off[[i]][tyT0_off, g,])), mcmcIter), type = "l", col = "purple") 

            
            abline(h = round(LogT0_off_real[tyT0_real,g], 3), col = "red", lwd = 1.5)
            lines(x = seq(1,mcmcIter), y = rep(lb, mcmcIter), type = "l", lty = 2, col = "blue") 
            lines(x = seq(1,mcmcIter), y = rep(ub, mcmcIter), type = "l", lty = 2, col = "blue") 
        }
}




plotT0_off_chain <- function(T0_off_chain, t0_off_real, geni, typeCellT0_off){
    for(g in geni){
        plotT0_off_chain_gene(T0_off_chain, t0_off_real, g, typeCellT0_off)
    }
}

plotT0_off_chain_gene <- function(T0_off_chain, t0_off_real, g, typeCellT0_off){
  for(tyT0_off in 1:dim(T0_off_chain)[1]){
      plotT0_off_chain_gene_tyT0(T0_off_chain, t0_off_real, g, tyT0_off, typeCellT0_off)
  }  
}


plotT0_off_chain_gene_tyT0 <- function(T0_off_chain, t0_off_real, g, tyT0_off, typeCellT0_off){
  T0_off_chain[tyT0_off,g,which(T0_off_chain[tyT0_off,g,] == Inf)] <- 1000
  plot(T0_off_chain[tyT0_off,g,], type="l", ylab = "T0_off", main = "", ylim = range(T0_off_chain[tyT0_off,g,], t0_off_real[tyT0_off,g]))
  title(main = paste("T0_off", round(t0_off_real[tyT0_off, g],3), "g", g, ", tyT0 ", tyT0_off) ,cex.main = 1)
  abline(h = round(t0_off_real[tyT0_off,g], 3), col = "red", lwd = 1.5)
}

plotU0_off_chain_gene <- function(u0_off_chain, u0_off_real, g, typeCellT0_off){
  for(tyT0_off in 1:dim(u0_off_chain)[1]){
      plotU0_off_chain_gene_tyT0(u0_off_chain, u0_off_real, g, tyT0_off, typeCellT0_off)
  }  
}

plotU0_off_chain_gene_tyT0 <- function(u0_off_chain, u0_off_real, g, tyT0_off, typeCellT0_off){
  plot(u0_off_chain[tyT0_off,g,], type="l", ylab = "u0_off", main = "",  ylim = range(u0_off_chain[tyT0_off,g,], u0_off_real[tyT0_off,g]))
  title(main = paste("u0_off", round(u0_off_real[tyT0_off, g],3), "g", g, ", tyT0 ", tyT0_off) ,cex.main = 1)
  abline(h = round(u0_off_real[tyT0_off,g], 3), col = "red", lwd = 1.5)
}


plotS0_off_chain_gene <- function(u0_off_chain, u0_off_real, g, typeCellT0_off){
  for(tyT0_off in 1:dim(u0_off_chain)[1]){
      plotU0_off_chain_gene_tyT0(u0_off_chain, u0_off_real, g, tyT0_off, typeCellT0_off)
  }  
}


plotS0_off_chain_gene_tyT0 <- function(s0_off_chain, s0_off_real, g, tyT0_off, typeCellT0_off){
  plot(s0_off_chain[tyT0_off,g,], type="l", ylab = "s0_off", main = "",  ylim = range(s0_off_chain[tyT0_off,g,], s0_off_real[tyT0_off,g]))
  title(main = paste("s0_off", round(s0_off_real[tyT0_off, g],3), "g", g, ", tyT0 ", tyT0_off) ,cex.main = 1)
  abline(h = round(s0_off_real[tyT0_off,g], 3), col = "red", lwd = 1.5)
}




plotTime_chain <- function(T_chain, t_real, geni, subtypeCell){
  for(g in geni){
      plotTime_chain_gene(T_chain, t_real, g, subtypeCell)
  }
}

plotTime_chain_gene <- function(T_chain, t_real, g, subtypeCell){
  for(sty in 1:dim(T_chain)[1]){
    plotTime_chain_gene_sty(T_chain, t_real, g, sty, subtypeCell)
  }  
}

plotTime_chain_gene_sty <- function(T_chain, t_real, g, sty, subtypeCell){
  subcellTy <- which(subtypeCell == sty)

  if(dim(T_chain)[1] == dim(t_real)[1]){
    plot(T_chain[sty,g,], type="l", ylab = "T", main = "", ylim = range(T_chain[sty, g, ], t_real[sty, g]))
    title(main = paste( "T", round(t_real[sty, g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
    abline(h = round(t_real[sty,g], 3), col = "red", lwd = 1.5)
  }else{
    plot(T_chain[sty,g,], type="l", ylab = "T", main = "", ylim = range(T_chain[sty, g, ], t_real[subcellTy[1], g]))
    title(main = paste( "T", round(t_real[subcellTy[1], g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
    abline(h = round(t_real[subcellTy[1],g], 3), col = "red", lwd = 1.5)
  }
  
}


plotPos_u_chain <- function(pos_u_chain, pos_u_real, geni, subtypeCell){
  for(g in geni){
      plotPos_u_chain_gene(pos_u_chain, pos_u_real, g, subtypeCell)
  }
}

plotPos_u_chain_gene <- function(pos_u_chain, pos_u_real, g, subtypeCell){
  for(sty in 1:dim(pos_u_chain)[1]){
    plotPos_u_chain_gene_sty(pos_u_chain, pos_u_real, g, sty, subtypeCell)
  }  
}

plotPos_u_chain_gene_sty <- function(pos_u_chain, pos_u_real, g, sty, subtypeCell){
  subcellTy <- which(subtypeCell == sty)

  if(dim(pos_u_chain)[1] == dim(pos_u_real)[1]){
    plot(pos_u_chain[sty,g,], type="l", ylab = "U", main = "", ylim = range(pos_u_chain[sty, g, ], pos_u_real[sty, g]))
    title(main = paste( "U", round(pos_u_real[sty, g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
    abline(h = round(pos_u_real[sty,g], 3), col = "red", lwd = 1.5)
  }else{
    plot(pos_u_chain[sty,g,], type="l", ylab = "U", main = "", ylim = range(pos_u_chain[sty, g, ], pos_u_real[subcellTy[1], g]))
    title(main = paste( "U", round(pos_u_real[subcellTy[1], g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
    abline(h = round(pos_u_real[subcellTy[1],g], 3), col = "red", lwd = 1.5)
  }
  

}


plotPos_s_chain <- function(pos_s_chain, pos_s_real, geni, subtypeCell){
  for(g in geni){
      plotPos_s_chain_gene(pos_s_chain, pos_s_real, g, subtypeCell)
  }
}

plotPos_s_chain_gene <- function(pos_s_chain, pos_s_real, g, subtypeCell){
  for(sty in 1:dim(pos_s_chain)[1]){
    plotPos_s_chain_gene_sty(pos_s_chain, pos_s_real, g, sty, subtypeCell)
  }  
}

plotPos_s_chain_gene_sty <- function(pos_s_chain, pos_s_real, g, sty, subtypeCell){
  subcellTy <- which(subtypeCell == sty)  

    if(dim(pos_s_chain)[1] == dim(pos_s_real)[1]){
      plot(pos_s_chain[sty,g,], type="l", ylab = "U", main = "", ylim = range(pos_s_chain[sty, g, ], pos_s_real[sty, g]))
      title(main = paste( "S", round(pos_s_real[sty, g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
      abline(h = round(pos_s_real[sty,g], 3), col = "red", lwd = 1.5)
    }else{
      plot(pos_s_chain[sty,g,], type="l", ylab = "U", main = "", ylim = range(pos_s_chain[sty, g, ], pos_s_real[subcellTy[1], g]))
      title(main = paste( "S", round(pos_s_real[subcellTy[1], g], 3), ", g", g, ", sty", sty) ,cex.main = 1)
      abline(h = round(pos_s_real[subcellTy[1],g], 3), col = "red", lwd = 1.5)
    }
  

}

plotMu_chain_gene_ty <- function(dictMu, tau_real, nameGenes, g, nameTypeCells, typeCell, typeCellReal, ty){
      
      ty_real = unique(typeCellReal[which(typeCell == ty)])
      cellTyReal = which(typeCellReal == ty_real)
        for(i in 1:length(dictTau)){
            plot(dictMu[[i]][ty,g,], type="l", ylab = "Tau",main = "")
            title(main = paste( "Mu", round(tau_real[cellTyReal[1], g], 3), "g", g, ",ty", ty, ",tyReal", ty_real) ,cex.main = 1)
            lines(x = seq(1,mcmcIter), y = rep((mean(dictMu[[i]][ty, g,])), mcmcIter), type = "l", col = "green") 
            lines(x = seq(1,mcmcIter), y = rep((median(dictMu[[i]][ty, g,])), mcmcIter), type = "l", col = "purple") 
            abline(h = round(tau_real[cellTyReal[1],g], 3), col = "red", lwd = 1.5)
        }
}

plotSigma_chain_gene_ty <- function(dictSigma, nameGenes, g, nameTypeCells, typeCell, typeCellReal, ty){
      
      ty_real = unique(typeCellReal[which(typeCell == ty)])
      cellTyReal = which(typeCellReal == ty_real)
        for(i in 1:length(dictTau)){
            plot(dictSigma[[i]][ty,g,], type="l", ylab = "Tau",main = "")
            title(main = paste( "Sigma, g", g, ",ty", ty, ",tyReal", ty_real) ,cex.main = 1)
            lines(x = seq(1,mcmcIter), y = rep((mean(dictSigma[[i]][ty, g,])), mcmcIter), type = "l", col = "green") 
            lines(x = seq(1,mcmcIter), y = rep((median(dictSigma[[i]][ty, g,])), mcmcIter), type = "l", col = "purple") 
            
        }
}


plotk_chain <- function(dictk, k_real, nameGenes, geni, nameTypeCells, typeCell, typeCellReal){
    for(g in geni){
        plotk_chain_gene(dictk, k_real, nameGenes, g, nameTypeCells, typeCell, typeCellReal)
    }
}

# ANCORA DA AGGIUNGERE I NOMI GIUSTI DEI TIPI DI CELLULA
plotk_chain_gene <- function(dictk, k_real, nameGenes, g, nameTypeCells, typeCell, typeCellReal){
    for(ty in 1:dim(dictk[[1]])[1]){
      plotk_chain_gene_ty(dictk, k_real, nameGenes, g, nameTypeCells, typeCell, typeCellReal, ty)
    }  
}

plotk_chain_gene_ty <- function(dictk, k_real, nameGenes, g, nameTypeCells, typeCell, typeCellReal, ty){
  
      ty_real = unique(typeCellReal[which(typeCell == ty)])
      cellTyReal = which(typeCellReal == ty_real)
        
      for(i in 1:length(dictk)){


            plot(dictk[[i]][ty,g,], type="l", ylab = "k",main = "")
            title(main = paste("k", k_real[cellTyReal[1], g], ",g", g, ",ty", ty,",tyReal", ty_real ) ,cex.main = 1)
            abline(h = k_real[cellTyReal[1],g], col = "red", lwd = 1.5)

            # lines(x = seq(1,mcmcIter), y = rep((mean(dictk[[i]][ty, g,])), mcmcIter), type = "l", col = "green") 
            # lines(x = seq(1,mcmcIter), y = rep((median(dictk[[i]][ty, g,])), mcmcIter), type = "l", col = "purple") 
        }

}

plotLogEta_chain <- function(dictLogEta, eta_real, nameGenes, geni){
    for(g in geni){
        plotLogEta_chain_gene(dictLogEta,eta_real,  nameGenes, g)
    }
}

plotLogEta_chain_gene <- function(dictLogEta, eta_real, nameGenes, g){
    for(i in 1:length(dictLogEta)){
        plot(dictLogEta[[i]][g,], type="l", ylab = "LogEta",main = "")
        title(main = paste( "LogEta", round(log(eta_real[g]), 3), "g", g, nameGenes[g]) ,cex.main = 1)
        lines(x = seq(1,mcmcIter), y = rep((mean(dictLogEta[[i]][g,])), mcmcIter), type = "l", col = "green")  
        lines(x = seq(1,mcmcIter), y = rep((median(dictLogEta[[i]][g,])), mcmcIter), type = "l", col = "purple") 
        abline(h = round(log(eta_real[g]), 3), col = "red", lwd = 1.5)

    }
}



plotEta_chain <- function(Eta_chain, eta_real, geni){
  for(g in geni){
      plotEta_chain_gene(Eta_chain, eta_real,  g)
  }
}

plotEta_chain_gene <- function(Eta_chain, eta_real, g){
  plot(Eta_chain[g,], type="l", ylab = "Eta",main = "", ylim = range(Eta_chain[g,], eta_real[g]))
  title(main = paste("Eta", round(eta_real[g], 3), "g", g) ,cex.main = 1)
  abline(h = round(eta_real[g], 3), col = "red", lwd = 1.5)
}

plotCatt_chain <- function(Catt_chain, catt_real, cells){
  for(c in cells){
      plotCatt_chain_cell(Catt_chain, catt_real, c)
  }
}


plotCatt_chain_cell <- function(Catt_chain, catt_real, c){
  plot(Catt_chain[c,], type = "l", ylab = "Catt", main = paste(round(catt_real[c], 3),"c", c))
  abline(h = round(catt_real[c], 3), col = "red", lwd = 1.5)
}


plotLogitCatt_chain <- function(dictLogitCatt, catt_real, nameCells){
        for(c in 1:dim(dictLogitCatt[[1]])[1]){
            plotLogitCatt_chain_cell(dictLogitCatt, catt_real, nameCells, c)
        }
}


plotLogitCatt_chain_cell <- function(dictLogitCatt, catt_real, nameCells, c){
    for(i in 1:length(dictLogitCatt)[1]){
        plot(dictLogitCatt[[i]][c,], type = "l", ylab = "LogitCatt", main = paste( round(logit(catt_real[c]), 3),"c", c, nameCells[c]))
        lines(x = seq(1,mcmcIter), y = rep((mean(logit(dictCatt[[i]][c,]))), mcmcIter), type = "l", col = "green")  
        lines(x = seq(1,mcmcIter), y = rep((median(logit(dictCatt[[i]][c,]))), mcmcIter), type = "l", col = "purple")  
        abline(h = round(logit(catt_real[c]), 3), col = "red", lwd = 1.5)

    }
}


phasePlotMEAN <- function(dictLogSS, dictTime, dictk, typeC, nameGenes, geni){
    for(g in geni){
        phasePlotMEAN_gene(dictLogSS, dictTime, dictk, typeC, nameGenes, g)
    }
}

phasePlotMEAN_gene_ty <- function(dictLogSS, dictT0_off, dictTau, dictk, alpha_real, beta_real, gamma_real, tau_real, t0_off_real, k_real, nameGenes, g, ty,  typeCell, typeCellT0_off, SS_ON){  
    pchPoss = seq(15,(15 + length(dictLogSS))) 
    
    for(i in 1:length(dictLogSS)){
        meanAlphaOFF <- mean(exp(dictLogSS[[i]][1,g,])*exp(dictLogSS[[i]][4,g,]))
        meanAlphaON  <- mean(exp(dictLogSS[[i]][4,g,])*(exp(dictLogSS[[i]][3,g,]) + exp(dictLogSS[[i]][1,g,])))
        meanBeta     <- mean(exp(dictLogSS[[i]][4,g,]))
        meanGamma    <- mean(exp(dictLogSS[[i]][4,g,])*exp(dictLogSS[[i]][1,g,])/exp(dictLogSS[[i]][2,g,]))
          
        meank        <- moda(dictk[[i]][ty, g, ])
        
        # facciamo la media condizionata ai valori di k 
        meanTau <- mean(dictTau[[i]][ty, g, which(dictk[[i]][ty, g, ] == meank)])
        cellTy <- which(typeCell == ty)
        tyT0_off <- unique(typeCellT0_off[cellTy])
        
        # se ci sono degli infiniti li sostituiamo con numeri alti 
        T0_off_chain <- dictT0_off[[i]]
        lengthInf <- length(T0_off_chain[tyT0_off, g, which(dictT0_off[[i]][tyT0_off, g, ] == Inf)])
        if(lengthInf > dim(T0_off_chain)[3]/2){
          meanT0_off <- Inf
          meanT0_off_forPlot <- 500
        }else{
          T0_off_chain[tyT0_off, g, which(T0_off_chain[tyT0_off, g,] > 100)] <- NA
          meanT0_off <- mean(T0_off_chain[tyT0_off, g,], na.rm = TRUE)
          meanT0_off_forPlot <- meanT0_off
        }
        
        
        # CURVA CON ENTRAMBI GLI STEADY STATE
        if(SS_ON){
            # ramo on con medie
            u_plot <- u(seq(0, 500, 0.1), t0_off = Inf, t0_on = 0, u0_off = meanAlphaON/meanBeta, u0_on = meanAlphaOFF/meanBeta, k = rep(2, length(seq(0, 500, 0.1))), c(meanAlphaOFF, meanAlphaON), meanBeta)
            s_plot <- s(seq(0, 500, 0.1), t0_off = Inf, t0_on = 0, u0_off = meanAlphaON/meanBeta, u0_on = meanAlphaOFF/meanBeta, s0_off = meanAlphaON/meanGamma, s0_on = meanAlphaOFF/meanGamma, k =   rep(2, length(seq(0, 500, 0.1))), c(meanAlphaOFF, meanAlphaON), meanBeta, meanGamma)
            plot(s_plot, u_plot, type = "l", xlab = "s (spliced)", ylab = "u (unspliced)", main = paste("Gene", g, "ty", ty, "meanT0_off", round(meanT0_off,3)), lty = 2, lwd = 2, ylim = c(meanAlphaOFF/meanBeta, meanAlphaON/meanBeta), xlim = c(meanAlphaOFF/meanGamma, meanAlphaON/meanGamma), col = "green")
      
            # ramo on con dati reali
            u_plot <- u(seq(0, 500, 0.1), t0_off = Inf, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], k = rep(2, length(seq(0, 500, 0.1))), alpha_real[g,], beta_real[g])
            s_plot <- s(seq(0, 500, 0.1), t0_off = Inf, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], s0_off = alpha_real[g,2]/gamma_real[g], s0_on = alpha_real[g,1]/gamma_real[g], k =   rep(2, length(seq(0, 500, 0.1))), alpha_real[g,], beta_real[g], gamma_real[g])
            lines(s_plot, u_plot, type = "l",lty = 2, col = "red")

            # ramo off
            u_plot <- u(seq(0, 500, 0.1), t0_off = 0, t0_on = 0, u0_off = meanAlphaON/meanBeta, u0_on = meanAlphaOFF/meanBeta, k = rep(0, length(seq(0, 500, 0.1))), c(meanAlphaOFF, meanAlphaON), meanBeta)
            s_plot <- s(seq(0, 500, 0.1), t0_off = 0, t0_on = 0, u0_off = meanAlphaON/meanBeta, u0_on = meanAlphaOFF/meanBeta, s0_off = meanAlphaON/meanGamma, s0_on = meanAlphaOFF/meanGamma, k = rep(0, length(seq(0, 500, 0.1))), c(meanAlphaOFF, meanAlphaON), meanBeta, meanGamma)
            lines(s_plot, u_plot, lty = 2, lwd = 2, col = "green")

            # ramo off con dati  reali
            u_plot <- u(seq(0, 500, 0.1), t0_off = 0, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], k = rep(0, length(seq(0, 500, 0.1))),alpha_real[g,], beta_real[g])
            s_plot <- s(seq(0, 500, 0.1), t0_off = 0, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], s0_off = alpha_real[g,2]/gamma_real[g], s0_on = alpha_real[g,1]/gamma_real[g], k = rep(0, length(seq(0, 500, 0.1))), alpha_real[g,], beta_real[g], gamma_real[g])
            lines(s_plot, u_plot, lty = 2, col = "red")


        }
        


        # curve con SS che terminano prima 
        u0_off <- u0(t0_off = meanT0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(meanAlphaOFF, meanAlphaON), meanBeta)
        s0_off <- s0(t0_off = meanT0_off, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(meanAlphaOFF, meanAlphaON), meanBeta, meanGamma)

        u0_off_real <- u0(t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = alpha_real[g,], beta_real[g])
        s0_off_real <- s0(t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = alpha_real[g,], beta_real[g], gamma_real[g])

        u_plotON <- u(seq(0, meanT0_off_forPlot, 0.1), t0_off = meanT0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = rep(2, length(seq(0,  meanT0_off_forPlot, 0.1))), c(meanAlphaOFF, meanAlphaON), meanBeta)
        s_plotON <- s(seq(0, meanT0_off_forPlot, 0.1), t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = rep(2, length(seq(0,  meanT0_off_forPlot, 0.1))), alpha = c(meanAlphaOFF, meanAlphaON), beta = meanBeta, gamma = meanGamma)

        u_plotON_REAL <- u(seq(0, t0_off_real[tyT0_off,g], 0.1), t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, k = rep(2, length(seq(0,  t0_off_real[tyT0_off,g], 0.1))), alpha_real[g,], beta_real[g])

        s_plotON_REAL <- s(seq(0, t0_off_real[tyT0_off,g], 0.1), t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = rep(2, length(seq(0,  t0_off_real[tyT0_off,g], 0.1))), alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        u_plotOFF <- u_withTau(seq(0, 500, 0.1), u0_off = u0_off, u0_on = NA, k = rep(0, length(seq(0, 500, 0.1))), alpha = c(meanAlphaOFF, meanAlphaON), beta = meanBeta)
        s_plotOFF <- s_withTau(seq(0, 500, 0.1), u0_off = u0_off, u0_on = NA, s0_off = s0_off, s0_on = NA, k = rep(0, length(seq(0, 500, 0.1))), alpha = c(meanAlphaOFF, meanAlphaON), beta = meanBeta, gamma = meanGamma)

        u_plotOFF_REAL <- u_withTau(seq(0, 500, 0.1), u0_off = u0_off_real, u0_on = NA, k = rep(0, length(seq(0, 500, 0.1))), alpha =alpha_real[g,], beta = beta_real[g])
        s_plotOFF_REAL <- s_withTau(seq(0, 500, 0.1), u0_off = u0_off_real, u0_on = NA, s0_off = s0_off_real, s0_on = NA, k = rep(0, length(seq(0, 500, 0.1))), alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        # ramo on simulato
        if(SS_ON){
          lines(s_plotON, u_plotON , type = "l", col = "green", lwd = 1) # , xlab = "u", ylab = "s", main = paste("MEAN Gene ", g, nameGenes[g], " ty ", ty, xlim = c(meanAlphaOFF/meanBeta, meanAlphaON/meanBeta), ylim = c(meanAlphaOFF/meanGamma, meanAlphaON/meanGamma)))?
        }else{
          plot(s_plotON, u_plotON, type = "l", col = "green", lwd = 1 , ylab = "u", xlab = "s", main = paste("Gene", g, " ty ", ty, "meanT0_off", round(meanT0_off, 3)), ylim = c(min(c(u_plotOFF, u_plotOFF_REAL, u_plotON, u_plotON_REAL)), max(c(u_plotOFF, u_plotOFF_REAL, u_plotON, u_plotON_REAL))), xlim = c(min(c(s_plotOFF, s_plotOFF_REAL, s_plotON, s_plotON_REAL)), max(c(s_plotOFF, s_plotOFF_REAL, s_plotON, s_plotON_REAL))))
        }
        
        
        # ramo on reale
        lines(s_plotON_REAL, u_plotON_REAL, type = "l", col = "red", lwd = 1)

        # ramo off simulato
        lines(s_plotOFF, u_plotOFF, type = "l", col = "green", lwd = 1)

        # ramo off reale
        lines(s_plotOFF_REAL, u_plotOFF_REAL, type = "l", col = "red", lwd = 1)
        
    

        u_point <- u_withTau(meanTau, u0_off = u0_off, u0_on = meanAlphaOFF/meanBeta, k = meank, alpha = c(meanAlphaOFF, meanAlphaON), beta = meanBeta)
        s_point <- s_withTau(meanTau, u0_off = u0_off, u0_on = meanAlphaOFF/meanBeta, s0_off = s0_off, s0_on = meanAlphaOFF/meanGamma, k = meank, alpha = c(meanAlphaOFF, meanAlphaON), beta = meanBeta, gamma = meanGamma)
#         print(c(ty,g, meanT0_off, t0_off_real[tyT0_off,g], u_point, s_point))

        points(s_point, u_point, col = "green", pch = 19, cex = 2)
      
        points(s_withTau(tau_real[cellTy[1],g], u0_off = u0_off_real, u0_on = alpha_real[g,1]/beta_real[g], s0_off = s0_off_real, s0_on = alpha_real[g,1]/gamma_real[g], k = k_real[cellTy[1], g], alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g]), u_withTau(tau_real[cellTy[1],g], u0_off = u0_off_real, u0_on = alpha_real[g,1]/beta_real[g], k = k_real[cellTy[1],g], alpha = alpha_real[g,], beta = beta_real[g]), col = "red", pch = 18, cex = 2)


        # legend("topleft", legend = c("Real curve", "Estimated curve", "Real tau", "Estimated tau"),  lwd=2.5, lty=c(1,1,0,0), pch=c(NA,NA,16,8), col = c("black", "red", "black", "black"))
      }
}


phasePlotMEDIAN <- function(dictLogSS, dictTime, dictk, typeC, nameGenes, geni){
    for(g in geni){
        phasePlotMEDIAN_gene(dictLogSS, dictTime, dictk, typeC, nameGenes, g)
    }
}

# phasePlotMEDIAN_gene <- function(dictLogSS, dictTime, dictk, typeC, nameGenes, g){  
#     pchPoss = seq(15,(15 + length(dictLogSS))) 
    
#     for(i in 1:length(dictLogSS)){
#         medianAlphaOFF <- median(exp(dictLogSS[[i]][1,g,])*exp(dictLogSS[[i]][4,g,]))
#         medianAlphaON  <- median(exp(dictLogSS[[i]][4,g,])*(exp(dictLogSS[[i]][3,g,]) + exp(dictLogSS[[i]][1,g,])))
#         medianBeta     <- median(exp(dictLogSS[[i]][4,g,]))
#         medianGamma    <- median(exp(dictLogSS[[i]][4,g,])*exp(dictLogSS[[i]][1,g,])/exp(dictLogSS[[i]][2,g,]))
    
#         mediank        <- apply(dictk[[i]][, g, ], MARGIN = 1, moda)
#         tyOFF        <- which(mediank == 0)
#         tyON         <- which(mediank == 2)
#         medianTime     <- rep(NA, length(mediank))
#         # facciamo la media condizionata ai valori di k 
#         for(ty in 1:length(mediank)){
#             if(mediank[ty] == 0){
#                 tmp <- 0
#             }else{
#                 tmp <- 2
#             }
#             medianTime[ty] <- median(dictTime[[i]][ty, g, which(dictk[[i]][ty, g, ] == tmp) ])
#         }
         
#         if(i == 1){
#             plot(u(seq(0, 5000, 0.1), k = rep(2, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta), s(seq(0, 5000, 0.1), k = rep(2, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma), type = "l", col = i, lwd = 1, xlab = "u", ylab = "s", main = paste("MEDIAN Gene ", g, nameGenes[g]))
#         }else{
#             lines(u(seq(0, 5000, 0.1), k = rep(2, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta), s(seq(0, 5000, 0.1), k = rep(2, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma), type = "l", col = i, lwd = 1, xlab = "u", ylab = "s")
#         }
#         lines(u(seq(0, 5000, 0.1), k = rep(0, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta), s(seq(0, 5000, 0.1), k = rep(0, length(seq(0, 5000, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma), type = "l", col = i, lwd = 1)
        
#         points(u(medianTime, mediank, c(medianAlphaOFF, medianAlphaON), medianBeta), s(medianTime, mediank, c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma), col = rep(1:length(unique(typeC))), pch = pchPoss[i], cex = 1.5)
      
     
#         # legend("topleft", legend = c("Real curve", "Estimated curve", "Real tau", "Estimated tau"),  lwd=2.5, lty=c(1,1,0,0), pch=c(NA,NA,16,8), col = c("black", "red", "black", "black"))
#       }
# }

phasePlotMEDIAN_gene_ty <- function(dictLogSS, dictT0_off, dictTau, dictk, alpha_real, beta_real, gamma_real, tau_real, t0_off_real, k_real, nameGenes, g, ty,  typeCell, typeCellT0_off, SS_ON){  
    pchPoss = seq(15,(15 + length(dictLogSS))) 
    
    for(i in 1:length(dictLogSS)){
        medianAlphaOFF <- median(exp(dictLogSS[[i]][1,g,])*exp(dictLogSS[[i]][4,g,]))
        medianAlphaON  <- median(exp(dictLogSS[[i]][4,g,])*(exp(dictLogSS[[i]][3,g,]) + exp(dictLogSS[[i]][1,g,])))
        medianBeta     <- median(exp(dictLogSS[[i]][4,g,]))
        medianGamma    <- median(exp(dictLogSS[[i]][4,g,])*exp(dictLogSS[[i]][1,g,])/exp(dictLogSS[[i]][2,g,]))
          
        mediank        <- moda(dictk[[i]][ty, g, ])
        
        # facciamo la media condizionata ai valori di k 
        medianTau <- median(dictTau[[i]][ty, g, which(dictk[[i]][ty, g, ] == mediank)])
        cellTy <- which(typeCell == ty)
        tyT0_off <- unique(typeCellT0_off[cellTy])
        
        # se ci sono degli infiniti li sostituiamo con numeri alti 
        T0_off_chain <- dictT0_off[[i]]
        lengthInf <- length(T0_off_chain[tyT0_off, g, which(dictT0_off[[i]][tyT0_off, g, ] == Inf)])
        if(lengthInf > dim(T0_off_chain)[3]/2){
          medianT0_off <- Inf
          medianT0_off_forPlot <- 100
        }else{
          T0_off_chain[tyT0_off, g, which(T0_off_chain[tyT0_off, g,] > 100)] <- NA
          medianT0_off <- median(T0_off_chain[tyT0_off, g,], na.rm = TRUE)
          medianT0_off_forPlot <- medianT0_off
        }
        
        
        # CURVA CON ENTRAMBI GLI STEADY STATE
        if(SS_ON){
            # ramo on con medie
            u_plot <- u(seq(0, 100, 0.1), t0_off = Inf, t0_on = 0, u0_off = medianAlphaON/medianBeta, u0_on = medianAlphaOFF/medianBeta, k = rep(2, length(seq(0, 100, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta)
            s_plot <- s(seq(0, 100, 0.1), t0_off = Inf, t0_on = 0, u0_off = medianAlphaON/medianBeta, u0_on = medianAlphaOFF/medianBeta, s0_off = medianAlphaON/medianGamma, s0_on = medianAlphaOFF/medianGamma, k =   rep(2, length(seq(0, 100, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma)
            plot(s_plot, u_plot, type = "l", xlab = "s (spliced)", ylab = "u (unspliced)", main = paste("Gene", g, "ty", ty, "medianT0_off", round(medianT0_off,3)), lty = 2, lwd = 2, ylim = c(medianAlphaOFF/medianBeta, medianAlphaON/medianBeta), xlim = c(medianAlphaOFF/medianGamma, medianAlphaON/medianGamma), col = "green")

            # ramo on con dati reali
            u_plot <- u(seq(0, 100, 0.1), t0_off = Inf, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], k = rep(2, length(seq(0, 100, 0.1))), alpha_real[g,], beta_real[g])
            s_plot <- s(seq(0, 100, 0.1), t0_off = Inf, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], s0_off = alpha_real[g,2]/gamma_real[g], s0_on = alpha_real[g,1]/gamma_real[g], k =   rep(2, length(seq(0, 100, 0.1))), alpha_real[g,], beta_real[g], gamma_real[g])

            lines(s_plot, u_plot, type = "l",lty = 2, col = "red")

            # ramo off
            u_plot <- u(seq(0, 100, 0.1), t0_off = 0, t0_on = 0, u0_off = medianAlphaON/medianBeta, u0_on = medianAlphaOFF/medianBeta, k = rep(0, length(seq(0, 100, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta)
            s_plot <- s(seq(0, 100, 0.1), t0_off = 0, t0_on = 0, u0_off = medianAlphaON/medianBeta, u0_on = medianAlphaOFF/medianBeta, s0_off = medianAlphaON/medianGamma, s0_on = medianAlphaOFF/medianGamma, k = rep(0, length(seq(0, 100, 0.1))), c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma)
            lines(s_plot, u_plot, lty = 2, lwd = 2, col = "green")

            # ramo off con dati  reali
            u_plot <- u(seq(0, 100, 0.1), t0_off = 0, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], k = rep(0, length(seq(0, 100, 0.1))),alpha_real[g,], beta_real[g])
            s_plot <- s(seq(0, 100, 0.1), t0_off = 0, t0_on = 0, u0_off = alpha_real[g,2]/beta_real[g], u0_on = alpha_real[g,1]/beta_real[g], s0_off = alpha_real[g,2]/gamma_real[g], s0_on = alpha_real[g,1]/gamma_real[g], k = rep(0, length(seq(0, 100, 0.1))), alpha_real[g,], beta_real[g], gamma_real[g])
            lines(s_plot, u_plot, lty = 2, col = "red")


        }
        


        # curve con SS che terminano prima 
        u0_off <- u0(t0_off = medianT0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(medianAlphaOFF, medianAlphaON), medianBeta)
        s0_off <- s0(t0_off = medianT0_off, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(medianAlphaOFF, medianAlphaON), medianBeta, medianGamma)

        u0_off_real <- u0(t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = alpha_real[g,], beta_real[g])
        s0_off_real <- s0(t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = alpha_real[g,], beta_real[g], gamma_real[g])

        u_plotON <- u(seq(0, medianT0_off_forPlot, 0.005), t0_off = medianT0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = rep(2, length(seq(0,  medianT0_off_forPlot, 0.005))), c(medianAlphaOFF, medianAlphaON), medianBeta)
        s_plotON <- s(seq(0, medianT0_off_forPlot, 0.005), t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = rep(2, length(seq(0,  medianT0_off_forPlot, 0.005))), alpha = c(medianAlphaOFF, medianAlphaON), beta = medianBeta, gamma = medianGamma)

        u_plotON_REAL <- u(seq(0, t0_off_real[tyT0_off,g], 0.005), t0_off = t0_off_real[tyT0_off,g], t0_on = 0, u0_off = NA, u0_on = NA, k = rep(2, length(seq(0,  t0_off_real[tyT0_off,g], 0.005))), alpha_real[g,], beta_real[g])

        s_plotON_REAL <- s(seq(0, t0_off_real[tyT0_off,g], 0.005), t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = rep(2, length(seq(0,  t0_off_real[tyT0_off,g], 0.005))), alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        u_plotOFF <- u_withTau(seq(0, 100, 0.005), u0_off = u0_off, u0_on = NA, k = rep(0, length(seq(0, 100, 0.005))), alpha = c(medianAlphaOFF, medianAlphaON), beta = medianBeta)
        s_plotOFF <- s_withTau(seq(0, 100, 0.005), u0_off = u0_off, u0_on = NA, s0_off = s0_off, s0_on = NA, k = rep(0, length(seq(0, 100, 0.005))), alpha = c(medianAlphaOFF, medianAlphaON), beta = medianBeta, gamma = medianGamma)

        u_plotOFF_REAL <- u_withTau(seq(0, 100, 0.005), u0_off = u0_off_real, u0_on = NA, k = rep(0, length(seq(0, 100, 0.005))), alpha =alpha_real[g,], beta = beta_real[g])
        s_plotOFF_REAL <- s_withTau(seq(0, 100, 0.005), u0_off = u0_off_real, u0_on = NA, s0_off = s0_off_real, s0_on = NA, k = rep(0, length(seq(0, 100, 0.005))), alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g])

        # ramo on simulato
        if(SS_ON){
          lines(s_plotON, u_plotON , type = "l", col = "green", lwd = 1) # , xlab = "u", ylab = "s", main = paste("MEAN Gene ", g, nameGenes[g], " ty ", ty, xlim = c(meanAlphaOFF/meanBeta, meanAlphaON/meanBeta), ylim = c(meanAlphaOFF/meanGamma, meanAlphaON/meanGamma)))?
        }else{
          plot(s_plotON, u_plotON, type = "l", col = "green", lwd = 1 , ylab = "u", xlab = "s", main = paste("Gene", g, " ty ", ty, "medianT0_off", round(medianT0_off, 3)), ylim = c(min(c(u_plotOFF, u_plotOFF_REAL, u_plotON, u_plotON_REAL)), max(c(u_plotOFF, u_plotOFF_REAL, u_plotON, u_plotON_REAL))), xlim = c(min(c(s_plotOFF, s_plotOFF_REAL, s_plotON, s_plotON_REAL)), max(c(s_plotOFF, s_plotOFF_REAL, s_plotON, s_plotON_REAL))))
        }
        
        # ramo on reale
        lines(s_plotON_REAL, u_plotON_REAL, type = "l", col = "red", lwd = 1)

        # ramo off simulato
        lines(s_plotOFF, u_plotOFF, type = "l", col = "green", lwd = 1)

        # ramo off reale
        lines(s_plotOFF_REAL, u_plotOFF_REAL, type = "l", col = "red", lwd = 1)
        
    
        print(c(max(s_plotON_REAL), max(u_plotON_REAL)))
        print(c(s_plotOFF_REAL[1], u_plotOFF_REAL[1]))

        u_point <- u_withTau(medianTau, u0_off = u0_off, u0_on = medianAlphaOFF/medianBeta, k = mediank, alpha = c(medianAlphaOFF, medianAlphaON), beta = medianBeta)
        s_point <- s_withTau(medianTau, u0_off = u0_off, u0_on = medianAlphaOFF/medianBeta, s0_off = s0_off, s0_on = medianAlphaOFF/medianGamma, k = mediank, alpha = c(medianAlphaOFF, medianAlphaON), beta = medianBeta, gamma = medianGamma)
#         print(c(ty,g, meanT0_off, t0_off_real[tyT0_off,g], u_point, s_point))

        points(s_point, u_point, col = "green", pch = 19, cex = 2)
      
        points(s_withTau(tau_real[cellTy[1],g], u0_off = u0_off_real, u0_on = alpha_real[g,1]/beta_real[g], s0_off = s0_off_real, s0_on = alpha_real[g,1]/gamma_real[g], k = k_real[cellTy[1], g], alpha = alpha_real[g,], beta = beta_real[g], gamma = gamma_real[g]), u_withTau(tau_real[cellTy[1],g], u0_off = u0_off_real, u0_on = alpha_real[g,1]/beta_real[g], k = k_real[cellTy[1],g], alpha = alpha_real[g,], beta = beta_real[g]), col = "red", pch = 18, cex = 2)


        # legend("topleft", legend = c("Real curve", "Estimated curve", "Real tau", "Estimated tau"),  lwd=2.5, lty=c(1,1,0,0), pch=c(NA,NA,16,8), col = c("black", "red", "black", "black"))
      }
}


############################ SCATTER PLOT ################################
# per ogni campione a posteriori: 
# prendo 10 punti da uOFF e uON e plotto la curva

scatterU <- function(dictLogSS, nameGenes, geni, n){
    for(g in geni){
        scatterU_gene(dictLogSS, nameGenes, 0, g, n)
        scatterU_gene(dictLogSS, nameGenes, 2, g, n)
    }
}

scatterU_gene <- function(dictLogSS, nameGenes, k, g, n){
   
    maxUON <- 0
    maxSON <- 0
    for(i in 1:length(dictLogSS)){
        maxUON <- max(maxUON, exp(dictLogSS[[i]][1,g,]) + exp(dictLogSS[[i]][3,g,]))
        maxSON <- max(maxSON, (exp(dictLogSS[[i]][1,g,]) + exp(dictLogSS[[i]][3,g,]))*exp(dictLogSS[[i]][2,g,])/exp(dictLogSS[[i]][1,g,]))
    }
    plot(NA, xlim=c(0,maxUON), ylim=c(0,maxSON), bty='n', xlab='unspliced', ylab='spliced', main = paste("Scatter U gene ", g, nameGenes[g], "phase", k))
    pchPoss <- c(20, 18)
    for(i in 1:length(dictLogSS)){
        col <-   i # rainbow(n+2) #
        uOFF <- exp(dictLogSS[[i]][1,g,])
        uON  <- exp(dictLogSS[[i]][1,g,]) + exp(dictLogSS[[i]][3,g,])
        alphaOFF <- uOFF*exp(dictLogSS[[i]][4,g,])
        alphaON  <- uON *exp(dictLogSS[[i]][4,g,]) 
        beta     <- exp(dictLogSS[[i]][4,g,])
        gamma    <- alphaON/exp(dictLogSS[[i]][2,g,])
        
        
        for(iter in 1:length(uOFF)){
            breaks <- seq(uOFF[iter], uON[iter], length.out = (n+2))           
            for(ty in 1:length(dictk[[i]][, g, iter])){           
                tau <- tau_inv(breaks, rep(k, n + 2 ), c(alphaOFF[iter], alphaON[iter]), beta[iter], gamma[iter])
                points(breaks, s(tau, rep(k, n + 2 ),  c(alphaOFF[iter], alphaON[iter]), beta[iter], gamma[iter]), col =  alpha(col, 0.4) , pch = pchPoss[[i]], cex = 0.5)
                
            }
        }
    }
}


# per ogni campione a posteriori: 
# prendo 10 punti da sOFF e sON e plotto la curva
scatterS <- function(dictLogSS, nameGenes, geni, n){
    for(g in geni){
        scatterS_gene(dictLogSS, nameGenes, 0, g, n)
        scatterS_gene(dictLogSS, nameGenes, 2, g, n)
    }
}

scatterS_gene <- function(dictLogSS, nameGenes, k, g, n){      
    maxUON <- 0
    maxSON <- 0
    for(i in 1:length(dictLogSS)){
        maxUON <- max(maxUON, exp(dictLogSS[[i]][1,g,]) + exp(dictLogSS[[i]][3,g,]))
        maxSON <- max(maxSON, (exp(dictLogSS[[i]][1,g,]) + exp(dictLogSS[[i]][3,g,]))*exp(dictLogSS[[i]][2,g,])/exp(dictLogSS[[i]][1,g,]))
    }
    plot(NA, xlim=c(0,maxUON), ylim=c(0,maxSON), bty='n', xlab='unspliced', ylab='spliced', main = paste("Scatter S gene ", g, nameGenes[g], "phase", k))
    pchPoss <- c(20, 18)
    for(i in 1:length(dictLogSS)){
        col <-  i # rainbow(n+2) 
        uOFF <- exp(dictLogSS[[i]][1,g,])
        beta     <- exp(dictLogSS[[i]][4,g,])
        alphaOFF <- uOFF * beta
        gamma <- alphaOFF/exp(dictLogSS[[i]][2, g,])
        alphaON <- exp(dictLogSS[[i]][3, g,])*beta + alphaOFF

        sOFF <- exp(dictLogSS[[i]][2, g,])
        sON  <- alphaON/gamma
        
        
        for(iter in 1:length(sOFF)){
            breaks <- seq(sOFF[iter], sON[iter], length.out = (n+2))           
            for(ty in 1:(length(dictk[[i]][, g, iter]))){           
                u_breaks <- u_di_s(breaks[2:(n+1)], rep(k, n  ), c(alphaOFF[iter], alphaON[iter]), beta[iter], gamma[iter])
                points(c(alphaOFF[iter]/beta[iter], u_breaks, alphaON[iter]/beta[iter]), breaks, col = alpha(col, 0.4), pch = pchPoss[i])
            }
        }
    }
}



scatterTime <- function(dictTime, dictk, nameTypeCells){
    for(i in 1:length(dictk)){
        meank <- apply(dictk[[i]], MARGIN = c(1, 2), moda)
        # facciamo la media condizionata ai valori di k 
        tmp <- ifelse(meank == 0, 0, 2)
        tmpIter <- rep(tmp, dim(dictk[[i]])[3])
        attr(tmpIter, "dim") <- dim(dictk[[i]])
        
        filter <- ifelse(dictk[[i]] == tmpIter, 1, 0)
        meanTime <- apply(dictTime[[i]]*filter, MARGIN = c(1, 2), sum)/apply(filter, MARGIN = c(1, 2), sum)
        
                
        for(k in 1:(n_typeC - 1)){
            for(j in (i+1):n_typeC){
                plot(meanTime[k, ], col = k, cex = 0.5, pch = 20, main = paste( "Mean Time ", nameTypeCells[[as.character(k)]], "VS", nameTypeCells[[as.character(j)]]))
                points(meanTime[j, ], col = j, cex = 0.5, pch = 20)
            }
        }
        
    }    
}

scatterTime_geneSpecific <- function(dictTime, dictk, nameTypeCells, nameGenesPancreas, indexGenes){
    n_typeC <- dim(dictTime[[1]])[1]
    n_iter <- dim(dictTime[[1]])[3]
    for(g in indexGenes){
        for(i in 1:length(dictk)){
            vecTime <- as.vector(dictTime[[i]][, g, ])
            k0 <- which(dictk[[i]][,g,] == 0)
            k2 <- which(dictk[[i]][,g,] == 2)
            if(length(k0) > 0){
                plot(vecTime[k0], col = rep(rainbow(n_typeC), n_iter)[k0], cex = 0.7, pch = 20, main = paste( "ph 0, g", g, nameGenesPancreas[g]))
            }else{
                plot(NA, xlim=c(0,3), ylim=c(0,3), bty='n',xaxt='n', yaxt='n', xlab='', ylab='', main = paste( "ph 0, g", g, nameGenesPancreas[g]))
                text(0,2.5,paste("No observations"), pos = 4)      
            }
            
            if(length(k2) > 0){
                plot(vecTime[k2], col = rep(rainbow(n_typeC), n_iter)[k2], cex = 0.7, pch = 20, main = paste("ph 2 g", g, nameGenesPancreas[g]))
            }else{
                plot(NA, xlim=c(0,3), ylim=c(0,3), bty='n',xaxt='n', yaxt='n', xlab='', ylab='', main = paste( "ph 2, g", g, nameGenesPancreas[g]))
                text(0,2.5,paste("No observations"), pos = 4) 
            }   
        }    
    }
}

contingencyTable <- function(dictTime, dictk, n_typeC, nameTypeCells){
    for(i in 1:length(dictk)){
        contingengyTable_funct(dictTime, dictk, "mean", i, n_typeC, nameTypeCells)    
        contingengyTable_funct(dictTime, dictk, "median", i, n_typeC, nameTypeCells)
    }
}

contingengyTable_funct <- function(dictTime, dictk, type, i, n_typeC, nameTypeCells){
    if(type == "mean"){
        plot(NA, xlim = c(0, 1), ylim = c(0,1), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste("MEAN"))
        meank <- apply(dictk[[i]], MARGIN = c(1,2), moda)
        meankIter <- rep(meank, dim(dictk[[i]])[3])
        filter <- ifelse(dictk[[i]] == meankIter, 1, 0)
        meanTime <- rep(NA, dim(dictTime[[i]])[1]*dim(dictTime[[i]])[2])
        attr(meanTime, "dim") <- dim(dictTime[[i]])[1:2]
        for(g in 1:dim(dictTime[[i]])[2]){
            for(ty in 1:dim(dictTime[[i]])[1]){
                meanTime[ty,g] <- mean(dictTime[[i]][ty, g, which(filter[ty, g,] == 1)])
            }
            
        }
        # meanTime <- apply(dictTime[[i]]*filter, MARGIN = c(1, 2), sum)/apply(filter, MARGIN = c(1, 2), sum)
        
        par(mfrow = c(3, 1))
        for(k in 1:(n_typeC-1)){
            for(j in (k + 1):n_typeC){
                # facciamo tabella di contingenza in cui guardiamo in quali tipi di cellula il tempo e maggiore in un tipo e minore nell'altro e si trovano sullo stesso ramo 
                # Confrontiamo anche in quanti casi si trovano su stati diversi 
                c1 <- c(length(which(meank[k,] == 0 & meank[j,] == 0)), length(which(meank[k,] == 0 & meank[j,] == 2)))
                c2 <- c(length(which(meank[k,] == 2 & meank[j,] == 0)), length(which(meank[k,] == 2 & meank[j,] == 2)))
                df_k <- data.frame("j_1 = 0" = c1, "j_1 = 2" = c2)
                rownames(df_k) <- c("j_2 = 0", "j_2 = 2")
    
                
                #             entrambe k = 0 | k = 2
                # t_j > t_k 
                # t_j < t_k
                c1 <- c(length(which(meanTime[k,] >= meanTime[j,] & (meank[k,] == 0 & meank[j,] == 0))), length(which(meanTime[k,] < meanTime[j,] & (meank[k,] == 0 & meank[j,] == 0))))
                c2 <- c(length(which(meanTime[k,] >= meanTime[j,] & (meank[k,] == 2 & meank[j,] == 2))), length(which(meanTime[k,] < meanTime[j,] & (meank[k,] == 2 & meank[j,] == 2))))
    
                df_time <- data.frame("k=0" = c1, "k = 2" = c2)
                rownames(df_time) <- c("t_{j_1} > t_{j_2}", "t_{j_1} < t_{j_2}")
    
                grid.arrange(tableGrob(df_k), tableGrob(df_time), ncol = 1, nrow = 3, top = textGrob(paste("j_1 = ", k, nameTypeCells[[as.character(k)]], ", j_2 = ", j, nameTypeCells[[as.character(j)]]),gp=gpar(fontsize=20,font=3)))
    
    
            }
        }
    }else if(type == "median"){
        plot(NA, xlim = c(0, 1), ylim = c(0,1), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste("MEDIAN"))
        medianTime <- rep(NA, dim(dictTime[[i]])[1]*dim(dictTime[[i]])[2])
        attr(medianTime, "dim") <- dim(dictTime[[i]])[1:2]
        for(g in 1:dim(dictTime[[i]])[2]){
            for(ty in 1:dim(dictTime[[i]])[1]){
                medianTime[ty,g] <- median(dictTime[[i]][ty, g, which(filter[ty, g,] == 1)])
            }
            
        }
    
        for(k in 1:(n_typeC-1)){
            for(j in (k + 1):n_typeC){
                # facciamo tabella di contingenza in cui guardiamo in quali tipi di cellula il tempo e maggiore in un tipo e minore nell'altro e si trovano sullo stesso ramo 
                # Confrontiamo anche in quanti casi si trovano su stati diversi 
                c1 <- c(length(which(meank[k,] == 0 & meank[j,] == 0)), length(which(meank[k,] == 0 & meank[j,] == 2)))
                c2 <- c(length(which(meank[k,] == 2 & meank[j,] == 0)), length(which(meank[k,] == 2 & meank[j,] == 2)))
                df_k <- data.frame("j_1 = 0" = c1, "j_1 = 2" = c2)
                rownames(df_k) <- c("j_2 = 0", "j_2 = 2")
                #             entrambe k = 0 | k = 2
                # t_j > t_k 
                # t_j < t_k
                c1 <- c(length(which(medianTime[k,] >= medianTime[j,] & (meank[k,] == 0 & meank[j,] == 0))), length(which(medianTime[k,] < medianTime[j,] & (meank[k,] == 0 & meank[j,] == 0))))
                c2 <- c(length(which(medianTime[k,] >= medianTime[j,] & (meank[k,] == 2 & meank[j,] == 2))), length(which(medianTime[k,] < medianTime[j,] & (meank[k,] == 2 & meank[j,] == 2))))
    
                df_time <- data.frame("k=0" = c1, "k = 2" = c2)
                rownames(df_time) <- c("t_{j_1} > t_{j_2}", "t_{j_1} < t_{j_2}")
    
                grid.arrange(tableGrob(df_k), tableGrob(df_time), ncol = 1, nrow = 3, top = textGrob(paste("j_1 = ", k, nameTypeCells[[as.character(k)]], ", j_2 = ", j, nameTypeCells[[as.character(j)]]),gp=gpar(fontsize=20,font=3)))
            }
        }
    }    
}


reshapeU_MCMC <- function(dictLogSS, dictTime, dictk, n_typeC, nameGenes){
    out <- list()
    for(i in 1:length(dictTime)){
        u_MCMC <- rep(NA, prod(dim(dictTime[[i]])))
        attr(u_MCMC, "dim") <- dim(dictTime[[i]])
        
        for(g in 1:dim(dictTime[[i]])[2]){
            for(iter in 1:dim(dictTime[[i]])[3]){
                SS <- exp(dictLogSS[[i]][,g,iter])
                rates <- ss_to_rates(SS[4], SS[1], SS[2], SS[3])
                alpha_off <- rates[1]
                alpha_on <- rates[2]
                gamma <- rates[3]
                beta <- SS[4]
                u_MCMC[,g,iter] <- u(dictTime[[i]][,g,iter], dictk[[i]][,g,iter], c(alpha_off, alpha_on), beta)
            }
        }

        name <- names(dictTime)[[i]]
        tmp <- u_MCMC
        attr(tmp, "dim") <- c(n_typeC*dim(u_MCMC)[3], dim(u_MCMC)[2])
        for(iter in 1:dim(u_MCMC)[3]){
            tmp[((iter-1)*n_typeC + 1):(iter*n_typeC),] <- u_MCMC[,,iter]
        }
        out[[name]] <- data.frame(tmp)
        colnames(out[[name]]) <- nameGenesPancreas
        out[[name]]$type <- rep(seq(1:n_typeC), dim(u_MCMC)[3])
    }
    return(out)
}


reshapeS_MCMC <- function(dictLogSS, dictTime, dictk, n_typeC, nameGenes){
    out <- list()
    for(i in 1:length(dictTime)){
        s_MCMC <- rep(NA, prod(dim(dictTime[[i]])))
        attr(s_MCMC, "dim") <- dim(dictTime[[i]])
        
        for(g in 1:dim(dictTime[[i]])[2]){
            for(iter in 1:dim(dictTime[[i]])[3]){
                SS <- exp(dictLogSS[[i]][,g,iter])
                rates <- ss_to_rates(SS[4], SS[1], SS[2], SS[3])
                alpha_off <- rates[1]
                alpha_on <- rates[2]
                gamma <- rates[3]
                beta <- SS[4]
                s_MCMC[,g,iter] <- s(dictTime[[i]][,g,iter], dictk[[i]][,g,iter], c(alpha_off, alpha_on), beta, gamma)
            }
        }

        name <- names(dictTime)[[i]]
        tmp <- s_MCMC
        attr(tmp, "dim") <- c(n_typeC*dim(s_MCMC)[3], dim(s_MCMC)[2])
        for(iter in 1:dim(s_MCMC)[3]){
            tmp[((iter-1)*n_typeC + 1):(iter*n_typeC),] <- s_MCMC[,,iter]
        }
        out[[name]] <- data.frame(tmp)
        colnames(out[[name]]) <- nameGenesPancreas
        out[[name]]$type <- rep(seq(1:n_typeC), dim(u_MCMC)[3])
    }
    return(out)
}


createFileName <- function(simSS, simT0, simTau, simEta, simCatt){  
  nameFILE <- ""
  if(simSS){
    nameFILE <- paste(nameFILE, "-SS", sep = "")
  }
  if(simT0){
    nameFILE <- paste(nameFILE, "-T0", sep = "")
  }
  if(simTau){
    nameFILE <- paste(nameFILE, "-Tau", sep = "")
  }
  if(simEta){
    nameFILE <- paste(nameFILE, "-Eta", sep = "")
  }
  if(simCatt){
    nameFILE <- paste(nameFILE, "-Catt", sep = "")
  }
  

  if(startsWith(nameFILE ,"-")){
    nameFILE <- sub("-", "", nameFILE)
  }
  return(nameFILE)
}



loadMCMCdata <- function(simSS, simT0, simTau, simEta, simCatt, typeSIMreal, mcmcIter = NA, path = NA){
  nameFILE <- createFileName(simSS, simT0, simTau, simEta, simCatt)
  n_clusters <- getN_clusters(typeSIMreal)
  typeSIM <- paste("epsilon", n_clusters, sep = "")

  if(is.na(path)){
    error("Path to results need to be provided")
  }
  
  if(endsWith(typeSIMreal, "low")){
    path <- paste(path, "LOW counts/", typeSIMreal, "/", sep = "")
  }else{
    path <- paste(path, "HIGH counts/", typeSIMreal, "/", sep = "")
  }
  
  nomeSim <- getNomeSim(typeSIMreal)
  return(paste(path, nameFILE, "_", typeSIMreal, "_", nomeSim, "_",typeSIM, "_", mcmcIter, ".RData", sep = ""))

}

getNomeSim <- function(typeSIMreal){
  if(typeSIMreal == "A"){
    nomeSim <- "sim_3_TRUE_1"
  }else if(typeSIMreal == "B"){
    nomeSim <- "sim_3_FALSE_1"
  }else if(typeSIMreal == "D"){
    nomeSim <- "sim_3_TRUE_1"
  }else if(typeSIMreal == "E"){
    nomeSim <- "sim_3_FALSE_1"
  }else if(typeSIMreal == "F"){
    nomeSim <- "sim_3_TRUE_1"
  }else if(typeSIMreal == "G"){
    nomeSim <- "sim_3_FALSE_1"
  }else if(typeSIMreal == "H"){
    nomeSim <- "sim_3_TRUE_1"
  }else if(typeSIMreal == "I"){
    nomeSim <- "sim_3_FALSE_1"
  }else if(typeSIMreal == "L"){
    nomeSim <- "sim_3_TRUE_1_catt1Eta0"
  }else if(typeSIMreal == "A_low"){
    nomeSim <- "sim_3_TRUE_2"
  }else if(typeSIMreal == "B_low"){
    nomeSim <- "sim_3_FALSE_2"
  }else if(typeSIMreal == "D_low"){
    nomeSim <- "sim_3_TRUE_2"
  }else if(typeSIMreal == "E_low"){
    nomeSim <- "sim_3_FALSE_2"
  }else if(typeSIMreal == "F_low"){
    nomeSim <- "sim_3_TRUE_2"
  }else if(typeSIMreal == "G_low"){
    nomeSim <- "sim_3_FALSE_2"
  }else if(typeSIMreal == "H_low"){
    nomeSim <- "sim_3_TRUE_2"
  }else if(typeSIMreal == "I_low"){
    nomeSim <- "sim_3_FALSE_2"
  }else if(typeSIMreal == "L_low"){
    nomeSim <- "sim_3_TRUE_2_catt1Eta0"
  }
  return(nomeSim)
}

getN_clusters <- function(typeSIMreal){
  if(startsWith(typeSIMreal, "A") | startsWith(typeSIMreal, "B") | startsWith(typeSIMreal, "L")){
    nClusters <- 1
  }else if(startsWith(typeSIMreal, "D")| startsWith(typeSIMreal, "E")){
    nClusters <- 3
  }else if(startsWith(typeSIMreal, "F")| startsWith(typeSIMreal, "G")){
    nClusters <- 10
  }else if(startsWith(typeSIMreal, "H") | startsWith(typeSIMreal, "I")){
    nClusters <- "ALL"
  }

  return(nClusters)
}

deleteNA_MCMC <- function(chain){
  dimStart <- dim(chain)

  div <- TRUE
  if(is.null(dimStart)){
    div <- FALSE
    dimStart <- length(chain)
  }

  
  res <- chain[which(!is.na(chain), arr.ind = TRUE)]

  if(!div){
    denom <- 1
  }else{
    denom <- prod(dimStart[-length(dimStart)])
  }
  attr(res, "dim") <- c(dimStart[-length(dimStart)], length(res)/denom)

  return(res)
}


# library(textTinyR)
get_size<-function(adata, layer=NA){
  #   ---------  
  # Get counts per observation in a layer.
  #   Arguments
  #   ---------
  #   adata
  #       Annotated data matrix.
  #   layer
  #       Name of later for which to retrieve initial size.
  #   Returns
  #   -------
  #   np.ndarray
  #       Initial counts per observation in the specified layer.
  #   
  
  if(is.na(layer)){
    X <- adata$X     
  }  else{
    X <- adata$layers[layer]
  }
  row_sum <- rowSums(X)
  res <- data.frame(unname(row_sum), row.names = names(row_sum))
  colnames(res) <- "somma" 
  return(res$somma) # ritorna il totale per ogni rigA (numero di count di geni complessivo per ogni cellula)
}
  



set_initial_size<-function(adata, layers = NA){ 
    # Set current counts per observation of a layer as its initial size.
    # The initial size is only set if it does not already exist.
    # Arguments
    # ---------
    # adata
    #     Annotated data matrix.
    # layers
    #     Name of layers for which to calculate initial size.
    # Returns
    # -------
    # None
  if(is.na(layers)){
    layers = c("spliced", "unspliced")
  }
  l <- list()
  for(layer in layers){
    initial_size_layer = paste("initial_size_", layer, sep = "")
    if(layer %in% adata$layers$keys() & ! eval(initial_size_layer) %in% names(adata$obs))
      adata$obs[[eval(initial_size_layer)]] = get_size(adata, layer)
  
  }
  return(adata)
}

cluster_genesUMAP <- function(unspliced, spliced, typeCell, n_clusters){
  # Given the unspliced and spliced counts, a UMAP projection is performed starting from all the cells. Then, for each type of cells, we consider only the cells of that type and we apply the k-means algorithm starting from the two components given by UMAP. The number of clusters used in the algorithm is given as input. 
  # We return the clustering that we have obtained for each type. 
  # Arguments
  # --------
  # unspliced
  #   Matrix with unspliced counts
  # spliced
  #   Matrix with spliced counts
  # typeCell
  #   Vector of dimension n_cells with the labels of the cells
  # n_clusters
  #   Number of subclasses for each type
  # ---------
  # Output 
  #   Subtypes returned by the k-means algorithm

  ##############################################################
  # UMAP 
  ##############################################################
  seed = 1234
  set.seed(seed)
  data <- cbind(unspliced, spliced)
  data.umap <- umap(data)

  umap_df <- as.data.frame(data.umap$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")

  if(typeof(typeCell) == "character"){
    typeCell <- numericTypeCell(typeCell)
  }
  umap_df$typeCell <- as.character(typeCell)

  ##############################################################
  ##############################################################
  # KMEANS con distanza EUCLIDEA
  ##############################################################
  kmeans.umap <- list()



  for(ty in unique(typeCell)){
      set.seed(1234 + ty)
      cellTy <- which(typeCell == ty)
      kmeans.umap[[ty]] <- kmeans(umap_df[cellTy,1:2], centers = n_clusters, iter.max = 10, nstart = 25)
  }

  subtypeCell_UMAP <- c()
  for(ty in unique(typeCell)){ 
    subtypeCell_UMAP <- c(subtypeCell_UMAP, kmeans.umap[[ty]]$cluster + n_clusters*(ty - 1))
  }

  return(subtypeCell_UMAP)

}


numericTypeCell <- function(typeCell){
  # Convert the vector with the types of cells from a string vectors into a integer vector
  # Arguments
  # ---------
  # typeCell
  #   Vector of dimension n_cells with the different types of cells. Each type is associated with a string
  # ---------
  # Output 
  # Vector of dimension n_cells with the different types of cells. Each type is associated with an integer
  types <- unique(typeCell)
  i <- 1
  for(ty in types){
    typeCell[typeCell == ty] <- i
    i <- i + 1
  }

  typeCell <- as.integer(typeCell)
  return(typeCell)
}


loadRealData <- function(dataset, typeProcessing, pathData = NULL, names = FALSE, sparse = FALSE){
  if(is.null(pathData)){
    error("Provide the directory where the real data are stored")
  }
  
  pathData <- paste0(pathData, "/", typeProcessing)
  
  # genes'names
  var <- loadVar(pathData)
  # cells'names
  obs <- loadObs(pathData)  

  if(typeProcessing == "moments"){ # import pre-processed data
    unspliced <- loadMu(pathData, names, sparse)
    spliced   <- loadMs(pathData, names, sparse)
  }else{  # import raw counts
    unspliced <- loadUnspliced(pathData, names, sparse)
    spliced   <- loadSpliced(pathData, names, sparse)
  }

  # import the type of cells 
  typeCell <- obs$clusters

  res <- list("unspliced" = unspliced, "spliced" = spliced, "typeCell" = typeCell)
  return(res)

}


loadVar <- function(pathData){
  # genes'names
  var <- read.csv(file = paste(pathData, '/var.csv', sep =""))
  return(var)
}

loadObs <- function(pathData){
  # nomi delle cellule
  obs <- read.csv(file = paste(pathData, '/obs.csv', sep = "")) 
  return(obs)
}

loadAdata <- function(pathData, obs = NULL, var = NULL, sparse = FALSE){
  if(is.null(obs)){
    obs <- loadObs(pathData)
  }
  if(is.null(var)){
    var <- loadVar(pathData)
  }
  # import Data e transform them into a matrix
  adata <- read.table(paste(pathData, '/X.csv', sep = ""), header = FALSE, sep = ",",row.names = obs$index, col.names = var$index)
  adata <- as(adata, "matrix") 

  if(sparse){
    adata <- as(adata, "dgCMatrix")
  }

  return(adata)
}
  
loadUnspliced <- function(pathData, names = FALSE, sparse = FALSE){
  # impot unspliced raw data and transform them into a spare matrix, if required
  unspliced <- read.csv(file = paste(pathData, '/unspliced.csv', sep = ""))

  if(names){
    var <- loadVar(pathData)
    obs <- loadObs(pathData) 
    colnames(unspliced) <-  var$index
    rownames(unspliced) <- obs$index
  }
 
  unspliced <- as(unspliced, "matrix") 

  if(sparse){
    unspliced <- as(unspliced, "dgCMatrix")
  }

  return(unspliced)
}

 
loadSpliced <- function(pathData, names = FALSE, sparse = FALSE){
  # impot spliced raw data and transform them into a spare matrix, if required
  spliced <- read.csv(file = paste(pathData, '/spliced.csv', sep = ""))

  if(names){
    var <- loadVar(pathData)
    obs <- loadObs(pathData) 
    colnames(spliced) <-  var$index
    rownames(spliced) <- obs$index
  }
 
  spliced <- as(spliced, "matrix") 

  if(sparse){
    spliced <- as(spliced, "dgCMatrix")
  }

  return(spliced)
}

loadMu <- function(pathData, names = FALSE, sparse = FALSE){
  # impot unspliced pre-processed moments and transform them into a spare matrix, if required
  Mu <- read.csv(file = paste(pathData, '/Mu.csv', sep = ""))

  if(names){
    var <- loadVar(pathData)
    obs <- loadObs(pathData) 
    colnames(Mu) <-  var$index
    rownames(Mu) <- obs$index
  }
 
  Mu <- as(Mu, "matrix") 

  if(sparse){
    Mu <- as(Mu, "dgCMatrix")
  }

  return(Mu)
}

loadMs <- function(pathData, names = FALSE, sparse = FALSE){
  # impot spliced pre-processed moments and transform them into a spare matrix, if required
  Ms <- read.csv(file = paste(pathData, '/Ms.csv', sep = ""))

  if(names){
    var <- loadVar(pathData)
    obs <- loadObs(pathData) 
    colnames(Ms) <-  var$index
    rownames(Ms) <- obs$index
  }
 
  Ms <- as(Ms, "matrix") 

  if(sparse){
    Ms <- as(Ms, "dgCMatrix")
  }

  return(Ms)

}


real_in_CI <- function(real, CI){
  return((real >= CI[1,]) * (real <= CI[2, ]))
}


boxplotMCMC <- function(chain, real, namePar, nameIndex ){
  df <- data.frame(par = as.vector(t(chain)), index = rep(as.character(seq(1, dim(chain)[1])), each = dim(chain)[2]))

  boxplot <- list()

  if(dim(chain)[1] < 100){
    set <- seq(1, dim(chain)[1])
  }else{
    set <- seq(1, dim(chain)[1], 100)
  }
  for(j in 1:length(set)){
    i <- set[j]
    sub <- subset(df, index == as.character(i))
    
    boxplot[[j]] <- ggplot(sub, aes(x = index, y = par)) +
        geom_boxplot() +
        labs(title = paste(namePar, ",", nameIndex, i), x = paste(nameIndex, i), y = "") + 
        geom_point(x = 1, y = real[i], color = "red", size = 5) + 
        ylim(min(c(sub$par, real[i])), max(c(sub$par, real[i]))) 
    
  }

  return(boxplot)
}


boxplotMCMCinx <- function(chain, real, namePar, nameIndex ){
  if(!is.null(dim(chain)[1])){
    if(dim(chain)[1] < 50){
      set <- seq(1, dim(chain)[1])
    }else{
      set <- seq(1, 50)
    }
  }else{
    attr(chain, "dim") <- c(1, length(chain))
    set <- 1
  }
  
  df <- data.frame(par = as.vector(t(chain[set,])), index = rep(set, each = dim(chain)[2]))
  dfReal <- data.frame(par = real[set], index = set)
  ggplot(df, aes(group = index, y = par, x = index)) +
        geom_boxplot() + geom_point(data = dfReal, aes(x = index, y = par), col = "red", size = 5)+ 
        labs(title = namePar, x = nameIndex, y = "") 
 
}

boxplotMCMCwithCI <- function(chain, real, namePar, nameIndex, alpha = 0.05){
  if(dim(chain)[1] < 50){
    set <- seq(1, dim(chain)[1])
  }else{
    set <- seq(1, 50)
  }

  minChain <- apply(chain[set,], FUN = min, MARGIN = c(1))
  maxChain <- apply(chain[set,], FUN = max, MARGIN = c(1))
  qCI <- apply(chain[set,], FUN = quantile, MARGIN = c(1), probs = c(alpha/2, 1-alpha/2))
  q25Chain <- apply(chain[set,], FUN = quantile, MARGIN = c(1), probs = 0.25)
  q50Chain <- apply(chain[set,], FUN = quantile, MARGIN = c(1), probs = 0.50)
  q75Chain <- apply(chain[set,], FUN = quantile, MARGIN = c(1), probs = 0.75)

  
  df <- data.frame(par = as.vector(t(chain[set,])), index = rep(set, each = dim(chain)[2]))
  
  dfLim <- data.frame(index = set, y0 = minChain, y100 = maxChain, y25 = q25Chain, y50 = q50Chain, y75 = q75Chain, yCIlow = qCI[1,], yCIup = qCI[2,])
  
  dfReal <- data.frame(par = real[set], index = set)

  ggplot(dfLim, aes(x = index)) + geom_point(data = df, aes(x = index, y = par), col = "gray") +
  geom_boxplot(
   aes(group = as.factor(index), ymin = yCIlow, lower = y25, middle = y50, upper = y75, ymax = yCIup),
   stat = "identity", color = "black"
 ) + labs(title = paste(namePar, "boxplot with CI of level", alpha), x = nameIndex, y = "") +  geom_point(data = dfReal, aes(x = index, y = par), col = "red", size = 5)
 
}



RNA_velocity <- function(beta, gamma, pos_u, pos_s){
  # beta is a n_genes vector
  if(is.null(dim(beta))){
    gamma <- rep(gamma, dim(pos_u)[1])
    gamma <- matrix(gamma, nrow = dim(pos_u)[1], byrow = TRUE)
    beta <- rep(beta, dim(pos_u)[1])
    beta <- matrix(beta, nrow = dim(pos_u)[1], byrow = TRUE)
    
    res <- pos_u*beta - pos_s*gamma
  }else{
    dimRow <- dim(pos_u)[1]
    dimCol <- dim(beta)[1]
    SampleToSave <- dim(beta)[2]
    beta <- as.vector(beta)
    beta <- rep(beta, dimRow)
    beta <- matrix(beta, nrow = dimRow, byrow = TRUE)
    gamma <- as.vector(gamma)
    gamma <- rep(gamma, dimRow)
    gamma <- matrix(gamma, nrow = dimRow, byrow = TRUE)

    pos_u <- matrix(pos_u, nrow = dimRow, ncol = dimCol * SampleToSave)
    pos_s <- matrix(pos_s, nrow = dimRow, ncol = dimCol * SampleToSave)

    res <- pos_u*beta - pos_s*gamma
    attr(res, "dim") <- c(dimRow, dimCol, SampleToSave)
  }
  return(res)
}



nomeSIM_plot <- function(name){
  label <- ""
  if(grepl("SW1", name)){
      label <- paste(label,"Sw$_1$", sep = "")
  }else{
      label <- paste(label,"Sw$_{10}$", sep = "")
  }

  if(grepl("T1", name)){
      label <- paste(label,"$T_1$", sep = "-")
  }else if(grepl("T2", name)){
      label <- paste(label,"$T_3$", sep = "-")
  }else if(grepl("T3", name)){
      label <- paste(label,"$T_{10}$", sep = "-")
  }else if(grepl("T4", name)){
      label <- paste(label,"$T_{sc}$", sep = "-")
  }

  if(grepl("D1", name)){
      label <- paste(label,"Poi", sep = "-")
  }else if(grepl("D4", name)){
      label <- paste(label,"NB", sep = "-")
  }
  return(label)
}



nomeSIM_ggpairs <- function(name){
label <- ""
  if(grepl("SW1", name)){
      label <- paste(label,"Sw[1]", sep = "")
  }else{
      label <- paste(label,"Sw[10]", sep = "")
  }

  if(grepl("T1", name)){
      label <- paste(label,"T[1]", sep = "-")
  }else if(grepl("T2", name)){
      label <- paste(label,"T[3]", sep = "-")
  }else if(grepl("T3", name)){
      label <- paste(label,"T[10]", sep = "-")
  }else if(grepl("T4", name)){
      label <- paste(label,"T[sc]", sep = "-")
  }

  if(grepl("D1", name)){
      label <- paste(label,"Poi", sep = "-")
  }else if(grepl("D4", name)){
      label <- paste(label,"NB", sep = "-")
  }

  return(label)
}

