library(mvtnorm)
library(onlineFDR)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(Matrix)


result_list_phi_psi <- list()

cal.FDP <- function(x.real, x.result, t){
  if (sum(x.result[1:t]) != 0) {
    fdp <- sum((1 - x.real[1:t]) * x.result[1:t]) / sum(x.result[1:t])
  } else {
    fdp <- 0
  }
  return(fdp)
}

cal.MDP <- function(x.real, x.result, t){
  if (sum(x.real[1:t]) != 0) {
    mdp <- sum((x.real[1:t]) * (1 - x.result[1:t])) / sum(x.real[1:t])
  } else {
    mdp <- 0
  }
  return(mdp)
}

calculate_p_e <- function(seed = 2025, pi1 = 0.2, t = 500, rho = 0.5, lag= 30, muc = 3){
  set.seed(seed)
  sigma <- 1
  out.loc <- sample(c(0, 1), t, replace = TRUE, prob = c(1 - pi1, pi1))
  mu <- rep(0, t)
  mu[which(out.loc == 1)] <- muc
  
  Sigma <- diag(rep(1, t))
  for (i in 1:lag) {
    for (j in 1:(t - i)) {
      Sigma[j, j + i] <- rho^i
      Sigma[j + i, j] <- rho^i
    }
  }
  X <- rmvnorm(1, mean = mu, sigma = Sigma)
  
  p <- rep(1, t)
  e <- rep(0, t)
  vec_lag_0 <- rho^seq(lag, 1, by = -1)
  
  p[1] <- 1 - pnorm(X[1], 0, sd = 1)
  e[1] <- exp(muc * X[i] - muc^2 / 2)
  
  for(i in 2:t) {
    if (i <= lag + 1) {
      mu_cond <- t(as.matrix(vec_lag_0[(lag - i + 2):lag])) %*% solve(Sigma[1:(i - 1), 1:(i - 1)]) %*% mu[1:(i - 1)]
      sigma_cond <- 1 - t(as.matrix(vec_lag_0[(lag - i + 2):lag])) %*% solve(Sigma[1:(i - 1), 1:(i - 1)]) %*% as.matrix(vec_lag_0[(lag - i + 2):lag])
    } else {
      chol_Sigma <- Cholesky(Sigma[1:(i - 1), 1:(i - 1)], LDL = FALSE)
      I_i <- Diagonal(n = (i-1))         
      Sigma_inv <- solve(chol_Sigma, I_i) 
      Sigma_inv <- as.matrix(Sigma_inv)
      mu_cond <- t(as.matrix(c(rep(0, (i - lag - 1)), vec_lag_0))) %*% Sigma_inv %*% mu[1:(i - 1)]
      sigma_cond <- 1 - t(as.matrix(vec_lag_0)) %*% Sigma_inv[ (i - lag):(i - 1), (i - lag):(i - 1)] %*% as.matrix(vec_lag_0)
    }
    p[i] <- 1 - pnorm(X[i], mu_cond, sd = sqrt(sigma_cond))
    e[i] <- exp(muc * (X[i] - mu_cond) / sigma_cond - muc^2 / (2 * sigma_cond))
  }
  out <- list(loc=out.loc, p=p, e=e)
  return(out)
}

one.exper <- function(
    seed,
    alpha = 0.05,
    out.loc,
    p,
    e,
    t   = 500,
    rho = 0.5,
    reward.omega = 0.5,
    penalty.omega = 0.5,
    init.omega = 0.005,
    lag = 30,
    lambda = 0.1,
    muc = 3,
    d = 0.9,
    epsilon.a = 0.2,
    epsilon.r = 1
){
  set.seed(seed)
  
  sigma <- 1
  out.loc <- out.loc
  p <- p
  e <- e

  rej.eLOND    <- rep(0, t)
  rej.eSAFFRON     <- rep(0, t)
  rej.eLORD    <- rep(0, t)
  rej.pS_RAI   <- rep(0, t)
  rej.pL_RAI  <- rep(0, t)
  
  gamma.eLOND   <- rep(0, t)
  gamma.eSAFFRON    <- rep(0, t)
  gamma.eLORD   <- rep(0, t)
  gamma.pS_RAI  <- rep(0, t)
  gamma.pL_RAI <- rep(0, t)
  
  level.eLOND   <- rep(0, t)
  level.eSAFFRON    <- rep(0, t)
  level.eLORD   <- rep(0, t)
  level.pS_RAI  <- rep(0, t)
  level.pL_RAI <- rep(0, t)
  
  ## GAI
  start.LORD <- proc.time()
  rej.LORD <- LORD(as.vector(p), alpha = alpha)$R
  time.LORD <- (proc.time() - start.LORD)[["elapsed"]]
  
  start.SAFFRON <- proc.time()
  rej.SAFFRON <- SAFFRON(as.vector(p), alpha = alpha)$R
  time.SAFFRON <- (proc.time() - start.SAFFRON)[["elapsed"]]
  
  ## SupLORD
  start.suplord <- proc.time()
  rej.suplord <- supLORD(as.vector(p), delta = 0.05 , eps = 0.07620258 , r = 30, eta = 0.05, rho = 30)$R
  time.suplord <- (proc.time() - start.suplord)[["elapsed"]]
  
  ## e-LOND
  start.eLOND <- proc.time()
  j <- 1
  gamma.eLOND[j] <- 1 / (j * (j + 1))
  level.eLOND[j] <- alpha * gamma.eLOND[j]
  rej.eLOND[j] <- round(e[j] > 1 / level.eLOND[j])
  for (j in 2:t) {
    gamma.eLOND[j] <- 1 / (j * (j + 1))
    level.eLOND[j] <- alpha * gamma.eLOND[j] * (sum(rej.eLOND[1:(j - 1)]) + 1)
    rej.eLOND[j] <- round(e[j] > 1 / level.eLOND[j])
  }
  time.eLOND <- (proc.time() - start.eLOND)[["elapsed"]]
  
  ## eSAFFRON: e-LORD
  cum.rej.eLORD <- 0
  rem.rej.eLORD <- 0
  start.eLORD <- proc.time()
  j <- 1
  gamma.eLORD[j] <- init.omega
  level.eLORD[j] <- gamma.eLORD[j] * alpha
  rej.eLORD[j] <- round(e[j] > 1 / level.eLORD[j])
  cum.rej.eLORD <- rej.eLORD[j]
  rem.rej.eLORD <- alpha - level.eLORD[j]
  for (j in 2:t) {
    gamma.eLORD[j] <- gamma.eLORD[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.eLORD) * (1 - rej.eLORD[j - 1]) -
      init.omega * reward.omega^(cum.rej.eLORD) * rej.eLORD[j - 1]
    level.eLORD[j] <- gamma.eLORD[j] * rem.rej.eLORD * (cum.rej.eLORD + 1)
    if (is.na(level.eLORD[j])) {
      level.eLORD[j] <- 0
    }
    rej.eLORD[j] <- round(e[j] > 1 / level.eLORD[j])
    rem.rej.eLORD <- rem.rej.eLORD - level.eLORD[j] / (cum.rej.eLORD + 1)
    cum.rej.eLORD <- cum.rej.eLORD + rej.eLORD[j]
  }
  time.eLORD <- (proc.time() - start.eLORD)[["elapsed"]]
  
  ## pL_RAI: pL-RAI
  cum.rej.pL_RAI <- 0
  rem.rej.pL_RAI <- 0
  start.pL_RAI <- proc.time()
  j <- 1
  gamma.pL_RAI[j] <- init.omega
  level.pL_RAI[j] <- gamma.pL_RAI[j] * alpha
  rej.pL_RAI[j] <- round(p[j] <= level.pL_RAI[j])
  cum.rej.pL_RAI <- rej.pL_RAI[j]
  rem.rej.pL_RAI <- alpha - level.pL_RAI[j]
  for (j in 2:t) {
    gamma.pL_RAI[j] <- gamma.pL_RAI[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.pL_RAI) * (1 - rej.pL_RAI[j - 1]) -
      init.omega * reward.omega^(cum.rej.pL_RAI) * rej.pL_RAI[j - 1]
    level.pL_RAI[j] <- gamma.pL_RAI[j] * rem.rej.pL_RAI * (cum.rej.pL_RAI + 1)
    if (is.na(level.pL_RAI[j])) {
      level.pL_RAI[j] <- 0
    }
    rej.pL_RAI[j] <- round(p[j] <= level.pL_RAI[j])
    rem.rej.pL_RAI <- rem.rej.pL_RAI - level.pL_RAI[j] / (cum.rej.pL_RAI + 1)
    cum.rej.pL_RAI <- cum.rej.pL_RAI + rej.pL_RAI[j]
  }
  time.pL_RAI <- (proc.time() - start.pL_RAI)[["elapsed"]]
  
  ## eSAFFRON: e-SAFFRON
  cum.rej.eSAFFRON <- 0
  rem.rej.eSAFFRON <- 0
  start.eSAFFRON <- proc.time()
  j <- 1
  gamma.eSAFFRON[j] <- init.omega
  level.eSAFFRON[j] <- gamma.eSAFFRON[j] * alpha * (1 - lambda)
  rej.eSAFFRON[j] <- round(e[j] > 1 / level.eSAFFRON[j])
  cum.rej.eSAFFRON <- rej.eSAFFRON[j]
  rem.rej.eSAFFRON <- alpha * (1 - lambda) - level.eSAFFRON[j] * round(e[j] < 1 / lambda)
  for (j in 2:t) {
    gamma.eSAFFRON[j] <- gamma.eSAFFRON[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.eSAFFRON) * (1 - rej.eSAFFRON[j - 1]) -
      init.omega * reward.omega^(cum.rej.eSAFFRON) * rej.eSAFFRON[j - 1]
    level.eSAFFRON[j] <- gamma.eSAFFRON[j] * rem.rej.eSAFFRON * (cum.rej.eSAFFRON + 1)
    if (is.na(level.eSAFFRON[j])) {
      level.eSAFFRON[j] <- 0
    }
    rej.eSAFFRON[j] <- round(e[j] > 1 / level.eSAFFRON[j])
    rem.rej.eSAFFRON <- rem.rej.eSAFFRON - round(e[j] < 1 / lambda) * level.eSAFFRON[j] / (cum.rej.eSAFFRON + 1)
    cum.rej.eSAFFRON <- cum.rej.eSAFFRON + rej.eSAFFRON[j]
  }
  time.eSAFFRON <- (proc.time() - start.eSAFFRON)[["elapsed"]]
  
  ## pS_RAI: pS-RAI
  cum.rej.pS_RAI <- 0
  rem.rej.pS_RAI <- 0
  start.pS_RAI <- proc.time()
  j <- 1
  gamma.pS_RAI[j] <- init.omega
  level.pS_RAI[j] <- gamma.pS_RAI[j] * alpha * (1 - lambda)
  rej.pS_RAI[j] <- round(p[j] <= level.pS_RAI[j])
  cum.rej.pS_RAI <- rej.pS_RAI[j]
  rem.rej.pS_RAI <- alpha * (1 - lambda) - level.pS_RAI[j] * round(p[j] > lambda)
  for (j in 2:t) {
    gamma.pS_RAI[j] <- gamma.pS_RAI[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.pS_RAI) * (1 - rej.pS_RAI[j - 1]) -
      init.omega * reward.omega^(cum.rej.pS_RAI) * rej.pS_RAI[j - 1]
    level.pS_RAI[j] <- gamma.pS_RAI[j] * rem.rej.pS_RAI * (cum.rej.pS_RAI + 1)
    if (is.na(level.pS_RAI[j])) {
      level.pS_RAI[j] <- 0
    }
    rej.pS_RAI[j] <- round(p[j] <= level.pS_RAI[j])
    rem.rej.pS_RAI <- rem.rej.pS_RAI - round(p[j] > lambda) * level.pS_RAI[j] / (cum.rej.pS_RAI + 1)
    cum.rej.pS_RAI <- cum.rej.pS_RAI + rej.pS_RAI[j]
  }
  time.pS_RAI <- (proc.time() - start.pS_RAI)[["elapsed"]]
  
  j <- t
  FDPs <- c(
    cal.FDP(out.loc[1:j], rej.eLORD[1:j], j),
    cal.FDP(out.loc[1:j], rej.eSAFFRON[1:j], j),
    cal.FDP(out.loc[1:j], rej.pL_RAI[1:j], j),
    cal.FDP(out.loc[1:j], rej.pS_RAI[1:j], j),
    cal.FDP(out.loc[1:j], rej.eLOND[1:j], j),
    cal.FDP(out.loc[1:j], rej.LORD[1:j], j),
    cal.FDP(out.loc[1:j], rej.SAFFRON[1:j], j),
    cal.FDP(out.loc[1:j], rej.suplord[1:j], j)
  )
  MDPs <- c(
    cal.MDP(out.loc[1:j], rej.eLORD[1:j], j),
    cal.MDP(out.loc[1:j], rej.eSAFFRON[1:j], j),
    cal.MDP(out.loc[1:j], rej.pL_RAI[1:j], j),
    cal.MDP(out.loc[1:j], rej.pS_RAI[1:j], j),
    cal.MDP(out.loc[1:j], rej.eLOND[1:j], j),
    cal.MDP(out.loc[1:j], rej.LORD[1:j], j),
    cal.MDP(out.loc[1:j], rej.SAFFRON[1:j], j),
    cal.MDP(out.loc[1:j], rej.suplord[1:j], j)
  )
  
  output <- c(FDPs, MDPs,
              time.eLORD, time.eSAFFRON, time.pL_RAI, time.pS_RAI,
              time.eLOND, time.LORD, time.SAFFRON, time.suplord)
  
  result_colnames <- c("FDP_eLORD", "FDP_eSAFFRON", "FDP_pL_RAI", "FDP_pS_RAI", 
                       "FDP_eLOND", "FDP_LORD", "FDP_SAFFRON", "FDP_suplord",
                       "MDP_eLORD", "MDP_eSAFFRON", "MDP_pL_RAI", "MDP_pS_RAI", 
                       "MDP_eLOND", "MDP_LORD", "MDP_SAFFRON", "MDP_suplord")
  names(output)[1:16] <- result_colnames

  return(output)
}


multi.experi <- function(pi1 = 0.2,
                         reward.omega = 0.5,
                         penalty.omega = 0.5,
                         init.omega = 0.005,
                         precomp,
                         seeds     = 2026:2075,
                         alpha     = 0.05,
                         t         = 500,
                         rho       = 0.5,
                         lag       = 30,
                         lambda    = 0.1,
                         muc       = 3) {

  pi1           <- as.numeric(pi1)
  reward.omega  <- as.numeric(reward.omega)
  penalty.omega <- as.numeric(penalty.omega)
  init.omega    <- as.numeric(init.omega)
  
  r <- length(seeds)
  

  mat <- sapply(seeds, function(seed) {
    key <- paste0("pi1_", pi1, "_seed_", seed)
    pe  <- precomp$vary_pi1_rho03[[key]]

    loc_vec <- as.numeric(pe$loc)
    p_vec   <- as.numeric(pe$p)
    e_vec   <- as.numeric(pe$e)
    
    one.exper(
      seed         = seed,
      alpha        = alpha,
      out.loc      = loc_vec,
      p            = p_vec,
      e            = e_vec,
      t            = t,
      rho          = rho,
      reward.omega = reward.omega,
      penalty.omega= penalty.omega,
      init.omega   = init.omega,
      lag          = lag,
      lambda       = lambda,
      muc          = muc
    )
  })
  mat <- as.data.frame(mat)

  row_names <- c(

    "FDP_eLORD",   "FDP_eSAFFRON",   "FDP_pL_RAI", "FDP_pS_RAI",
    "FDP_eLOND",   "FDP_LORD",   "FDP_SAFFRON", "FDP_suplord",

    "MDP_eLORD",   "MDP_eSAFFRON",   "MDP_pL_RAI", "MDP_pS_RAI",
    "MDP_eLOND",   "MDP_LORD",   "MDP_SAFFRON", "MDP_suplord",

    "time_eLORD",  "time_eSAFFRON",  "time_pL_RAI", "time_pS_RAI",
    "time_eLOND",  "time_LORD",  "time_SAFFRON", "time_suplord"
  )
  
  stopifnot(length(row_names) == nrow(mat))  
  rownames(mat) <- row_names


  fdp_idx <- grep("^FDP_", rownames(mat))
  mdp_idx <- grep("^MDP_", rownames(mat))
  target_idx <- c(fdp_idx, mdp_idx)

  avg <- rowMeans(mat)

  se_idx <- grep("^(FDP|MDP)_", names(avg))

  se_vals <- apply(mat[se_idx, , drop = FALSE],
                   1, sd) / sqrt(r)
  names(se_vals) <- paste0("SE_", names(avg)[se_idx])
  
  return( c(avg, se_vals) )
}
