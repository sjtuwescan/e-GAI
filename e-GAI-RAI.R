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

  rej.elond    <- rep(0, t)
  rej.ours     <- rep(0, t)
  rej.ours0    <- rep(0, t)
  rej.ours.p   <- rep(0, t)
  rej.ours0.p  <- rep(0, t)
  
  gamma.elond   <- rep(0, t)
  gamma.ours    <- rep(0, t)
  gamma.ours0   <- rep(0, t)
  gamma.ours.p  <- rep(0, t)
  gamma.ours0.p <- rep(0, t)
  
  level.elond   <- rep(0, t)
  level.ours    <- rep(0, t)
  level.ours0   <- rep(0, t)
  level.ours.p  <- rep(0, t)
  level.ours0.p <- rep(0, t)
  
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
  start.elond <- proc.time()
  j <- 1
  gamma.elond[j] <- 1 / (j * (j + 1))
  level.elond[j] <- alpha * gamma.elond[j]
  rej.elond[j] <- round(e[j] > 1 / level.elond[j])
  for (j in 2:t) {
    gamma.elond[j] <- 1 / (j * (j + 1))
    level.elond[j] <- alpha * gamma.elond[j] * (sum(rej.elond[1:(j - 1)]) + 1)
    rej.elond[j] <- round(e[j] > 1 / level.elond[j])
  }
  time.elond <- (proc.time() - start.elond)[["elapsed"]]
  
  ## Ours: e-LORD
  cum.rej.ours0 <- 0
  rem.rej.ours0 <- 0
  start.ours0 <- proc.time()
  j <- 1
  gamma.ours0[j] <- init.omega
  level.ours0[j] <- gamma.ours0[j] * alpha
  rej.ours0[j] <- round(e[j] > 1 / level.ours0[j])
  cum.rej.ours0 <- rej.ours0[j]
  rem.rej.ours0 <- alpha - level.ours0[j]
  for (j in 2:t) {
    gamma.ours0[j] <- gamma.ours0[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.ours0) * (1 - rej.ours0[j - 1]) -
      init.omega * reward.omega^(cum.rej.ours0) * rej.ours0[j - 1]
    level.ours0[j] <- gamma.ours0[j] * rem.rej.ours0 * (cum.rej.ours0 + 1)
    if (is.na(level.ours0[j])) {
      level.ours0[j] <- 0
    }
    rej.ours0[j] <- round(e[j] > 1 / level.ours0[j])
    rem.rej.ours0 <- rem.rej.ours0 - level.ours0[j] / (cum.rej.ours0 + 1)
    cum.rej.ours0 <- cum.rej.ours0 + rej.ours0[j]
  }
  time.ours0 <- (proc.time() - start.ours0)[["elapsed"]]
  
  ## Ours0.p: conditional-p-LORD
  cum.rej.ours0.p <- 0
  rem.rej.ours0.p <- 0
  start.ours0.p <- proc.time()
  j <- 1
  gamma.ours0.p[j] <- init.omega
  level.ours0.p[j] <- gamma.ours0.p[j] * alpha
  rej.ours0.p[j] <- round(p[j] <= level.ours0.p[j])
  cum.rej.ours0.p <- rej.ours0.p[j]
  rem.rej.ours0.p <- alpha - level.ours0.p[j]
  for (j in 2:t) {
    gamma.ours0.p[j] <- gamma.ours0.p[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.ours0.p) * (1 - rej.ours0.p[j - 1]) -
      init.omega * reward.omega^(cum.rej.ours0.p) * rej.ours0.p[j - 1]
    level.ours0.p[j] <- gamma.ours0.p[j] * rem.rej.ours0.p * (cum.rej.ours0.p + 1)
    if (is.na(level.ours0.p[j])) {
      level.ours0.p[j] <- 0
    }
    rej.ours0.p[j] <- round(p[j] <= level.ours0.p[j])
    rem.rej.ours0.p <- rem.rej.ours0.p - level.ours0.p[j] / (cum.rej.ours0.p + 1)
    cum.rej.ours0.p <- cum.rej.ours0.p + rej.ours0.p[j]
  }
  time.ours0.p <- (proc.time() - start.ours0.p)[["elapsed"]]
  
  ## Ours: e-SAFFRON
  cum.rej.ours <- 0
  rem.rej.ours <- 0
  start.ours <- proc.time()
  j <- 1
  gamma.ours[j] <- init.omega
  level.ours[j] <- gamma.ours[j] * alpha * (1 - lambda)
  rej.ours[j] <- round(e[j] > 1 / level.ours[j])
  cum.rej.ours <- rej.ours[j]
  rem.rej.ours <- alpha * (1 - lambda) - level.ours[j] * round(e[j] < 1 / lambda)
  for (j in 2:t) {
    gamma.ours[j] <- gamma.ours[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.ours) * (1 - rej.ours[j - 1]) -
      init.omega * reward.omega^(cum.rej.ours) * rej.ours[j - 1]
    level.ours[j] <- gamma.ours[j] * rem.rej.ours * (cum.rej.ours + 1)
    if (is.na(level.ours[j])) {
      level.ours[j] <- 0
    }
    rej.ours[j] <- round(e[j] > 1 / level.ours[j])
    rem.rej.ours <- rem.rej.ours - round(e[j] < 1 / lambda) * level.ours[j] / (cum.rej.ours + 1)
    cum.rej.ours <- cum.rej.ours + rej.ours[j]
  }
  time.ours <- (proc.time() - start.ours)[["elapsed"]]
  
  ## Ours.p: conditional-p-SAFFRON
  cum.rej.ours.p <- 0
  rem.rej.ours.p <- 0
  start.ours.p <- proc.time()
  j <- 1
  gamma.ours.p[j] <- init.omega
  level.ours.p[j] <- gamma.ours.p[j] * alpha * (1 - lambda)
  rej.ours.p[j] <- round(p[j] <= level.ours.p[j])
  cum.rej.ours.p <- rej.ours.p[j]
  rem.rej.ours.p <- alpha * (1 - lambda) - level.ours.p[j] * round(p[j] > lambda)
  for (j in 2:t) {
    gamma.ours.p[j] <- gamma.ours.p[j - 1] + init.omega * penalty.omega^(j - 1 - cum.rej.ours.p) * (1 - rej.ours.p[j - 1]) -
      init.omega * reward.omega^(cum.rej.ours.p) * rej.ours.p[j - 1]
    level.ours.p[j] <- gamma.ours.p[j] * rem.rej.ours.p * (cum.rej.ours.p + 1)
    if (is.na(level.ours.p[j])) {
      level.ours.p[j] <- 0
    }
    rej.ours.p[j] <- round(p[j] <= level.ours.p[j])
    rem.rej.ours.p <- rem.rej.ours.p - round(p[j] > lambda) * level.ours.p[j] / (cum.rej.ours.p + 1)
    cum.rej.ours.p <- cum.rej.ours.p + rej.ours.p[j]
  }
  time.ours.p <- (proc.time() - start.ours.p)[["elapsed"]]
  
  j <- t
  FDPs <- c(
    cal.FDP(out.loc[1:j], rej.ours0[1:j], j),
    cal.FDP(out.loc[1:j], rej.ours[1:j], j),
    cal.FDP(out.loc[1:j], rej.ours0.p[1:j], j),
    cal.FDP(out.loc[1:j], rej.ours.p[1:j], j),
    cal.FDP(out.loc[1:j], rej.elond[1:j], j),
    cal.FDP(out.loc[1:j], rej.LORD[1:j], j),
    cal.FDP(out.loc[1:j], rej.SAFFRON[1:j], j),
    cal.FDP(out.loc[1:j], rej.suplord[1:j], j)
  )
  MDPs <- c(
    cal.MDP(out.loc[1:j], rej.ours0[1:j], j),
    cal.MDP(out.loc[1:j], rej.ours[1:j], j),
    cal.MDP(out.loc[1:j], rej.ours0.p[1:j], j),
    cal.MDP(out.loc[1:j], rej.ours.p[1:j], j),
    cal.MDP(out.loc[1:j], rej.elond[1:j], j),
    cal.MDP(out.loc[1:j], rej.LORD[1:j], j),
    cal.MDP(out.loc[1:j], rej.SAFFRON[1:j], j),
    cal.MDP(out.loc[1:j], rej.suplord[1:j], j)
  )
  
  output <- c(FDPs, MDPs,
              time.ours0, time.ours, time.ours0.p, time.ours.p,
              time.elond, time.LORD, time.SAFFRON, time.suplord)
  
  result_colnames <- c("FDP_eLORD", "FDP_eSAFFRON", "FDP_pL_RAI", "FDP_pS_RAI", 
                       "FDP_elond", "FDP_LORD", "FDP_SAFFRON", "FDP_suplord",
                       "MDP_eLORD", "MDP_eSAFFRON", "MDP_pL_RAI", "MDP_pS_RAI", 
                       "MDP_elond", "MDP_LORD", "MDP_SAFFRON", "MDP_suplord")
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

    "FDP_ours0",   "FDP_ours",   "FDP_ours0.p", "FDP_ours.p",
    "FDP_elond",   "FDP_LORD",   "FDP_SAFFRON", "FDP_suplord",

    "MDP_ours0",   "MDP_ours",   "MDP_ours0.p", "MDP_ours.p",
    "MDP_elond",   "MDP_LORD",   "MDP_SAFFRON", "MDP_suplord",

    "time_ours0",  "time_ours",  "time_ours0.p", "time_ours.p",
    "time_elond",  "time_LORD",  "time_SAFFRON", "time_suplord"
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
