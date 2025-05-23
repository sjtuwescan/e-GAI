rm(list = ls())
library(mvtnorm)
library(onlineFDR)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(Matrix)

cal.FDP <- function(x.real, x.result, t) {
  if (sum(x.result[1:t]) != 0) {
    sum((1 - x.real[1:t]) * x.result[1:t]) / sum(x.result[1:t])
  } else 0
}
cal.MDP <- function(x.real, x.result, t) {
  if (sum(x.real[1:t]) != 0) {
    sum(x.real[1:t] * (1 - x.result[1:t])) / sum(x.real[1:t])
  } else 0
}
cal.dFDP <- function(x.real, x.result, t, d) {
  if (sum(x.result[1:t]) != 0) {
    w <- d^(t - seq_len(t))
    sum(w * (1 - x.real[1:t]) * x.result[1:t]) / sum(w * x.result[1:t])
  } else 0
}
cal.dMDP <- function(x.real, x.result, t, d) {
  if (sum(x.real[1:t]) != 0) {
    w <- d^(t - seq_len(t))
    sum(w * x.real[1:t] * (1 - x.result[1:t])) / sum(w * x.real[1:t])
  } else 0
}

memLORD <- function(d, delta, alpha = 0.1, w0 = alpha/10){
  N <- length(d)
  gammaJ <- rep(0,N)
  gammaJ[1] <- 0.07720838*log(2)
  for (J in seq(N-1)){
    gammaJ[J+1] <- 0.07720838*log(max(2,J+1))/((J+1)*exp(sqrt(log(J+1))))
  }
  tau <- 0
  b <- alpha-w0
  alphai <- rep(0,N)
  R <- rep(0,N)

  alphai[1] <- w0*gammaJ[1]
  R[1] <- round(d[1]<=alphai[1])
  if (R[1]>0){tau <- 1;b <- alpha}
  if (N==1){return(data.frame(pval=d,alphai=alphai,R=R))}

  for (J in seq(2,N)){
    if (length(tau)==1&&tau==0){alphai[J] <- gammaJ[J]*w0}
    if (length(tau)==1&&tau>0){alphai[J] <- gammaJ[J]*w0*delta^(J-tau)+(alpha-w0)*gammaJ[J-tau]*delta^(J-tau)}
    if (length(tau)>1){
      sog <- sum(delta^(J-tau[2:length(tau)])*gammaJ[J-tau[2:length(tau)]])
      alphai[J] <- gammaJ[J]*w0*delta^(J-tau[1])+(alpha-w0)*gammaJ[J-tau[1]]*delta^(J-tau[1])+alpha*sog
    }
    
    R[J] <- round(d[J]<=alphai[J])
    

    if (R[J]>0){tau <- which(R==1);b <- alpha}
  }
  out <- data.frame(pval=d,alphai=alphai,R=R)
  return(out)
}

one.exper <- function(
    seed,
    alpha = 0.01,
    lambda = 0.1,
    t = 10000,
    muc = 4,
    pi1 = 0.01,
    sigma = 1,
    rho = -1/(t-1)+1e-5,
    eta = 0.01,
    t0 = t/2,
    r = 100,
    reward.omega = 0.5,
    penalty.omega = 0.5,
    init.omega = 1/t,
    d = 0.99,
    w0 = alpha/10,
    epsilon.a = 1e-6,
    epsilon.r = 0.1
) {
  set.seed(seed)
  sigma <- 1
  
  out.loc <- sample(c(0,1), t, replace=TRUE, prob=c(1-pi1,pi1))
  X <- rep(0,t)
  phi <- rep(0,t)
  p <- rep(0,t)
  e <- rep(0,t)
  if (out.loc[1]==1){
    X[1] <- muc+rnorm(1,mean = 0,sd=sigma)
  }else{X[1] <- rnorm(1,mean = 0,sd=sigma)}
  p[1] <- 1-pnorm(X[1],mean = 0,sd=sigma)
  e[1] <- exp((muc*(2*X[1]-muc))/(2*sigma^2))
  
  for (j in seq(2,t)){
    phi[j] <- 2/(1+exp(-eta*(j-t0)))-1
    if (out.loc[j]==1){
      X[j] <- phi[j]*X[j-1]+muc+rnorm(1,mean = 0,sd=sigma)
    }else{X[j] <- phi[j]*X[j-1]+rnorm(1,mean = 0,sd=sigma)}
    
    p[j] <- 1-pnorm(X[j],mean = phi[j]*X[j-1],sd=sigma)
    e[j] <- exp((muc*(2*X[j]-muc-2*phi[j]*X[j-1]))/(2*sigma^2))
  }

  methods <- c(
    "mem_e_LORD", "mem_e_SAFFRON", "mem_LORDpp",
    "e_LORD", "e_SAFFRON",
    "mem_e_LORD_reset", "mem_e_SAFFRON_reset", "mem_LORDpp_reset",
    "mem_cp_LORD", "mem_cp_SAFFRON",
    "mem_cp_LORD_reset", "mem_cp_SAFFRON_reset"
  )
  rej_list  <- vector("list", length(methods))
  names(rej_list) <- methods

  init_vectors <- function() {
    list(
      gamma = numeric(t),
      level = numeric(t),
      rej   = integer(t),
      dcum  = 0,
      cum   = 0,
      rem   = 0,
      resets = 0
    )
  }
  
  # ------------------- mem-e-LORD -------------------
  m0 <- init_vectors()
  m0$gamma[1] <- init.omega
  m0$level[1] <- init.omega * alpha
  m0$rej[1]   <- as.integer(e[1] > 1 / m0$level[1])
  m0$dcum     <- m0$rej[1]
  m0$cum      <- m0$rej[1]
  m0$rem      <- alpha - m0$level[1]
  for (j in 2:t) {
    m0$gamma[j] <- m0$gamma[j-1] + init.omega * penalty.omega^(j-1 - m0$cum) * (1 - m0$rej[j-1]) -
      init.omega * reward.omega^m0$cum * m0$rej[j-1]
    m0$level[j] <- m0$gamma[j] * m0$rem * (d * m0$dcum + 1)
    if (is.na(m0$level[j])) m0$level[j] <- 0
    m0$rej[j]   <- as.integer(e[j] > 1 / m0$level[j])
    m0$rem      <- m0$rem - m0$level[j] / (d * m0$dcum + 1)
    m0$dcum     <- d * m0$dcum + m0$rej[j]
    m0$cum      <- m0$cum + m0$rej[j]
  }
  rej_list$mem_e_LORD <- m0$rej
  
  # ----------------- mem-e-SAFFRON ------------------
  m1 <- init_vectors()
  m1$gamma[1] <- init.omega
  m1$level[1] <- init.omega * alpha * (1 - lambda)
  m1$rej[1]   <- as.integer(e[1] > 1 / m1$level[1])
  m1$dcum     <- m1$rej[1]
  m1$cum      <- m1$rej[1]
  m1$rem      <- alpha*(1-lambda) - m1$level[1] * as.integer(e[1] < 1/lambda)
  for (j in 2:t) {
    m1$gamma[j] <- m1$gamma[j-1] + init.omega * penalty.omega^(j-1 - m1$cum) * (1 - m1$rej[j-1]) -
      init.omega * reward.omega^m1$cum * m1$rej[j-1]
    m1$level[j] <- m1$gamma[j] * m1$rem * (d * m1$dcum + 1)
    if (is.na(m1$level[j])) m1$level[j] <- 0
    m1$rej[j]   <- as.integer(e[j] > 1 / m1$level[j])
    m1$rem      <- m1$rem - as.integer(e[j] < 1/lambda) * m1$level[j] / (d * m1$dcum + 1)
    m1$dcum     <- d * m1$dcum + m1$rej[j]
    m1$cum      <- m1$cum + m1$rej[j]
  }
  rej_list$mem_e_SAFFRON <- m1$rej
  
  # ------------------ mem-LORD++ --------------------
  res2 <- memLORD(p, d, alpha)
  rej_list$mem_LORDpp <- res2$R
  
  # ------------------- e-LORD -----------------------
  e0 <- init_vectors()
  e0$gamma[1] <- init.omega
  e0$level[1] <- init.omega * alpha
  e0$rej[1]   <- as.integer(e[1] > 1 / e0$level[1])
  e0$cum      <- e0$rej[1]
  e0$rem      <- alpha - e0$level[1]
  for (j in 2:t) {
    e0$gamma[j] <- e0$gamma[j-1] + init.omega * penalty.omega^(j-1 - e0$cum) * (1 - e0$rej[j-1]) -
      init.omega * reward.omega^e0$cum * e0$rej[j-1]
    e0$level[j] <- e0$gamma[j] * e0$rem * (e0$cum + 1)
    if (is.na(e0$level[j])) e0$level[j] <- 0
    e0$rej[j]   <- as.integer(e[j] > 1 / e0$level[j])
    e0$rem      <- e0$rem - e0$level[j] / (e0$cum + 1)
    e0$cum      <- e0$cum + e0$rej[j]
  }
  rej_list$e_LORD <- e0$rej
  
  # ----------------- e-SAFFRON ----------------------
  e1 <- init_vectors()
  e1$gamma[1] <- init.omega
  e1$level[1] <- init.omega * alpha * (1 - lambda)
  e1$rej[1]   <- as.integer(e[1] > 1 / e1$level[1])
  e1$cum      <- e1$rej[1]
  e1$rem      <- alpha*(1-lambda) - e1$level[1] * as.integer(e[1] < 1/lambda)
  for (j in 2:t) {
    e1$gamma[j] <- e1$gamma[j-1] + init.omega * penalty.omega^(j-1 - e1$cum) * (1 - e1$rej[j-1]) -
      init.omega * reward.omega^e1$cum * e1$rej[j-1]
    e1$level[j] <- e1$gamma[j] * e1$rem * (e1$cum + 1)
    if (is.na(e1$level[j])) e1$level[j] <- 0
    e1$rej[j]   <- as.integer(e[j] > 1 / e1$level[j])
    e1$rem      <- e1$rem - as.integer(e[j] < 1/lambda) * e1$level[j] / (e1$cum + 1)
    e1$cum      <- e1$cum + e1$rej[j]
  }
  rej_list$e_SAFFRON <- e1$rej
  
  # ---------------- mem-e-LORD-reset ----------------
  mr0 <- m0  # start from memory-e-LORD state
  mr0$resets <- 0
  for (j in seq_len(t)) {
    if (mr0$level[j]/alpha < epsilon.a && mr0$dcum < epsilon.r) {
      # reset
      mr0$gamma[j] <- init.omega - init.omega * penalty.omega^j * (1 - mr0$rej[j]) - init.omega * reward.omega * mr0$rej[j]
      mr0$rem      <- alpha
      mr0$dcum     <- 0
      mr0$cum      <- 0
      mr0$resets   <- mr0$resets + 1
    }
  }
  rej_list$mem_e_LORD_reset <- mr0$rej
  
  # --------------- mem-e-SAFFRON-reset -------------
  mr1 <- m1
  mr1$resets <- 0
  for (j in seq_len(t)) {
    if (mr1$level[j]/alpha < epsilon.a && mr1$dcum < epsilon.r) {
      mr1$gamma[j] <- init.omega - init.omega * penalty.omega^j * (1 - mr1$rej[j]) - init.omega * reward.omega * mr1$rej[j]
      mr1$rem      <- alpha*(1-lambda)
      mr1$dcum     <- 0
      mr1$cum      <- 0
      mr1$resets   <- mr1$resets + 1
    }
  }
  rej_list$mem_e_SAFFRON_reset <- mr1$rej
  
  # ------------- mem-LORD++-reset ------------------
  rr <- res2
  rr$reset_count <- 0
  detect <- rr$alphai
  t0 <- t; offset <- 0
  while (any(detect/alpha < epsilon.a) && offset < t0) {
    bad <- which(detect/alpha < epsilon.a)[1]
    if (sum(d^(bad - seq_len(bad)) * rr$R[(offset+1):(offset+bad)]) < epsilon.r) {
      offset <- offset + bad
      rr2 <- memLORD(p[(offset+1):t0], d, alpha)
      rr$R[(offset+1):t0]     <- rr2$R
      rr$alphai[(offset+1):t0] <- rr2$alphai
      rr$reset_count <- rr$reset_count + 1
      detect <- rr$alphai[(offset+1):t0]
    } else break
  }
  rej_list$mem_LORDpp_reset <- rr$R
  
  # --------------- mem-cp-LORD ---------------------
  c0 <- init_vectors()
  c0$gamma[1] <- init.omega
  c0$level[1] <- init.omega * alpha
  c0$rej[1]   <- as.integer(p[1] < c0$level[1])
  c0$dcum     <- c0$rej[1]
  c0$cum      <- c0$rej[1]
  c0$rem      <- alpha - c0$level[1]
  for (j in 2:t) {
    c0$gamma[j] <- c0$gamma[j-1] + init.omega * penalty.omega^(j-1 - c0$cum) * (1 - c0$rej[j-1]) -
      init.omega * reward.omega^c0$cum * c0$rej[j-1]
    c0$level[j] <- c0$gamma[j] * c0$rem * (d * c0$dcum + 1)
    if (is.na(c0$level[j])) c0$level[j] <- 0
    c0$rej[j]   <- as.integer(p[j] < c0$level[j])
    c0$rem      <- c0$rem - c0$level[j] / (d * c0$dcum + 1)
    c0$dcum     <- d * c0$dcum + c0$rej[j]
    c0$cum      <- c0$cum + c0$rej[j]
  }
  rej_list$mem_cp_LORD <- c0$rej
  
  # ------------- mem-cp-SAFFRON --------------------
  c1 <- init_vectors()
  c1$gamma[1] <- init.omega
  c1$level[1] <- init.omega * alpha * (1 - lambda)
  c1$rej[1]   <- as.integer(p[1] < c1$level[1])
  c1$dcum     <- c1$rej[1]
  c1$cum      <- c1$rej[1]
  c1$rem      <- alpha*(1-lambda) - c1$level[1] * as.integer(p[1] > lambda)
  for (j in 2:t) {
    c1$gamma[j] <- c1$gamma[j-1] + init.omega * penalty.omega^(j-1 - c1$cum) * (1 - c1$rej[j-1]) -
      init.omega * reward.omega^c1$cum * c1$rej[j-1]
    c1$level[j] <- c1$gamma[j] * c1$rem * (d * c1$dcum + 1)
    if (is.na(c1$level[j])) c1$level[j] <- 0
    c1$rej[j]   <- as.integer(p[j] < c1$level[j])
    c1$rem      <- c1$rem - as.integer(p[j] > lambda) * c1$level[j] / (d * c1$dcum + 1)
    c1$dcum     <- d * c1$dcum + c1$rej[j]
    c1$cum      <- c1$cum + c1$rej[j]
  }
  rej_list$mem_cp_SAFFRON <- c1$rej
  
  # --------- mem-cp-LORD-reset --------------------
  cr0 <- c0; cr0$resets <- 0
  for (j in seq_len(t)) {
    if (cr0$level[j]/alpha < epsilon.a && cr0$dcum < epsilon.r) {
      cr0$gamma[j] <- init.omega - init.omega*penalty.omega^j*(1-cr0$rej[j]) - init.omega*reward.omega*cr0$rej[j]
      cr0$rem      <- alpha
      cr0$dcum     <- 0
      cr0$cum      <- 0
      cr0$resets   <- cr0$resets + 1
    }
  }
  rej_list$mem_cp_LORD_reset <- cr0$rej
  
  # ------- mem-cp-SAFFRON-reset ------------------
  cr1 <- c1; cr1$resets <- 0
  for (j in seq_len(t)) {
    if (cr1$level[j]/alpha < epsilon.a && cr1$dcum < epsilon.r) {
      cr1$gamma[j] <- init.omega - init.omega*penalty.omega^j*(1-cr1$rej[j]) - init.omega*reward.omega*cr1$rej[j]
      cr1$rem      <- alpha*(1-lambda)
      cr1$dcum     <- 0
      cr1$cum      <- 0
      cr1$resets   <- cr1$resets + 1
    }
  }
  rej_list$mem_cp_SAFFRON_reset <- cr1$rej

  FDPs <- sapply(methods, function(m) cal.dFDP(out.loc, rej_list[[m]], t, d))
  MDPs <- sapply(methods, function(m) cal.dMDP(out.loc, rej_list[[m]], t, d))

  output <- c(FDPs, MDPs)
  names(output) <- c(
    paste0("dFDP_", methods),
    paste0("dMDP_", methods)
  )
  return(output)
}

multi.exper <- function(    seed = 2025,
                            alpha = 0.05,
                            lambda = 0.1,
                            t = 20000,
                            muc = 3,
                            pi1 = 0.01,
                            sigma = 1,
                            rho = -1/(t-1) + 1e-5,
                            eta = 0.01,
                            t0 = t/2,
                            r = 100,
                            reward.omega = 0.5,
                            penalty.omega = 0.5,
                            init.omega = 1/t,
                            d = 0.99,
                            w0 = alpha/10,
                            epsilon.a = 1e-6,
                            epsilon.r = 0.1) {

  res_list <- lapply(seq_len(r), function(i) {
    tryCatch({
      one.exper(seed         = seed + i,
                alpha        = alpha,
                lambda       = lambda,
                t            = t,
                muc          = muc,
                pi1          = pi1,
                sigma        = sigma,
                rho          = rho,
                eta          = eta,
                t0           = t0,
                r            = r,
                reward.omega = reward.omega,
                penalty.omega= penalty.omega,
                init.omega   = init.omega,
                d            = d,
                w0           = w0,
                epsilon.a    = epsilon.a,
                epsilon.r    = epsilon.r)
    }, error = function(e) {
      message(sprintf("one.exper error %d experiment：%s", i, e$message))
      NULL
    })
  })
  

  res_list <- Filter(Negate(is.null), res_list)


  experiment_results <- do.call(cbind, res_list)
  
  avg <- rowMeans(experiment_results)
  
  se_idx <- grep("^(dFDP|dMDP)_", names(avg))
  
  n_eff <- ncol(experiment_results)
  se_vals <- apply(experiment_results[se_idx, , drop = FALSE], 1, sd) / sqrt(n_eff)
  names(se_vals) <- paste0("SE_", names(avg)[se_idx])
  
  c(avg, se_vals)
}


full_value_names <- col_names <- c(
  "dFDP_mem_e_LORD",
  "dFDP_mem_e_SAFFRON",
  "dFDP_mem_LORDpp",
  "dFDP_e_LORD",
  "dFDP_e_SAFFRON",
  "dFDP_mem_e_LORD_reset",
  "dFDP_mem_e_SAFFRON_reset",
  "dFDP_mem_LORDpp_reset",
  "dFDP_mem_cp_LORD",
  "dFDP_mem_cp_SAFFRON",
  "dFDP_mem_cp_LORD_reset",
  "dFDP_mem_cp_SAFFRON_reset",
  "dMDP_mem_e_LORD",
  "dMDP_mem_e_SAFFRON",
  "dMDP_mem_LORDpp",
  "dMDP_e_LORD",
  "dMDP_e_SAFFRON",
  "dMDP_mem_e_LORD_reset",
  "dMDP_mem_e_SAFFRON_reset",
  "dMDP_mem_LORDpp_reset",
  "dMDP_mem_cp_LORD",
  "dMDP_mem_cp_SAFFRON",
  "dMDP_mem_cp_LORD_reset",
  "dMDP_mem_cp_SAFFRON_reset",
  "SE_dFDP_mem_e_LORD",
  "SE_dFDP_mem_e_SAFFRON",
  "SE_dFDP_mem_LORDpp",
  "SE_dFDP_e_LORD",
  "SE_dFDP_e_SAFFRON",
  "SE_dFDP_mem_e_LORD_reset",
  "SE_dFDP_mem_e_SAFFRON_reset",
  "SE_dFDP_mem_LORDpp_reset",
  "SE_dFDP_mem_cp_LORD",
  "SE_dFDP_mem_cp_SAFFRON",
  "SE_dFDP_mem_cp_LORD_reset",
  "SE_dFDP_mem_cp_SAFFRON_reset",
  "SE_dMDP_mem_e_LORD",
  "SE_dMDP_mem_e_SAFFRON",
  "SE_dMDP_mem_LORDpp",
  "SE_dMDP_e_LORD",
  "SE_dMDP_e_SAFFRON",
  "SE_dMDP_mem_e_LORD_reset",
  "SE_dMDP_mem_e_SAFFRON_reset",
  "SE_dMDP_mem_LORDpp_reset",
  "SE_dMDP_mem_cp_LORD",
  "SE_dMDP_mem_cp_SAFFRON",
  "SE_dMDP_mem_cp_LORD_reset",
  "SE_dMDP_mem_cp_SAFFRON_reset"
)




