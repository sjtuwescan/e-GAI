pS <- function(d, alpha = 0.05, gammai, w0=0.0001, rho_reward = 0.5, rho_penalty = 0.5, lambda = 0.1, random = TRUE,
                     display_progress = FALSE, date.format = "%Y-%m-%d") {
  
  d <- checkPval(d)
  
  if (is.data.frame(d)) {
    d <- checkdf(d, random, date.format)
    pval <- d$pval
  } else if (is.vector(d)) {
    pval <- d
  } else {
    stop("d must either be a dataframe or a vector of p-values.")
  }
  
  N <- length(pval)
  
  if (alpha <= 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (lambda <= 0 || lambda > 1) {
    stop("lambda must be between 0 and 1.")
  }
  out <- pS_faster(pval, 
                         lambda = lambda,
                         alpha = alpha,
                         w0 = w0,
                         rho_reward = rho_reward, 
                         rho_penalty = rho_penalty,
                         display_progress = display_progress)
  out$R <- as.numeric(out$R)
  if(is.data.frame(d) && !is.null(d$id)) {
    out$id <- d$id
  }
  out
  
}

pS_faster <- function(pval, lambda, alpha, w0, rho_reward, rho_penalty, display_progress = TRUE) {
  N <- length(pval)
  # 初始化变量
  alphai <- numeric(N)
  wt <- numeric(N)
  R <- logical(N)
  
  wt[1] <- w0
  alphai[1] <- wt[1]*alpha
  R[1] <- (pval[1] <= alphai[1])
  
  Cjplus <- integer(N)
  cand <- integer(N)
  candsum <- 0
  K <- sum(R)
  
  for (i in 2:N) { # R 的下标从 1 开始，因此从 2 开始循环
    wt[i] <- wt[i-1] - w0*(R[i-1]*rho_reward^K - (1-R[i-1])*rho_penalty^(i-1-K))
    cand[i - 1] <- (pval[i-1] <= lambda)
    candsum <- candsum + (1-cand[i - 1])/(K-R[i-1]+1)*alphai[i-1]
    alphaitilde <- wt[i]*(alpha*(1-lambda)-candsum)*(K+1)
    alphai[i] <- min(lambda, alphaitilde)
    if (pval[i] <= alphai[i]) {
      R[i] <- TRUE
      K <- K + 1
    }
  }
  
  data.frame(pval = pval, alphai = alphai, R = R)
}
