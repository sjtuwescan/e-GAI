pL <- function(d, alpha = 0.05, gammai, w0=0.0001, rho_reward = 0.5, rho_penalty = 0.5, random = TRUE,
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
  
  out <- pL_faster(pval, 
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

pL_faster <- function(pval, alpha, w0, rho_reward, rho_penalty, display_progress = TRUE) {
  N <- length(pval)
  # 初始化变量
  alphai <- numeric(N)
  wt <- numeric(N)
  R <- logical(N)
  
  wt[1] <- w0
  alphai[1] <- wt[1]*alpha
  R[1] <- (pval[1] <= alphai[1])
  
  Cjplus <- integer(N)
  candsum <- 0
  K <- sum(R)
  
  
  for (i in 2:N){
    wt[i] <- wt[i-1] - w0*(R[i-1]*rho_reward^K - (1-R[i-1])*rho_penalty^(i-1-K))
    #wt[i] <- wt[i-1] - w0*(R[i-1]*(K+1)^{-1.6} - (1-R[i-1])*(i-K)^{-1.6})
    candsum <- candsum + 1/(K-R[i-1]+1)*alphai[i-1]
    alphaitilde <- wt[i]*(alpha-candsum)*(K+1)
    alphai[i] <- alphaitilde
    if (pval[i] <= alphai[i]) {
      R[i] <- TRUE
      K <- K + 1
    }
  }
  
  data.frame(pval = pval, alphai = alphai, R = R)
}
