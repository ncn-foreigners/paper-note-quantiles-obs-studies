## function taken from the paper Sant'Anna, P. H., Song, X., & Xu, Q. (2022). Covariate distribution balance via propensity scores. Journal of Applied Econometrics, 37(6), 1093-1120.
## in particular from here: http://qed.econ.queensu.ca/jae/datasets/santanna001/

cbps_sim_data <- function(n, dgp, scale.ps = 1, y_nonlin = FALSE){
  #----------------------------------------------------------------------------
  #----------------------------------------------------------------------------
  # Kand and Schafer design - Correct model
  #----------------------------------------------------------------------------
  #----------------------------------------------------------------------------
  if (dgp==1){
    # Generate X
    X1 <- stats::rnorm(n)
    X2 <- stats::rnorm(n)
    X3 <- stats::rnorm(n)
    X4 <- stats::rnorm(n)
    X <- cbind(X1,X2,X3,X4)
    
    # Beta of pscores
    b.ps <- c(-1, 0.5, -0.25, -0.1)
    # Generate pscore and Treatment Status
    b.ps.Times.x <- X %*% b.ps
    pr <- 1/(1+exp(- scale.ps * b.ps.Times.x))
    Treat <- stats::rbinom(n, size=1, prob = pr)
    
    # the potential outcomes
    # Beta of regression
    b.reg <- c(27.4, 13.7, 13.7, 13.7)
    b.reg.Times.x <- X %*% b.reg
    if (y_nonlin) b.reg.Times.x <- (2*X[, 1] + 2*X[, 2] + X[, 3])^2
    
    Y0 <- 200 - b.reg.Times.x + stats::rnorm(n)
    Y1 <- 210 + b.reg.Times.x + stats::rnorm(n)
    # Observed outcome
    Y <- Y1*Treat + Y0*(1-Treat)
    # Dataset
    datamatrix <- data.frame(cbind(Y, Treat, X, X, pr))
    
    colnames(datamatrix) <- c("Y", "Treat", 
                              paste("X", 1:dim(X)[2], sep = ""), 
                              paste("trueX", 1:dim(X)[2], sep = ""),
                              "truePS")
  }
  
  #----------------------------------------------------------------------------
  # Kand and Schafer design - Misspecified model
  #----------------------------------------------------------------------------
  if (dgp == 2){
    # Generate X
    X1 <- stats::rnorm(n)
    X2 <- stats::rnorm(n)
    X3 <- stats::rnorm(n)
    X4 <- stats::rnorm(n)
    X <- cbind(X1,X2,X3,X4)
    
    # Beta of pscores
    b.ps <- c(-1, 0.5, -0.25, -0.1)
    # Generate pscore and Treatment Status
    b.ps.Times.x <- X %*% b.ps
    pr <- 1/(1+exp(- scale.ps * b.ps.Times.x))
    Treat <- stats::rbinom(n, size=1, prob = pr)
    
    # the potential outcomes
    # Beta of regression
    b.reg <- c(27.4, 13.7, 13.7, 13.7)
    b.reg.Times.x <- X %*% b.reg
    
    Y0 <- 200 - b.reg.Times.x + stats::rnorm(n)
    Y1 <- 210 + b.reg.Times.x + stats::rnorm(n)
    
    # Observed outcome
    Y <- Y1*Treat + Y0*(1-Treat)
    
    # Dataset
    W1 <- exp(X1/2)
    W2 <- X2/(1+exp(X1))
    W3 <- (X1*X3/25 + 0.6)^3 
    W4 <- (X1 + X4 + 20)^2 
    W <- cbind(W1, W2, W3, W4)
    datamatrix <- data.frame(cbind(Y, Treat, W, X, pr))
    
    colnames(datamatrix) <- c("Y", "Treat", 
                              paste("X", 1:dim(X)[2], sep = ""),
                              paste("trueX", 1:dim(X)[2], sep = ""),
                              "truePS")
  }
  #----------------------------------------------------------------------------
  
  #----------------------------------------------------------------------------
  
  return(datamatrix)
}


inflc_cbps<- function(Treat, X, pscore, weights = NULL, method = "exact"){
  Treat <- as.vector(Treat)
  n <- length(Treat)
  X <- as.matrix(X)
  k <- ncol(X)
  PS <- as.vector(pscore)
  if(is.null(weights)) weights <- rep(1, n)
  if(!is.numeric(weights)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # Avoid diving by zero
  probs.min<- 1e-6
  PS <- base::pmin(1-probs.min, PS)
  PS <- base::pmax(probs.min, PS)
  PS <- base::as.vector(PS)
  #-----------------------------------------------------------------------------
  if (method == "exact"){
    
    aux <- as.vector(weights * (Treat - PS)/(PS*(1-PS))) 
    g <- aux * X  # n x k 
    aux.dot <- as.vector(weights * (((Treat - PS)^2))/(PS*(1-PS)))
    G <- base::crossprod( aux.dot * X, X) / n 
    infl  <- g %*% MASS::ginv(G)   # n x k  k x k
    
  } else if (method == "over"){
    
    aux1 <- as.vector(weights * (Treat - PS)) 
    aux2 <- as.vector(weights * (Treat - PS)/(PS*(1-PS))) 
    g <- cbind( aux1 * X, aux2 * X) # n x 2k
    aux1.dot <- as.vector(weights * PS*(1 - PS)) 
    aux2.dot <- as.vector( - weights * (((Treat - PS)^2))/(PS*(1-PS)))
    G <- t(cbind(base::crossprod( aux1.dot * X, X), base::crossprod( aux2.dot * X, X))) / n  # 2k x k
    
    Omega <-  rbind(cbind(base::crossprod(weights *PS*(1 - PS)* X, X), 
                          base::crossprod(weights *X, X)),
                    cbind(base::crossprod(weights *X, X), 
                          base::crossprod(weights *X/(PS*(1 - PS)), X))) /n   # 2k x 2k
    
    infl <- -( g %*% MASS::ginv(Omega) %*% G) %*% MASS::ginv( t(G) %*% MASS::ginv(Omega) %*% G) 
    
  }
  
  return(infl)
}

dist_balance <- function(covs, weights, chunk = 1000){
  #-----------------------------------------------------------------------------
  # Ensure covs are stored as matrix
  covs <- as.matrix(covs)
  # Ensure weights are stored as numetic
  weights <- as.numeric(weights)
  #-----------------------------------------------------------------------------
  # Define variables to be used in the loop
  # Number of covariates
  k_dim = dim(covs)[2]
  # number of observations
  n_obs = dim(covs)[1]
  
  if(n_obs != length(weights)){
    stop("Weights must have same number of obs. as covs")
  }
  
  # Initialize `Rw` row vector (n_obs dimension)
  cdf_balance <- rep(0,n_obs)
  
  # We split n columns into l tiles, each with chunk columns
  l <- floor(n_obs/chunk) + 1
  
  #--------------------------------
  for (i in 1:l) {
    start <- min(chunk * (i - 1) + 1, n_obs)
    end <- min(chunk * i, n_obs)
    indicators <- base::outer(covs[,1], covs[start:end,1], "<=")
    if (k_dim>1) {
      for ( bb in 2:k_dim){
        indicators <- indicators * 
          base::outer(covs[, bb], covs[start:end, bb], "<=")
      }
    }
    
    cdf_balance[start:end] <- base::colMeans(weights * indicators)
  }
  
  ks <- sqrt(n_obs) * max(cdf_balance)
  cvm <- sum(cdf_balance^2)
  
  ret <- list(ks_bal = ks,
              cvm_bal = cvm)
  ret
}

