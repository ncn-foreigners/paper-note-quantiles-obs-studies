## function taken from the paper Sant'Anna, P. H., Song, X., & Xu, Q. (2022). Covariate distribution balance via propensity scores. Journal of Applied Econometrics, 37(6), 1093-1120.
## in particular from here: http://qed.econ.queensu.ca/jae/datasets/santanna001/

cbps_sim_data <- function(n, dgp, scale.ps = 1){
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

