## function taken from the paper Sant'Anna, P. H., Song, X., & Xu, Q. (2022). Covariate distribution balance via propensity scores. Journal of Applied Econometrics, 37(6), 1093-1120.
## in particular from here: http://qed.econ.queensu.ca/jae/datasets/santanna001/

dgps_ips <- function(n, dgp, scale.ps = 1){
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


##  https://stackoverflow.com/questions/58395772/r-data-table-sample-by-group-with-different-sampling-proportion
group_sampler <- function(data, group_col, sample_sizes) {
  data[, .SD[sample(.N, sample_sizes[.GRP])], keyby = group_col]
}


## for sim 1

eb_sim <- function(data, type, q_probs=probs_quar) {
  
  vars <- paste0("x", 1:6)
  vars_num <- vars[-6]
  df_sim <- copy(data[design == type])
  #df_sim <- group_sampler(data[design == type], "flag", c(control_n, treatment_n))
  df_sim_t <- df_sim[flag == TRUE] ## treated
  df_sim_c <- df_sim[flag == FALSE] ## control
  
  ## ebal
  invisible(
    capture.output(
      ebal_res <- ebalance(Treatment = df_sim$flag, X = df_sim[, ..vars])
    )
  )
  
  ## jointcal
  jcal_res_eb <- joint_calib(formula_quantiles = ~ x1 + x2 + x3 + x4 + x5,
                             formula_totals = ~ x1 + x2 + x3 + x4 + x5 + x6,
                             data = df_sim_c,
                             N = nrow(df_sim_t),
                             pop_quantiles = lapply(df_sim_t[, ..vars_num], quantile, probs = q_probs),
                             pop_totals = colSums(df_sim_t[, ..vars]),
                             method = "eb")
  
  jcal_res_raking <- joint_calib(formula_quantiles = ~ x1 + x2 + x3 + x4 + x5,
                                 formula_totals = ~ x1 + x2 + x3 + x4 + x5 + x6,
                                 data = df_sim_c,
                                 N = nrow(df_sim_t),
                                 pop_quantiles = lapply(df_sim_t[, ..vars_num], quantile, probs = q_probs),
                                 pop_totals = colSums(df_sim_t[, ..vars]),
                                 method = "raking")
  
  jcal_res_raking_qonly <- joint_calib(formula_quantiles = ~ x1 + x2 + x3 + x4 + x5,
                                       data = df_sim_c,
                                       N = nrow(df_sim_t),
                                       pop_quantiles = lapply(df_sim_t[, ..vars_num], quantile, probs = q_probs),
                                       method = "raking")
  
    df_result <- rbind(
      data.frame(est="eb", des = type,
                 df_sim_t[, lapply(.SD, mean), .SDcols = patterns("y")] - 
                   df_sim_c[, lapply(.SD, weighted.mean, w=ebal_res$w), .SDcols = patterns("y")]),
      data.frame(est="qrak", des = type,
                 df_sim_t[, lapply(.SD, mean), .SDcols = patterns("y")] - 
                   df_sim_c[, lapply(.SD, weighted.mean, w=jcal_res_raking$g), .SDcols = patterns("y")]),
      data.frame(est="qeb", des = type,
                 df_sim_t[, lapply(.SD, mean), .SDcols = patterns("y")] - 
                   df_sim_c[, lapply(.SD, weighted.mean, w=jcal_res_eb$g), .SDcols = patterns("y")]),
      data.frame(est="qrak_only", des = type,
                 df_sim_t[, lapply(.SD, mean), .SDcols = patterns("y")] - 
                   df_sim_c[, lapply(.SD, weighted.mean, w=jcal_res_raking_qonly$g), .SDcols = patterns("y")]),
      data.frame(est="naive", des = type,
                 df_sim_t[flag==TRUE, lapply(.SD, mean), .SDcols = patterns("y")] - 
                   df_sim_c[, lapply(.SD, mean), .SDcols = patterns("y")])
    )
  
  df_result
}

