library(WeightIt) ## 1.2.1
library(cobalt)
library(IPS) ## remotes::install_github("pedrohcgs/IPS")

library(data.table)
library(laeken)

library(doSNOW)
library(progress)
library(foreach)
library(doRNG)

source("codes/functions.R")

cores <- 8
sims <- 125*cores

n_sample <- 1000
probs1 <- c(0.25, 0.5, 0.75)
probs2 <- seq(0.1, 0.9,0.1)
tau <- c(0.10, 0.25, 0.5, 0.75, 0.9)
bw <- "nrd0" 

fm1 <- Treat ~ X1 + X2 + X3 + X4
fm2 <- Y ~ Treat

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

cl <- makeCluster(cores)
registerDoSNOW(cl)

registerDoRNG(2024)

results_simulation2 <- foreach(
  k=1:sims, 
  .combine = rbind,  
  .packages = c("IPS", "WeightIt", "cobalt", "data.table"), 
  .options.snow = opts) %dopar% {
    
  # data generation process 1 -----------------------------------------------
  data_sim1 <- cbps_sim_data(n = n_sample, dgp=1)
  xcov <- model.matrix(~X1 + X2 + X3 + X4, data_sim1)
  
  ## just identified
  ## CBPS just identified with weightit
  cbps_just <- CBPS::CBPS(formula = Treat ~ X1 + X2 + X3 + X4, data = data_sim1, ATT = 0, standardize = F, method = "exact")
  
  cbps_just1 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps")
  cbps_just1_fit <- lm_weightit(fm2, data = data_sim1, weightit = cbps_just1)
  cbps_just1_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, cbps_just1$ps, method = "exact")
  cbps_just1_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat,  x= xcov,  
                             ps = cbps_just1$ps, beta.lin.rep = cbps_just1_lin_pred, tau=tau, bw = bw)
  ## DPS with quartiles
  dps_just1 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps", quantile = list(probs1))
  dps_just1_fit <- lm_weightit(fm2, data = data_sim1, weightit = dps_just1)
  dps_just1_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, dps_just1$ps, method = "exact")
  dps_just1_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat,  x= xcov,  
                            ps = dps_just1$ps, beta.lin.rep = dps_just1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with deciles
  dps_just2 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps", over = F, quantile = list(probs2))
  dps_just2_fit <- lm_weightit(fm2, data = data_sim1, weightit = dps_just2)
  dps_just2_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, dps_just2$ps, method = "exact")
  dps_just2_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat, x= xcov,  
                            ps = dps_just2$ps, beta.lin.rep = dps_just2_lin_pred, tau=tau, bw = bw)
  
  ## CBPS over-identified
  cbps_over1 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps", over = T)
  cbps_over1_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, cbps_over1$ps, method = "over")
  cbps_over1_ate <- IPS::ATE(data_sim1$Y, data_sim1$Treat, xcov, cbps_over1$ps, cbps_over1_lin_pred)
  cbps_over1_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat,  x= xcov,  
                             ps = cbps_over1$ps, beta.lin.rep = cbps_over1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with quartiles
  dps_over1 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps", over = T, quantile = list(probs1))
  dps_over1_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, dps_over1$ps, method = "over")
  dps_over1_ate <- IPS::ATE(data_sim1$Y, data_sim1$Treat, xcov, dps_over1$ps, dps_over1_lin_pred)
  dps_over1_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat,  x= xcov,  
                             ps = dps_over1$ps, beta.lin.rep = dps_over1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with deciles
  dps_over2 <- weightit(formula = fm1, data = data_sim1, estimand = "ATE", method = "cbps", over = T, quantile = list(probs2))
  dps_over2_lin_pred <- inflc_cbps(data_sim1$Treat, xcov, dps_over2$ps, method = "over")
  dps_over2_ate <- IPS::ATE(data_sim1$Y, data_sim1$Treat, xcov, dps_over2$ps, dps_over2_lin_pred)
  dps_over2_qte <- IPS::QTE(y = data_sim1$Y, d = data_sim1$Treat,  x= xcov,  
                             ps = dps_over2$ps, beta.lin.rep = dps_over2_lin_pred, tau=tau, bw = bw)
  
  ## IPS ind
  fit.ind <- IPS::IPS_ind(d = data_sim1$Treat, x = xcov, maxit = 25000, beta.initial = cbps_just$coefficients)
  ate_IPS_ind <- IPS::ATE(data_sim1$Y,  data_sim1$Treat, xcov, fit.ind$fitted.values, fit.ind$lin.rep)
  qte_IPS_ind <- IPS::QTE(y=data_sim1$Y, d=data_sim1$Treat, x=xcov, 
                          ps=fit.ind$fitted.values, beta.lin.rep=fit.ind$lin.rep, tau=tau, bw = bw)
  ## IPS exp
  fit.exp <- IPS::IPS_exp(d = data_sim1$Treat, x = xcov, beta.initial = cbps_just$coefficients, maxit = 25000)
  ate_IPS_exp <- IPS::ATE(data_sim1$Y, data_sim1$Treat, xcov, fit.exp$fitted.values, fit.exp$lin.rep)
  qte_IPS_exp <- IPS::QTE(y=data_sim1$Y, d=data_sim1$Treat, x=xcov, 
                          ps=fit.exp$fitted.values,beta.lin.rep=fit.exp$lin.rep, tau=tau, bw = bw)
  
  ## IPS proj 
  fit.proj <- IPS::IPS_proj(d = data_sim1$Treat,x = xcov,beta.initial = cbps_just$coefficients, maxit = 25000)
  ate_IPS_proj <- IPS::ATE(data_sim1$Y, data_sim1$Treat, xcov, fit.proj$fitted.values, fit.proj$lin.rep)
  qte_IPS_proj <- IPS::QTE(y=data_sim1$Y, d=data_sim1$Treat, x=xcov, 
                           ps=fit.proj$fitted.values, beta.lin.rep=fit.proj$lin.rep, tau=tau, bw = bw)
  
  ## quality of weights
  ks_rms <- bal.init(x = data_sim1[,-1],
                     treat = data_sim1$Treat,
                     stat = "ks.rms",
                     estimand = "ATE")
  
  k_dist <- bal.init(x = data_sim1[,-1],
                     treat = data_sim1$Treat,
                     stat = "kernel.dist",
                     estimand = "ATE")
  
  ### quality as in Sant'anna et al. 
  for_quality_santa <- list(CBPS_j =cbps_just1$ps, DPS1_j =dps_just1$ps, DPS2_j =dps_just2$ps,
                            CBPS_o =cbps_over1$ps, DPS1_o =dps_over1$ps, DPS2_o = dps_over2$ps, 
                            IPS_ind = fit.ind$fitted.values, IPS_exp = fit.exp$fitted.values, IPS_proj = fit.proj$fitted.values)
  
  for_quality_santa_ww <- lapply(for_quality_santa, function(x) {
    w1 <- data_sim1$Treat/x
    w0 <- (1-data_sim1$Treat)/(1 - x)
    w1/mean(w1) - w0/mean(w0)
  })
  
  k_dim <- NCOL(xcov)
  
  indicators <- base::outer(xcov[,1], xcov[,1], "<=")
  for ( bb in 2:k_dim){
    indicators <- indicators *
      base::outer(xcov[, bb], xcov[, bb], "<=")
  }
  
  balance_cbps <- sapply(for_quality_santa_ww, \(x) abs(base::colMeans(x * indicators)))
  ks_tests <- sqrt(n_sample) * apply(balance_cbps, 2, max)
  cvm_tests <- colSums(balance_cbps^2)
  
  ### quality from weightit
  for_quality <- list(CBPS_j=cbps_just1,
                      DPS1_j=dps_just1,
                      DPS2_j=dps_just2,
                      CBPS_o=cbps_over1,
                      DPS1_o=dps_over1,
                      DPS2_o=dps_over2,
                      IPS_ind=as.weightit(x = data_sim1$Treat/fit.ind$fitted.values + (data_sim1$Treat==0)/(1 - fit.ind$fitted.values), 
                                          treat = data_sim1$Treat, estimand = "ATE", covs=xcov[,-1],
                                          ps = fit.ind$fitted.values),
                      IPS_exp=as.weightit(x = data_sim1$Treat/fit.exp$fitted.values + (data_sim1$Treat==0)/(1 - fit.exp$fitted.values), 
                                          treat = data_sim1$Treat, estimand = "ATE", covs=xcov[,-1],
                                          ps = fit.exp$fitted.values),
                      IPS_proj=as.weightit(x = data_sim1$Treat/fit.proj$fitted.values + (data_sim1$Treat==0)/(1 - fit.proj$fitted.values), 
                                           treat = data_sim1$Treat, estimand = "ATE", covs=xcov[,-1],
                                           ps = fit.proj$fitted.values))
  
  weightit_quality <- lapply(for_quality, bal.tab, stats = c("m", "v", "ks"), thresholds = c(m = .05))
  weightit_bal <- lapply(weightit_quality, "[[", "Balance")
  weightit_bal <- lapply(weightit_bal, "[[", "KS.Adj")
  weightit_bal <- lapply(weightit_bal, tail, n= 4)
  weightit_bal <- do.call("rbind", weightit_bal) |> as.data.table()
  names(weightit_bal) <- paste0("X", 1:4)
  weightit_bal[, method := names(for_quality)]
  weightit_bal[, est := "balance"]
  
  ## overall
  weightit_kernel <- sapply(for_quality, function(x) bal.compute(k_dist, weights = get.w(x)))
  weightit_ks <- sapply(for_quality, function(x) bal.compute(ks_rms, weights = get.w(x)))
  
  weightit_bal[, over_kernel := weightit_kernel]
  weightit_bal[, over_ks := weightit_ks]
  weightit_bal[, sant_ks := ks_tests]
  weightit_bal[, sant_cvm := cvm_tests]
  
  ## ESS
  weightit_ess <- sapply(weightit_quality, \(x) x$Observations$Control[2])
  weightit_bal[, ess := weightit_ess]
  weightit_bal[, d:=1]
  
  ## save results
  ## ate and se (do poprawy)
  ate_results1 <- rbind(
    data.frame(method = "CBPS_j",  d=1,tau=NA, eff = coef(summary(cbps_just1_fit))[2,1], eff_se = coef(summary(cbps_just1_fit))[2,2]),
    data.frame(method = "CBPS_o",  d=1,tau=NA, eff = cbps_over1_ate$ate, eff_se = cbps_over1_ate$ate.se[1,1]),
    data.frame(method = "DPS1_j",  d=1,tau=NA, eff = coef(summary(dps_just1_fit))[2,1],  eff_se = coef(summary(dps_just1_fit))[2,2]),
    data.frame(method = "DPS1_o",  d=1,tau=NA, eff = dps_over1_ate$ate,  eff_se = dps_over1_ate$ate.se[1,1]),
    data.frame(method = "DPS2_j",  d=1,tau=NA, eff = coef(summary(dps_just2_fit))[2,1],  eff_se = coef(summary(dps_just2_fit))[2,2]),
    data.frame(method = "DPS2_o",  d=1,tau=NA, eff = dps_over2_ate$ate,  eff_se = dps_over2_ate$ate.se[1,1]),
    data.frame(method = "IPS_ind", d=1,tau=NA, eff = ate_IPS_ind$ate,    eff_se = ate_IPS_ind$ate.se[1,1]),
    data.frame(method = "IPS_exp", d=1,tau=NA, eff = ate_IPS_exp$ate,    eff_se = ate_IPS_exp$ate.se[1,1]),
    data.frame(method = "IPS_proj",d=1,tau=NA, eff = ate_IPS_proj$ate,   eff_se = ate_IPS_proj$ate.se[1,1])
  )
  ## qte and se
  qte_results1 <- rbind(
    data.frame(method = "CBPS_j",  d=1,tau=tau, eff = cbps_just1_qte$qte, eff_se = cbps_just1_qte$qte.se),
    data.frame(method = "CBPS_o",  d=1,tau=tau, eff = cbps_over1_qte$qte, eff_se = cbps_over1_qte$qte.se),
    data.frame(method = "DPS1_j",  d=1,tau=tau, eff = dps_just1_qte$qte,  eff_se = dps_just1_qte$qte.se),
    data.frame(method = "DPS1_o",  d=1,tau=tau, eff = dps_over1_qte$qte,  eff_se = dps_over1_qte$qte.se),
    data.frame(method = "DPS2_j",  d=1,tau=tau, eff = dps_just2_qte$qte,  eff_se = dps_just2_qte$qte.se),
    data.frame(method = "DPS2_o",  d=1,tau=tau, eff = dps_over2_qte$qte,  eff_se = dps_over2_qte$qte.se),
    data.frame(method = "IPS_ind", d=1,tau=tau, eff = qte_IPS_ind$qte,    eff_se = qte_IPS_ind$qte.se),
    data.frame(method = "IPS_exp", d=1,tau=tau, eff = qte_IPS_exp$qte,    eff_se = qte_IPS_exp$qte.se),
    data.frame(method = "IPS_proj",d=1,tau=tau, eff = qte_IPS_proj$qte,   eff_se = qte_IPS_proj$qte.se)
  )
  
  result_d1 <- rbind(as.data.table(ate_results1), as.data.table(qte_results1), weightit_bal, fill = T)
  

  # data generation process 2 --------------------------------------------------------------

  data_sim2 <- cbps_sim_data(n = n_sample, dgp=2)
  xcov <- model.matrix(~X1 + X2 + X3 + X4, data_sim2)
  
  ## just identified
  ## CBPS just identified with weightit
  cbps_just <- CBPS::CBPS(formula = Treat ~ X1 + X2 + X3 + X4, data = data_sim2, ATT = 0, standardize = F, method = "exact")
  
  cbps_just1 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps")
  cbps_just1_fit <- lm_weightit(fm2, data = data_sim2, weightit = cbps_just1)
  cbps_just1_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, cbps_just1$ps, method = "exact")
  cbps_just1_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat,  x= xcov,  
                             ps = cbps_just1$ps, beta.lin.rep = cbps_just1_lin_pred, tau=tau, bw = bw)
  ## DPS with quartiles
  dps_just1 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps", quantile = list(probs1))
  dps_just1_fit <- lm_weightit(fm2, data = data_sim2, weightit = dps_just1)
  dps_just1_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, dps_just1$ps, method = "exact")
  dps_just1_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat,  x= xcov,  
                            ps = dps_just1$ps, beta.lin.rep = dps_just1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with deciles
  dps_just2 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps", over = F, quantile = list(probs2))
  dps_just2_fit <- lm_weightit(fm2, data = data_sim2, weightit = dps_just2)
  dps_just2_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, dps_just2$ps, method = "exact")
  dps_just2_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat, x= xcov,  
                            ps = dps_just2$ps, beta.lin.rep = dps_just2_lin_pred, tau=tau, bw = bw)
  
  ## CBPS over-identified
  cbps_over1 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps", over = T)
  cbps_over1_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, cbps_over1$ps, method = "over")
  cbps_over1_ate <- IPS::ATE(data_sim2$Y, data_sim2$Treat, xcov, cbps_over1$ps, cbps_over1_lin_pred)
  cbps_over1_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat,  x= xcov,  
                             ps = cbps_over1$ps, beta.lin.rep = cbps_over1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with quartiles
  dps_over1 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps", over = T, quantile = list(probs1))
  dps_over1_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, dps_over1$ps, method = "over")
  dps_over1_ate <- IPS::ATE(data_sim2$Y, data_sim2$Treat, xcov, dps_over1$ps, dps_over1_lin_pred)
  dps_over1_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat,  x= xcov,  
                            ps = dps_over1$ps, beta.lin.rep = dps_over1_lin_pred, tau=tau, bw = bw)
  
  ## DPS with deciles
  dps_over2 <- weightit(formula = fm1, data = data_sim2, estimand = "ATE", method = "cbps", over = T, quantile = list(probs2))
  dps_over2_lin_pred <- inflc_cbps(data_sim2$Treat, xcov, dps_over2$ps, method = "over")
  dps_over2_ate <- IPS::ATE(data_sim2$Y, data_sim2$Treat, xcov, dps_over2$ps, dps_over2_lin_pred)
  dps_over2_qte <- IPS::QTE(y = data_sim2$Y, d = data_sim2$Treat,  x= xcov,  
                            ps = dps_over2$ps, beta.lin.rep = dps_over2_lin_pred, tau=tau, bw = bw)
  
  ## IPS ind
  fit.ind <- IPS::IPS_ind(d = data_sim2$Treat, x = xcov, maxit = 25000, beta.initial = cbps_just$coefficients)
  ate_IPS_ind <- IPS::ATE(data_sim2$Y,  data_sim2$Treat, xcov, fit.ind$fitted.values, fit.ind$lin.rep)
  qte_IPS_ind <- IPS::QTE(y=data_sim2$Y, d=data_sim2$Treat, x=xcov, 
                          ps=fit.ind$fitted.values, beta.lin.rep=fit.ind$lin.rep, tau=tau, bw = bw)
  ## IPS exp
  fit.exp <- IPS::IPS_exp(d = data_sim2$Treat, x = xcov, beta.initial = cbps_just$coefficients, maxit = 25000)
  ate_IPS_exp <- IPS::ATE(data_sim2$Y, data_sim2$Treat, xcov, fit.exp$fitted.values, fit.exp$lin.rep)
  qte_IPS_exp <- IPS::QTE(y=data_sim2$Y, d=data_sim2$Treat, x=xcov, 
                          ps=fit.exp$fitted.values,beta.lin.rep=fit.exp$lin.rep, tau=tau, bw = bw)
  
  ## IPS proj 
  fit.proj <- IPS::IPS_proj(d = data_sim2$Treat,x = xcov,beta.initial = cbps_just$coefficients, maxit = 25000)
  ate_IPS_proj <- IPS::ATE(data_sim2$Y, data_sim2$Treat, xcov, fit.proj$fitted.values, fit.proj$lin.rep)
  qte_IPS_proj <- IPS::QTE(y=data_sim2$Y, d=data_sim2$Treat, x=xcov, 
                           ps=fit.proj$fitted.values, beta.lin.rep=fit.proj$lin.rep, tau=tau, bw = bw)
  
  ## quality of weights
  ks_rms <- bal.init(x = data_sim2[,-1],
                     treat = data_sim2$Treat,
                     stat = "ks.rms",
                     estimand = "ATE")
  
  k_dist <- bal.init(x = data_sim2[,-1],
                     treat = data_sim2$Treat,
                     stat = "kernel.dist",
                     estimand = "ATE")
  
  ### quality as in Sant'anna et al. 
  for_quality_santa <- list(CBPS_j =cbps_just1$ps, DPS1_j =dps_just1$ps, DPS2_j =dps_just2$ps,
                            CBPS_o =cbps_over1$ps, DPS1_o =dps_over1$ps, DPS2_o = dps_over2$ps, 
                            IPS_ind = fit.ind$fitted.values, IPS_exp = fit.exp$fitted.values, IPS_proj = fit.proj$fitted.values)
  
  for_quality_santa_ww <- lapply(for_quality_santa, function(x) {
    w1 <- data_sim2$Treat/x
    w0 <- (1-data_sim2$Treat)/(1 - x)
    w1/mean(w1) - w0/mean(w0)
  })
  
  k_dim <- NCOL(xcov)
  
  indicators <- base::outer(xcov[,1], xcov[,1], "<=")
  for ( bb in 2:k_dim){
    indicators <- indicators *
      base::outer(xcov[, bb], xcov[, bb], "<=")
  }
  
  balance_cbps <- sapply(for_quality_santa_ww, \(x) abs(base::colMeans(x * indicators)))
  ks_tests <- sqrt(n_sample) * apply(balance_cbps, 2, max)
  cvm_tests <- colSums(balance_cbps^2)
  
  ### quality from weightit
  
  for_quality <- list(CBPS_j=cbps_just1,
                      DPS1_j=dps_just1,
                      DPS2_j=dps_just2,
                      CBPS_o = cbps_over1,
                      DPS1_o=dps_over1,
                      DPS2_o=dps_over2,
                      IPS_ind=as.weightit(x = data_sim2$Treat/fit.ind$fitted.values + (data_sim2$Treat==0)/(1 - fit.ind$fitted.values), 
                                          treat = data_sim2$Treat, estimand = "ATE", covs=xcov[,-1],
                                          ps = fit.ind$fitted.values),
                      IPS_exp=as.weightit(x = data_sim2$Treat/fit.exp$fitted.values + (data_sim2$Treat==0)/(1 - fit.exp$fitted.values), 
                                          treat = data_sim2$Treat, estimand = "ATE", covs=xcov[,-1],
                                          ps = fit.exp$fitted.values),
                      IPS_proj=as.weightit(x = data_sim2$Treat/fit.proj$fitted.values + (data_sim2$Treat==0)/(1 - fit.proj$fitted.values), 
                                           treat = data_sim2$Treat, estimand = "ATE", covs=xcov[,-1],
                                           ps = fit.proj$fitted.values))
  
  weightit_quality <- lapply(for_quality, bal.tab, stats = c("m", "v", "ks"), thresholds = c(m = .05))
  weightit_bal <- lapply(weightit_quality, "[[", "Balance")
  weightit_bal <- lapply(weightit_bal, "[[", "KS.Adj")
  weightit_bal <- lapply(weightit_bal, tail, n= 4)
  weightit_bal <- do.call("rbind", weightit_bal) |> as.data.table()
  names(weightit_bal) <- paste0("X", 1:4)
  weightit_bal[, method := names(for_quality)]
  weightit_bal[, est := "balance"]
  
  ## overall
  weightit_kernel <- sapply(for_quality, function(x) bal.compute(k_dist, weights = get.w(x)))
  weightit_ks <- sapply(for_quality, function(x) bal.compute(ks_rms, weights = get.w(x)))
  
  weightit_bal[, over_kernel := weightit_kernel]
  weightit_bal[, over_ks := weightit_ks]
  weightit_bal[, sant_ks := ks_tests]
  weightit_bal[, sant_cvm := cvm_tests]
  
  ## ESS
  weightit_ess <- sapply(weightit_quality, \(x) x$Observations$Control[2])
  weightit_bal[, ess := weightit_ess]
  weightit_bal[, d:=2]
  
  ## save results
  ## ate and se (do poprawy)
  ate_results2 <- rbind(
    data.frame(method = "CBPS_j",  d=2,tau=NA, eff = coef(summary(cbps_just1_fit))[2,1], eff_se = coef(summary(cbps_just1_fit))[2,2]),
    data.frame(method = "CBPS_o",  d=2,tau=NA, eff = cbps_over1_ate$ate, eff_se = cbps_over1_ate$ate.se[1,1]),
    data.frame(method = "DPS1_j",  d=2,tau=NA, eff = coef(summary(dps_just1_fit))[2,1],  eff_se = coef(summary(dps_just1_fit))[2,2]),
    data.frame(method = "DPS1_o",  d=2,tau=NA, eff = dps_over1_ate$ate,  eff_se = dps_over1_ate$ate.se[1,1]),
    data.frame(method = "DPS2_j",  d=2,tau=NA, eff = coef(summary(dps_just2_fit))[2,1],  eff_se = coef(summary(dps_just2_fit))[2,2]),
    data.frame(method = "DPS2_o",  d=2,tau=NA, eff = dps_over2_ate$ate,  eff_se = dps_over2_ate$ate.se[1,1]),
    data.frame(method = "IPS_ind", d=2,tau=NA, eff = ate_IPS_ind$ate,    eff_se = ate_IPS_ind$ate.se[1,1]),
    data.frame(method = "IPS_exp", d=2,tau=NA, eff = ate_IPS_exp$ate,    eff_se = ate_IPS_exp$ate.se[1,1]),
    data.frame(method = "IPS_proj",d=2,tau=NA, eff = ate_IPS_proj$ate,   eff_se = ate_IPS_proj$ate.se[1,1])
  )
  ## qte and se
  qte_results2 <- rbind(
    data.frame(method = "CBPS_j",  d=2,tau=tau, eff = cbps_just1_qte$qte, eff_se = cbps_just1_qte$qte.se),
    data.frame(method = "CBPS_o",  d=2,tau=tau, eff = cbps_over1_qte$qte, eff_se = cbps_over1_qte$qte.se),
    data.frame(method = "DPS1_j",  d=2,tau=tau, eff = dps_just1_qte$qte,  eff_se = dps_just1_qte$qte.se),
    data.frame(method = "DPS1_o",  d=2,tau=tau, eff = dps_over1_qte$qte,  eff_se = dps_over1_qte$qte.se),
    data.frame(method = "DPS2_j",  d=2,tau=tau, eff = dps_just2_qte$qte,  eff_se = dps_just2_qte$qte.se),
    data.frame(method = "DPS2_o",  d=2,tau=tau, eff = dps_over2_qte$qte,  eff_se = dps_over2_qte$qte.se),
    data.frame(method = "IPS_ind", d=2,tau=tau, eff = qte_IPS_ind$qte,    eff_se = qte_IPS_ind$qte.se),
    data.frame(method = "IPS_exp", d=2,tau=tau, eff = qte_IPS_exp$qte,    eff_se = qte_IPS_exp$qte.se),
    data.frame(method = "IPS_proj",d=2,tau=tau, eff = qte_IPS_proj$qte,   eff_se = qte_IPS_proj$qte.se)
  )
  
  result_d2 <- rbind(as.data.table(ate_results2), as.data.table(qte_results2), weightit_bal, fill = T)
  
  result <- rbind(result_d1, result_d2)[, rep := k]
  result

}

stopCluster(cl)

saveRDS(results_simulation2, 
        file = "results/sim2-cbps-1000.rds")

