library(kbal)  ## it takes couple of hours
library(mvnfast)
library(WeightIt) ## 1.2.1
library(marginaleffects)
library(cobalt)

library(data.table)
library(laeken)

library(doSNOW)
library(progress)
library(foreach)
library(doRNG)

source("codes/functions.R")

cores <- 8
sims <- 2*70*cores ## 1000 iterations

n <- 2000 ## initial sample size
n_sample <- 1000 ## required sample size n_0=n_1=n (#n control = #n treatment)
mu <- c(0, 0, 0) 
Sigma <- matrix(c(2, 1, -1,
                  1, 1, -0.5,
                  -1,-0.5, 1), 
                nrow=3, 
                ncol=3)

tau <- c(0.10, 0.25, 0.5, 0.75, 0.9)
probs1 <- c(0.25, 0.5, 0.75)
probs2 <- seq(0.1, 0.9,0.1)
fm1 <- treatment ~ x1 + x2 + x3 + x4 + x5 + x6
fm2_y1 <- y1 ~ treatment
fm2_y2 <- y2 ~ treatment
fm2_y3 <- y3 ~ treatment

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

cl <- makeCluster(cores)
registerDoSNOW(cl)

registerDoRNG(2024)

results_simulation1 <- foreach(
  k=1:sims, 
  .combine = rbind,  
  .packages = c("mvnfast", "kbal", "WeightIt", "cobalt", "kbal", "data.table", "laeken", "marginaleffects"), 
  .options.snow = opts) %dopar% {
    
    x123 <- rmvn(n=n, mu=mu, sigma=Sigma)
    x1 <- x123[,1]
    x2 <- x123[,2]
    x3 <- x123[,3]
    x4 <- runif(n,-3,3)
    x5 <- rchisq(n,1)
    x6 <- rbinom(n,1,0.5)
    ep1 <- rnorm(n,0,30)
    ep2 <- rnorm(n,0,100)
    ep3 <- rchisq(n, 5)
    ep3 <- ((ep3 - mean(ep3)) / sd(ep3)) * sqrt(67.6) + 0.5
    d1 <- as.numeric((x1 + 2*x2 -2*x3 -x4 -0.5*x5 +x6 + ep1) > 0)
    d2 <- as.numeric((x1 + 2*x2 -2*x3 -x4 -0.5*x5 +x6 + ep2) > 0)
    d3 <- as.numeric((x1 + 2*x2 -2*x3 -x4 -0.5*x5 +x6 + ep3) > 0)
    eta <- rnorm(n)
    
    y1 <- x1 + x2 + x3 - x4 + x5 + x6 + eta
    y2 <- x1 + x2 + 0.2*x3*x4 - sqrt(x5) + eta
    y3 <- (x1 + x2 + x5)^2 + eta
    
    df <- data.frame(id =1:n, x1,x2,x3,x4,x5,x6,d1,d2,d3,y1,y2,y3) |> setDT()
    df_long <- melt(df, id.vars = c(1:7,11:13), value.name = "treatment", variable.name = "design")
    
    ## sample by 1000 from reach group
    df_long <- df_long[, .SD[sample(.N, n_sample, replace = T)], by = .(design, treatment)]
    y_results <- list()
    
    for (des in c("d1", "d2", "d3")) {
      
      df_long_des <- df_long[design == des]
      
      ks_rms <- bal.init(x = model.matrix(~x1+x2+x3+x4+x5+x6-1, data = df_long_des),
                     treat = df_long_des$treatment,
                     stat = "ks.rms",
                     estimand = "ATT")
      
      k_dist <- bal.init(x = model.matrix(~x1+x2+x3+x4+x5+x6-1, data = df_long_des),
                        treat = df_long_des$treatment,
                        stat = "kernel.dist",
                        estimand = "ATT")

      ebal_res <- weightit(formula = fm1,
                           data = df_long_des, estimand = "ATT", method = "ebal")
      
      kbal_res <- weightit(formula = fm1,
                           data = df_long_des, method = kbal.fun, estimand = "ATT",
                           include.obj = TRUE)
      
      deb_res1 <- weightit(formula = fm1,
                           data = df_long_des, estimand = "ATT", method = "ebal",
                           quantile = list(x1=probs1, x2=probs1,x3=probs1,x4=probs1,x5=probs1))
      
      deb_res2 <- weightit(formula = fm1,
                           data = df_long_des, estimand = "ATT", method = "ebal",
                           quantile = list(x1=probs2, x2=probs2,x3=probs2,x4=probs2,x5=probs2))
      
      if (any(deb_res2$weights == 0)) {
        print("some errors occured")
        next
      }
      
      df_long_des[, ":="(w_ebal = ebal_res$weights, 
                         w_kbal = kbal_res$weights,
                         w_deb1  = deb_res1$weights, 
                         w_deb2  = deb_res2$weights)]
      
      weightit_obj <- list(ebal=ebal_res, 
                           kbal=kbal_res, 
                           deb1=deb_res1, 
                           deb2=deb_res2)
      
      ebal_res_y1 <- lm_weightit(fm2_y1, data = df_long_des, weightit = ebal_res)
      ebal_res_y2 <- lm_weightit(fm2_y2, data = df_long_des, weightit = ebal_res)
      ebal_res_y3 <- lm_weightit(fm2_y3, data = df_long_des, weightit = ebal_res)
      
      kbal_res_y1 <- lm_weightit(fm2_y1, data = df_long_des, weightit = kbal_res)
      kbal_res_y2 <- lm_weightit(fm2_y2, data = df_long_des, weightit = kbal_res)
      kbal_res_y3 <- lm_weightit(fm2_y3, data = df_long_des, weightit = kbal_res)
      
      deb_res1_y1 <- lm_weightit(fm2_y1, data = df_long_des, weightit = deb_res1)
      deb_res1_y2 <- lm_weightit(fm2_y2, data = df_long_des, weightit = deb_res1)
      deb_res1_y3 <- lm_weightit(fm2_y3, data = df_long_des, weightit = deb_res1)
      
      deb_res2_y1 <- lm_weightit(fm2_y1, data = df_long_des, weightit = deb_res2)
      deb_res2_y2 <- lm_weightit(fm2_y2, data = df_long_des, weightit = deb_res2)
      deb_res2_y3 <- lm_weightit(fm2_y3, data = df_long_des, weightit = deb_res2)
      
      ## prediction with
      ## save into one data.frame
      weightit_results <- list(ebal_res_y1, ebal_res_y2, ebal_res_y3,
                               kbal_res_y1, kbal_res_y2, kbal_res_y3,
                               deb_res1_y1, deb_res1_y2, deb_res1_y3,
                               deb_res2_y1, deb_res2_y2, deb_res2_y3)
      
      ## actually i do not need this but leave as an option
      weightit_results_att <- lapply(weightit_results, avg_comparisons, variables = "treatment")
      
      coefs <- sapply(weightit_results_att, \(x) coef(x))
      cis <- lapply(weightit_results_att, \(x) c(x$conf.low, x$conf.high))
      cis <- do.call("rbind", cis)
      
      weightit_df <- data.table(treat = coefs, ci_low = cis[, 1], ci_upp = cis[, 2], est = "att", tau = NA)
      
      weightit_df[, method := c("ebal_y1", "ebal_y2", "ebal_y3",
                                "kbal_y1", "kbal_y2", "kbal_y3",
                                "deb1_y1", "deb1_y2", "deb1_y3",
                                "deb2_y1", "deb2_y2", "deb2_y3")]
      
      weightit_df[, c("method", "var"):=tstrsplit(method, "_")] 
      
      ## quantiles
      y1_res_q <- df_long_des[, lapply(.SD, 
                                       FUN = function(x) weightedQuantile(y1[treatment==1], x[treatment==1], tau) - 
                                         weightedQuantile(y1[treatment==0], x[treatment==0],tau)), 
                              .SDcols = patterns("w_")][, ":="(var="y1", est="qtt", tau=tau)]
      
      y2_res_q <- df_long_des[, lapply(.SD, 
                                       FUN = function(x) weightedQuantile(y2[treatment==1], x[treatment==1], tau) - 
                                         weightedQuantile(y2[treatment==0], x[treatment==0], tau)), 
                              .SDcols = patterns("w_")][, ":="(var="y2", est="qtt", tau=tau)]
      
      y3_res_q <- df_long_des[, lapply(.SD, 
                                       FUN = function(x) weightedQuantile(y3[treatment==1], x[treatment==1], tau) - 
                                         weightedQuantile(y3[treatment==0], x[treatment==0],tau)), 
                              .SDcols = patterns("w_")][, ":="(var="y3", est="qtt", tau=tau)]
      
      res_quantiles <- rbind(y1_res_q, y2_res_q, y3_res_q) |>
        melt(id.vars = c("var", "est", "tau"), value.name = "treat", variable.name = "method")
      
      
      ## quality information
      ## ks.adj
      weightit_quality <- lapply(weightit_obj, bal.tab, stats = c("m", "v", "ks"), thresholds = c(m = .05))
      weightit_bal <- lapply(weightit_quality, "[[", "Balance")
      weightit_bal <- lapply(weightit_bal, "[[", "KS.Adj")
      weightit_bal <- do.call("rbind", weightit_bal) |> as.data.table()
      names(weightit_bal) <- paste0("x", 1:6)
      weightit_bal[, method := names(weightit_obj)]
      weightit_bal[, est := "balance"]
      
      ## overall
      weightit_kernel <- sapply(weightit_obj, function(x) bal.compute(k_dist, weights = get.w(x)))
      weightit_ks <- sapply(weightit_obj, function(x) bal.compute(ks_rms, weights = get.w(x)))
      
      weightit_bal[, over_kernel := weightit_kernel]
      weightit_bal[, over_ks := weightit_ks]
      
      ## ESS
      weightit_ess <- sapply(weightit_quality, \(x) x$Observations$Control[2])
      weightit_bal[, ess := weightit_ess]
      
      y_results[[des]] <- rbind(weightit_df, res_quantiles, weightit_bal, fill = T)
      
    }
    
    rbindlist(y_results, idcol = "design")[, rep := k]
    
  }

stopCluster(cl)

saveRDS(results_simulation1, 
        file = "results/sim1-ebal-updated-1000.rds")

