## setup
install.packages(c("remotes", "CBPS", "cobalt", "xtable"))
remotes::install_github("pedrohcgs/IPS") ## integrated CBPS
remotes::install_github("ngreifer/WeightIt") ## WeightIt development package version: 0.14.2.9002
remotes::install_github("ncn-foreigners/jointCalib")  ## 0.1.2 that contains joint_calib_cbps function

## cran packages
library("CBPS")
library("cobalt")
library("xtable")
## github packages
library("WeightIt")
library("jointCalib")
library("IPS")

source("codes/functions.R")

load("data/data401k.RData")

data401k <- subset(data401k, inc>0)
data401k <- within(data401k, {
  age <- (age - 25) / (64 - 25) 
  inc <- (inc + 2652) / (242124 + 2652)
  fsize <-  fsize / 13
  educ  <-  educ / 18
  loginc <- log(inc) 
  incsq <- inc^2
  logincsq <- loginc^2
  agesq <- age^2
  educsq <- educ^2
  fsizesq <- fsize^2
})
                      

fm <- e401 ~ (inc + loginc + age +fsize + educ+ hown + marr+ twoearn+ db+ pira)^2 + 
            incsq + logincsq + agesq + fsizesq + educsq

fm_y_bal <- e401 ~ inc + loginc + age + fsize + educ + hown + marr + twoearn +  db + pira

fm_x <- ~ (inc + loginc + age + fsize + educ + hown + marr + twoearn +  db + pira)^2 + 
          incsq + logincsq + agesq + fsizesq + educsq

fm_c <- ~ inc + loginc + age + fsize + educ

## quantiles
tau_v <- c(0.1, 0.25, 0.5, 0.75, 0.9)

## projection results
ips_proj <- readRDS("dataips-projection.Rds")

# CBPS --------------------------------------------------------------------

cbop0 <- CBPS(fm, 
              ATT = 0,
              standardize = F, 
              data = data401k, 
              twostep = T, 
              method = "exact")

res_cbop1 <- weightit(formula = fm, 
                      data = data401k,
                      method = "cbps", 
                      estimand = "ATE", 
                      over = FALSE,
                      twostep = T)

data401k$www_cbps1 <- res_cbop1$weights

res_cbop2 <- joint_calib_cbps(formula_quantiles =  fm_c,  
                              formula_means = fm_x,
                              treatment = ~ e401, 
                              data = data401k,
                              #probs = c(0.25, 0.5, 0.75),
                              probs = list(inc = seq(0.1,0.9,0.1),
                                           loginc = seq(0.1,0.9,0.1),
                                           age = seq(0.1,0.9,0.1),
                                           educ = c(0.25, 0.5, 0.75, 0.9),
                                           fsize = c(0.25, 0.5, 0.75, 0.9)),
                              standardize = F,
                              twostep = T ,
                              method = "exact")

data401k$www_cbps2  <- ifelse(data401k$e401==1, 1/res_cbop2$fitted.values, 1/(1-res_cbop2$fitted.values))
data401k$www_proj  <- ifelse(data401k$e401==1, 1/ips_proj$fitted.values, 1/(1-ips_proj$fitted.values))

tab_cbps1 <- bal.tab(x = fm_y_bal,
                     data = data401k,
                     weights = "www_cbps1",
                     s.d.denom = "pooled",
                     stats = c("m", "v", "ks"))

tab_cbps2 <- bal.tab(x = fm_y_bal,
                     data = data401k,
                     weights = "www_cbps2",
                     s.d.denom = "pooled",
                     stats = c("m", "v", "ks"))

tab_cbps3 <- bal.tab(x = fm_y_bal,
                     data = data401k,
                     weights = "www_proj",
                     s.d.denom = "pooled",
                     stats = c("m", "v", "ks"))

# average treatment effect ------------------------------------------------
## net_tfa
ate_cbps_exact <- IPS::ATE(data401k$net_tfa, 
                           cbop0$y,
                           cbop0$x, 
                           cbop0$fitted.values,  
                           inflc_cbps(cbop0$y, 
                                      cbop0$x, 
                                      cbop0$fitted.values, 
                                      method = "exact"))
ate_cbps_prop <- IPS::ATE(data401k$net_tfa, 
                          res_cbop2$y,
                          res_cbop2$x, 
                          res_cbop2$fitted.values,  
                          inflc_cbps(res_cbop2$y, res_cbop2$x, res_cbop2$fitted.values, method = "exact"))

ate_ips_proj <- IPS::ATE(data401k$net_tfa, 
                         cbop0$y,
                         cbop0$x, 
                         ips_proj$fitted.values, 
                         ips_proj$lin.rep)




# quantile treatment effect -----------------------------------------------

qte_cbps_exact <- IPS::QTE(data401k$net_tfa, 
                           d = cbop0$y,
                           x = cbop0$x, 
                           ps = cbop0$fitted.values,
                           inflc_cbps(cbop0$y,
                                      cbop0$x,
                                      cbop0$fitted.values,
                                      method = "exact"), 
                           tau = tau_v, bw = "nrd0", whs = NULL)

qte_cbps_prop <- IPS::QTE(data401k$net_tfa, 
                          d = res_cbop2$y,
                          x = res_cbop2$x, 
                          ps = res_cbop2$fitted.values,
                          inflc_cbps(res_cbop2$y,
                                     res_cbop2$x,
                                     res_cbop2$fitted.values,
                                     method = "exact"), 
                          tau = tau_v, bw = "nrd0", whs = NULL)

qte_ips_proj <- IPS::QTE(data401k$net_tfa,
                         d = cbop0$y,
                         x = cbop0$x,
                         ps = ips_proj$fitted.values,
                         beta.lin.rep = ips_proj$lin.rep,
                         tau = tau_v, bw = "nrd0", whs = NULL)



## tw
ate_cbps_exact_tw <- IPS::ATE(data401k$tw, 
                             cbop0$y,
                             cbop0$x, 
                             cbop0$fitted.values,  
                             inflc_cbps(cbop0$y, 
                                        cbop0$x, 
                                        cbop0$fitted.values, 
                                        method = "exact"))
ate_cbps_prop_tw <- IPS::ATE(data401k$tw, 
                      res_cbop2$y,
                      res_cbop2$x, 
                      res_cbop2$fitted.values,  
                      inflc_cbps(res_cbop2$y, res_cbop2$x, res_cbop2$fitted.values, method = "exact"))

ate_ips_proj_tw <- IPS::ATE(data401k$tw, 
                     cbop0$y,
                     cbop0$x, 
                     ips_proj$fitted.values, 
                     ips_proj$lin.rep)


rbind(unlist(ate_cbps_exact_tw[1:2]),
      unlist(ate_ips_proj_tw[1:2]),
      unlist(ate_cbps_prop_tw[1:2]))
      

# quantile treatment effect -----------------------------------------------

qte_cbps_exact_tw <- IPS::QTE(data401k$tw, 
                         d = cbop0$y,
                         x = cbop0$x, 
                         ps = cbop0$fitted.values,
                         inflc_cbps(cbop0$y,
                                    cbop0$x,
                                    cbop0$fitted.values,
                                    method = "exact"), 
                         tau = tau_v, bw = "nrd0", whs = NULL)

qte_cbps_prop_tw <- IPS::QTE(data401k$tw, 
                         d = res_cbop2$y,
                         x = res_cbop2$x, 
                         ps = res_cbop2$fitted.values,
                         inflc_cbps(res_cbop2$y,
                                    res_cbop2$x,
                                    res_cbop2$fitted.values,
                                    method = "exact"), 
                         tau = tau_v, bw = "nrd0", whs = NULL)

qte_ips_proj_tw <- IPS::QTE(data401k$tw,
                       d = cbop0$y,
                       x = cbop0$x,
                       ps = ips_proj$fitted.values,
                       beta.lin.rep = ips_proj$lin.rep,
                       tau = tau_v, bw = "nrd0", whs = NULL)

data.frame(
  cbps = qte_cbps_exact_tw$qte,
  proj = qte_ips_proj_tw$qte,
  dps = qte_cbps_prop_tw$qte)

data.frame(
  cbps = qte_cbps_exact_tw$qte.se,
  proj = qte_ips_proj_tw$qte.se,
  dps = qte_cbps_prop_tw$qte.se)


# balance checks ----------------------------------------------------------

xbal_names <- c("inc", "loginc" , "age" ,"fsize" , "educ", "hown" , "marr", "twoearn", "db", "pira")
xbal <- cbop0$x[,xbal_names]
n <- nrow(xbal)

w_cbps1 <- cbop0$y/cbop0$fitted.values
w_cbps0 <- (1-cbop0$y)/(1 - cbop0$fitted.values)
w_cbps1 <- w_cbps1/mean(w_cbps1)
w_cbps0 <- w_cbps0/mean(w_cbps0)
w_cbps_ate <- w_cbps1 - w_cbps0

w_dps1 <- cbop0$y/res_cbop2$fitted.values
w_dps0 <- (1-cbop0$y)/(1 - res_cbop2$fitted.values)
w_dps1 <- w_dps1/mean(w_dps1)
w_dps0 <- w_dps0/mean(w_dps0)
w_dps_ate <- w_dps1 - w_dps0

w_proj1 <- cbop0$y/ips_proj$fitted.values
w_proj0 <- (1-cbop0$y)/(1 - ips_proj$fitted.values)
w_proj1 <- w_proj1/mean(w_proj1)
w_proj0 <- w_proj0/mean(w_proj0)
w_proj_ate <- w_proj1 - w_proj0

k_dim <- ncol(xbal)

indicators <- base::outer(xbal[,1], xbal[,1], "<=")
for ( bb in 2:k_dim){
  indicators <- indicators * 
    base::outer(xbal[, bb], xbal[, bb], "<=")
}

cdf_balance_cbps <- abs(base::colMeans(w_cbps_ate * indicators))
cdf_balance_proj <- abs(base::colMeans(w_proj_ate * indicators))
cdf_balance_dps <- abs(base::colMeans(w_dps_ate * indicators))

cdf_balance_proj_1 <- abs(base::colMeans((w_proj1 - 1) * indicators))
cdf_balance_cbps_1 <- abs(base::colMeans((w_cbps1 - 1) * indicators))
cdf_balance_dps_1 <- abs(base::colMeans((w_dps1 - 1) * indicators))

cdf_balance_proj_0 <- abs(base::colMeans((w_proj0 - 1) * indicators))
cdf_balance_cbps_0 <- abs(base::colMeans((w_cbps0 - 1) * indicators))
cdf_balance_dps_0 <- abs(base::colMeans((w_dps0 - 1) * indicators))


ks_proj <-  sqrt(n) * c(
  max(cdf_balance_proj), 
  max(cdf_balance_proj_1), 
  max(cdf_balance_proj_0)
)

ks_cbps <-  sqrt(n) * c(
  max(cdf_balance_cbps), 
  max(cdf_balance_cbps_1), 
  max(cdf_balance_cbps_0)
)

ks_dps <-  sqrt(n) * c(
  max(cdf_balance_dps), 
  max(cdf_balance_dps_1), 
  max(cdf_balance_dps_0)
)

cvm_proj <-  c(
  sum(cdf_balance_proj^2), 
  sum(cdf_balance_proj_1^2), 
  sum(cdf_balance_proj_0^2)
)

cvm_cbps <-  c(
  sum(cdf_balance_cbps^2), 
  sum(cdf_balance_cbps_1^2), 
  sum(cdf_balance_cbps_0^2)
)

cvm_dps <-  c(
  sum(cdf_balance_dps^2), 
  sum(cdf_balance_dps_1^2), 
  sum(cdf_balance_dps_0^2)
)

out_p_bal <- matrix(0, ncol = 6, nrow = 3)
out_p_bal <- data.frame(out_p_bal)

out_p_bal[1,] <- 100 * c(ks_cbps/sqrt(n), sqrt(cvm_cbps/n))
out_p_bal[2,] <- 100 * c(ks_proj/sqrt(n), sqrt(cvm_proj/n))
out_p_bal[3,] <- 100 * c(ks_dps/sqrt(n), sqrt(cvm_dps/n))

out_p_bal <- t(out_p_bal)

rownames(out_p_bal) <- c("ks", "ks_1", "ks_0", "cvm", "cvm_1", "cvm_0")
colnames(out_p_bal) <- c("CBPS-just", "Proj", "DPS") 

out_p_bal


# reporting ---------------------------------------------------------------


res1_q <- data.frame(
  cbps = qte_cbps_exact$qte,
  proj = qte_ips_proj$qte,
  dps = qte_cbps_prop$qte) |>
  as.matrix()

res1_q_se <- data.frame(
  cbps = qte_cbps_exact$qte.se,
  proj = qte_ips_proj$qte.se,
  dps = qte_cbps_prop$qte.se) |>
  as.matrix()

tab_emp_net <- cbind(unlist(ate_cbps_exact[1:2]),
                     unlist(ate_ips_proj[1:2]),
                     unlist(ate_cbps_prop[1:2]))

tab_emp_net <- rbind(tab_emp_net,
                     res1_q[1,], res1_q_se[1,],
                     res1_q[2,], res1_q_se[2,],
                     res1_q[3,], res1_q_se[3,],
                     res1_q[4,], res1_q_se[4,],
                     res1_q[5,], res1_q_se[5,])



res1_q_tw <- data.frame(
  cbps = qte_cbps_exact_tw$qte,
  proj = qte_ips_proj_tw$qte,
  dps = qte_cbps_prop_tw$qte) |>
  as.matrix()

res1_q_tw_se <- data.frame(
  cbps = qte_cbps_exact_tw$qte.se,
  proj = qte_ips_proj_tw$qte.se,
  dps = qte_cbps_prop_tw$qte.se) |>
  as.matrix()

tab_tw <- cbind(unlist(ate_cbps_exact_tw[1:2]),
                     unlist(ate_ips_proj_tw[1:2]),
                     unlist(ate_cbps_prop_tw[1:2]))

tab_tw <- rbind(tab_tw,
                     res1_q_tw[1,], res1_q_tw_se[1,],
                     res1_q_tw[2,], res1_q_tw_se[2,],
                     res1_q_tw[3,], res1_q_tw_se[3,],
                     res1_q_tw[4,], res1_q_tw_se[4,],
                     res1_q_tw[5,], res1_q_tw_se[5,])



tab_res <- rbind(cbind(tab_emp_net, tab_tw), cbind(out_p_bal, out_p_bal))

print(xtable(tab_res, digits = 0),
      format.args = list(big.mark = ","))

       