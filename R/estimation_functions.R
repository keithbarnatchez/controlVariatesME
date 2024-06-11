#' Control variates method -- r package
#'
#' @param df A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariate names
#' @param Y A string containing the outcome name
#' @param A A string containing the treatment name
#' @param Astar A vector of error-prone treatment measurement names
#' @param Ystar A vector of error-prone outcome measurement names
#' @param val_idx A vector of validation data indicators
#' @param estimand The estimand of interest (e.g. ATE, ATT, ATC)
#' @param probs A vector of probabilities for sampling into val data
#' @param sl.lib A vector of super learner libraries to use
#' @returns A list containing the ATE estimate and associated variance estimate
#' @examples
#' controlVariatesME(df, X='X', Y='Y', A='A', Astar='Astar', val_idx='S', sl.lib='SL.ranger')
#' 
#' 
#' @return A list containing the ATE est and associated variance
#' @export
#' 
controlVariatesME <- function(df, X,
                       Y,
                       A, Astar,
                       val_idx,
                       sl.lib,
                       estimand='ATE',
                       Ystar=NA
                       ) {

  
  # Keep track of val data
  df_val <- df[df[,val_idx]==1,] 
  
  # Step 1: obtain estimate in the validation data
  tau_val_mod <- get_tau_hat(df_val,X,Y,A,estimand,sl.lib)
  tau_hat_val <- tau_val_mod$estimates$RD['Estimate']
  v_hat <- var(tau_val_mod$obs_est$aipw_eif1 - 
                 tau_val_mod$obs_est$aipw_eif0)/nrow(df_val)
  
  # Step 2: construct control variates
  if (is.na(Ystar)) { # if no error-prone outcome
    cv_mods <- get_psi(df,X,Y,Astar,estimand,sl.lib)
  } else {
    cv_mods <- get_psi(df,X,Ystar,Astar,estimand,sl.lib)
  }
  tau_ep_main <- cv_mods$psi1$estimates$RD['Estimate'] 
  tau_ep_val <- cv_mods$psi2$estimates$RD['Estimate']
  
  # Step 3: estimate Gamma and V
  Sig_hat <- get_vcov_asym(tau_val_mod,cv_mods,df$val_idx)
  gamma_hat <- Sig_hat$gamma_hat ; V_hat <- Sig_hat$V_hat
  
  # Step 4: form the final estimate
  tau_cv <- tau_hat_val - (gamma_hat/V_hat)*(tau_ep_val-tau_ep_main)
  
  # Get variance estimate and 95% CI
  var_hat <- v_hat - gamma_hat^2/V_hat
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)
  
  # Return the ATE est and associated variance
  return(list(tau_cv=tau_cv,
              var_hat=var_hat,
              CI=c(ci_low,ci_high))
  )
  
}

demean <- function(x) {
  return(x - mean(x))
}

#' Get val. estimate
#'
#' @param df A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param Y A vector of outcomes
#' @param A A vector of treatments
#' @param estimand The estimand of interest (e.g. ATE, ATT, ATC)
#' @param sl.lib A vector of super learner libraries to use
#' 
#' @return A list containing the AIPW estimate of the estimand
#' @export
#' 
get_tau_hat <- function(df,X,Y,A,estimand,sl.lib) {

  res <- AIPW::AIPW$new(Y=df[,Y],
                  A=df[,A],
                  W=subset(df,select=X),
                  Q.SL.library = sl.lib,
                  g.SL.library = sl.lib,
                  k_split = 1,
                  verbose=FALSE)$fit()$summary()
  
  return(res) 
}

#' Function for getting asymptotic variance of control variates estimator
#'
#' @param tau_val_mod A list containing the AIPW estimate of the estimand
#' @param cv_mods A list containing the AIPW estimate of the control variates
#' @param val_idx A vector of validation data indicators
#' 
#' @return A list containing estimates of V and Gamma
#' @export
get_vcov_asym <- function(tau_val_mod,cv_mods,
                          val_idx) {
 
  # Get influence function estimates
  phi_main <- cv_mods[[1]]$obs_est$aipw_eif1 - cv_mods[[1]]$obs_est$aipw_eif0 # error prone main
  phi_val <- cv_mods[[2]]$obs_est$aipw_eif1 - cv_mods[[2]]$obs_est$aipw_eif0 # error prone val
  varphi_val <- tau_val_mod$obs_est$aipw_eif1 - tau_val_mod$obs_est$aipw_eif0 # unbiased val
  
  # Get sample sizes
  n_main <- length(phi_main) ; n_val <- length(phi_val)
  
  # Estimate Gamma
  gamma_hat <- (1 - n_val/n_main)*(1/n_val)*cov(cbind(demean(phi_val),demean(varphi_val)))[1,2]
  
  # Estimate V
  V_hat <- (1 - n_val/n_main)*(1/n_val)*mean(demean(phi_main)^2)
  
  return(list(gamma_hat=gamma_hat,
              V_hat=V_hat))
  
}

#' get_psi
#' 
#' @param df A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param Y A vector of outcomes
#' @param A A vector of treatments
#' @param estimand The estimand of interest (e.g. ATE, ATT, ATC)
#' @param sl.lib A vector of super learner libraries to use
#' 
#' @return A list containing the AIPW estimate of the control variates
#' @export
get_psi <- function(df,X,Y,Astar,estimand,sl.lib) {

  dfs <- list()
  dfs[[1]] <- df # full data
  dfs[[2]] <- df[(df[,val_idx]==1),] # validation data
  
  psi1 <- AIPW$new(Y=dfs[[1]][,Y],
                   A=dfs[[1]][,Astar],
                   W=subset(dfs[[1]],select=X),
                   Q.SL.library = sl.lib,
                   g.SL.library = sl.lib,
                   k_split = 1,
                   verbose=FALSE)$fit()$summary() 
  
  psi2 <- AIPW$new(Y=dfs[[2]][,Y],
                   A=dfs[[2]][,Astar],
                   W=subset(dfs[[2]],select=X),
                   Q.SL.library = sl.lib,
                   g.SL.library = sl.lib,
                   k_split = 1,
                   verbose=FALSE)$fit()$summary() 
  
  res_list = list(psi1=psi1,psi2=psi2)
  
  return(res_list)
  
}

