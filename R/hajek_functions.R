#' Control variates method
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
#' 
#' @return A list containing the ATE est and associated variance
#' @export
#' 
controlVariatesMEHajek <- function(df, X,
                              Y,
                              A, Astar,
                              val_idx,
                              probs,
                              sl.lib,
                              estimand='ATE',
                              Ystar=NA
                              ) {
  
  
  # # Sort df so that the first n observations have val_idx=1
  # df <- df[order(df$val_idx, decreasing=T),] 
  # 
  # Keep track of the val data
  df_val <- df[df[,val_idx]==1,] # %>% dplyr::filter(val_idx==1)
  
  # Step 1: obtain point estimate in the validation data
  tau_val <- get_tau_hat_hajek(df_val,X,Y,A,probs,sl.lib)
  tau_hat_val <- tau_val$est ; v_hat <- tau_val$var_est
  
  # Step 2: get the control variate
  tau_ep_val <- get_tau_hat_hajek(df_val,X,Y,Astar,probs,sl.lib)
  df$const_prob <- 1
  tau_ep_main <- get_tau_hat_hajek(df,X,Y,Astar,'const_prob',sl.lib)
  
  # Step 3: estimate Gamma and V
  Sig_hat <- get_vcov_asym_hajek(tau_val$eif,tau_ep_val$eif,tau_ep_main$eif)
  gamma_hat <- Sig_hat$gamma_hat ; V_hat <- Sig_hat$V_hat
  
  # Step 4: subtract off (may need to flip the subtraction sign)
  tau_cv <- tau_hat_val - (gamma_hat/V_hat)*(tau_ep_val$est-tau_ep_main$est)
  
  # Get variance estimate and 95% CI
  var_hat <- v_hat - (gamma_hat^2/V_hat)/nrow(df_val)
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)
  
  # browser()
  
  # Return the ATE est and associated variance
  return(list(tau_cv=tau_cv,
              var_hat=var_hat,
              tau_hat_val=tau_hat_val,
              var_hat_val=v_hat,
              CI=c(ci_low,ci_high))
  )
  
}


#' Implement Hajék estimator -- rescale EIF of val. only estimator by function
#' of inverse sampling probs 
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
get_tau_hat_hajek <- function(df,X,Y,A,probs,sl.lib) {

  probs <- df[,probs]
  
  res <- AIPW::AIPW$new(Y=df[,Y],
                        A=df[,A],
                        W=subset(df,select=X),
                        Q.SL.library = sl.lib,
                        g.SL.library = sl.lib,
                        k_split = 1,
                        verbose=FALSE)$fit()$summary()
  
  # Extract EIF of val only ate est
  mod <- list()
  mod$est <- res$estimates$RD['Estimate']
  mod$se <- res$estimates$RD['SE']
  mod$eif1 <- res$obs_est$aipw_eif1 ; mod$eif0 <- res$obs_est$aipw_eif0
  mod$eif <- mod$eif1 - mod$eif0
  
  # Obtain final estimate via Hajék rescaling
  rescaled_eif <- probs^(-1)*mod$eif/mean(probs^(-1))
  # rescaled_eif <- mod$eif
  est <- mean(rescaled_eif) # mean(probs^(-1)*mod$eif)/mea (probs^(-1))
  
  # Get the SE estimate
  var_est <- var ( rescaled_eif ) / nrow(df)
  
  return(list(
    est=est,
    var_est=var_est,
    unscaled_eif=mod$eif,
    eif=probs^(-1)*mod$eif/mean(probs^(-1))
  )) 
}

#' Function for getting asymptotic variance of control variates estimator
#'
#' @param tau_val_mod A list containing the AIPW estimate of the estimand
#' @param cv_mods A list containing the AIPW estimate of the control variates
#' @param val_idx A vector of validation data indicators
#' 
#' @return A list containing estimates of V and Gamma
#' @export
get_vcov_asym_hajek <- function(varphi_val,phi_val,phi_main) {
  
  # Get sample sizes
  n_main <- length(phi_main) ; n_val <- length(phi_val)
  
  # Estimate Gamma
  gamma_hat <- (1 - n_val/n_main)*cov(cbind(demean(phi_val),
                                  demean(varphi_val)))[1,2]/n_val
  
  # Estimate V
  V_hat <- (1 - n_val/n_main)*var(phi_main)/n_val # mean(demean(phi_main)^2)/n_val
  
  return(list(gamma_hat=gamma_hat,
              V_hat=V_hat))
  
}

