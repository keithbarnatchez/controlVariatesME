 
#' Function for checking if a vector is binary 
#' 
#' @param x A vector of values
#' @return A boolean indicating if the vector is binary
#' @export
is_binary <- function(x) {
  
  
  unique_values <- unique(x)
  if (length(unique_values) <= 2 && all(unique_values %in% c(0, 1))) {
    return(TRUE)
  } else if (length(unique_values) == 1 && unique_values %in% c(0, 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' Function for estimating the outcome model
#'
#' @param data A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param Y A vector of outcomes
#' @param A A vector of treatments
#' @param val_idx A vector of validation data indicators
#' @param sl.lib A vector of super learner libraries to use
#' 
#' @return A list containing the estimated outcome model
#' @export
mu_hat_est <-function(data,X,Y,A,val_idx,
                      sl.lib) {

  # First, make sure to only fit on validation data
  data_fit <- data[data[,val_idx]==1,] # %>% filter(val_idx==1)
  
  Y = data_fit[,Y]
  outcovs = data_fit[,c(X,A)]
  outcovs_main = data[,c(X,A)]
  
  # Check is y is binary or not
  fam <- 'gaussian'
  if (is_binary(Y)) {fam <- 'binomial'}
  
  # Fit the outcome model E[Y | A=a,X,S=1]
  mu_mod <- SuperLearner(Y=Y,
                         X=outcovs,
                         SL.library=sl.lib,
                         family=fam)

  # Get predicted values under treatment
  XA1 <- outcovs_main ; XA1[,A] <- 1
  mu_hat1 <- predict(mu_mod,
                     XA1,
                     onlySL=TRUE)
  
  # Get predicted values under no treatment
  XA0 <- outcovs_main ; XA0[,A] <- 0
  mu_hat0 <- predict(mu_mod,
                     XA0,
                     onlySL=TRUE)
  
  return(list(mu_hat1=mu_hat1, mu_hat0 = mu_hat0))
}

#' Function for estimating the propensity score model
#'
#' @param data A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param A A vector of treatments
#' @param sl.lib A vector of super learner libraries to use
#' 
#' @return A vector containing the estimated propensity score model
#' @export
pi_hat_est <- function(data,X,A,val_idx,
                       sl.lib) {
  
  # Only fit with val data
  data_fit <- data[data[,val_idx]==1,]
  A <- data_fit[,A] ; Xmain <- data[,X] ; X <- data_fit[,X] 
  
  
  # Fit the ps model on the validation data
  pi_mod <- SuperLearner::SuperLearner(Y=A,
                                       X=X,
                                       SL.library =sl.lib,
                                       family='binomial')
  
  # Get predicted values over the whole dataset
  pi_hat <- predict(pi_mod,Xmain,onlySL=TRUE,type='response')
  
  return( pi_hat )
  
}

#' Function for estimating the probability of selection into the validation set
#'
#' @param data A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param val_idx A vector of validation indices
#' @param sl.lib A vector of super learner libraries to use
#' 
#' @return A vector containing the estimated probability of selection into the validation set
#' @export
kappa_hat_est <- function(data,
                          X,val_idx,
                          sl.lib) {

  
  # Fit the ps model on the validation data
  S <- data[,val_idx] ; Xmain <- data[,X]
  kappa_mod <- SuperLearner::SuperLearner(Y=S,
                                          X=Xmain,
                                          SL.library =sl.lib,
                                          family='binomial')
  
  # Get predicted values over the whole dataset
  kappa_hat <- predict(kappa_mod,onlySL=TRUE)
  
  return(kappa_hat)
  
}

#' Function for obtaining the generalizability efficient influence function
#' terms
#'
#' @param mu_hat1 A vector of predicted outcomes under treatment
#' @param mu_hat0 A vector of predicted outcomes under no treatment
#' @param pi_hat A vector of predicted propensity scores
#' @param kappa_hat A vector of predicted selection probabilities
#' @param data A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariates
#' @param Y A vector of outcomes
#' @param A A vector of treatments
#' @param val_idx A vector of validation indices
#' 
#' @return A vector containing the estimated EIF
#' @export
generalizability_eif <- function(mu_hat1, mu_hat0,
                                 pi_hat,
                                 kappa_hat,
                                 data,
                                 X,Y,A,val_idx) {

  
  S <- data[,val_idx] ; A <- data[,A] ; Y <- data[,Y]
  nmain <- length(S) ; nval <- sum(S)
  
  # EIF for E(E(Y|A=1,S=1,X))
  psihat1 <-  (S*A*(Y-mu_hat1)/(kappa_hat*pi_hat) +
                 mu_hat1)
  
  # EIF for E(E(Y|A=0,S=1,X))
  psihat0 <- (S*(1-A)*(Y-mu_hat0)/(kappa_hat*(1-pi_hat)) +
                mu_hat0)
  
  return(psihat1 - psihat0) # ATE EIF is diff of these 2 terms
  
}


#' Function for implementing generalizability estimator
#'
#' @param data A dataframe containing the outcome, treatment and covariates
#' @param X A vector of covariate names
#' @param Y A vector of outcome names
#' @param A A vector of treatments
#' @param val_idx A vector of validation indices
#' @param sl.lib A vector of super learner libraries to use
#'
#' @return A vector containing the estimated EIF
#' @export
generalizability_est <- function(data, 
                                 X,Y,A,val_idx,
                                 sl.lib,
                                 rho=NA) {
  
  # Get nuisance model estimates
  # Outcome model
  mu_hat_ests <- mu_hat_est(data, X, Y, A, val_idx, sl.lib)
  mu_hat1 <- mu_hat_ests$mu_hat1$pred
  mu_hat0 <- mu_hat_ests$mu_hat0$pred
  
  # Propensity score model
  pi_hat <- pi_hat_est(data,X,A,val_idx,sl.lib)$pred
  
  # Participation model
  # when rho is missing, implies participation probs are not a constant function
  # but instead depend on X. if rho supplied, can just use it as "estimate"
  # May make more sense to just always model kappa, even when rho is fixed, along
  # similar lines to IPW out-performing the oracle under proper model form
  kappa_hat <- rho
  if (is.na(rho)) { 
    kappa_hat <- kappa_hat_est(data,X,val_idx,sl.lib)$pred
  }
  
  # Form EIF
  genz_eif <- generalizability_eif(mu_hat1, mu_hat0,
                                   pi_hat,
                                   kappa_hat,
                                   data,
                                   X,Y,A,val_idx)
  
  # Use EIF to get ATE est and variance est
  ate_est <- mean(genz_eif)
  var_est <- var(genz_eif)/length(genz_eif) # sqrt(length(genz_eif))
  
  return(list(ATE=ate_est,
              EIF=genz_eif,
              v_hat=var_est))
  
}


#' Main function for implementing the "generalizability" version of the control
#' variates estimator
#' 
#' @param data A dataframe 
#' @param X A vector of covariates
#' @param Y A vector of outcomes
#' @param A A vector of treatments
#' @param val_idx A vector of validation indices
#' @param sl.lib A vector of super learner libraries to use
#' @param estimand The estimand of interest (either 'ATE' or 'ATT')
#' @param rho A vector of participation probabilities
#' 
#' @return A list containing the ATE estimate, corresponding variance estimate,
#' and estimates of Gamma and V from the asymptotic covariance matrix
#' @export
controlVariatesMEGen <- function(data,
                                 X,Y,A,Astar,val_idx,
                                 sl.lib,
                                 estimand='ATE',
                                 rho=NA) {

  n <- nrow(data)
  
  # First, get ATE estimate pre-variance reduction
  ate_eif_genz <- generalizability_est(data,
                                       X,Y,A,val_idx,
                                       sl.lib,
                                       rho)
  tau_hat_val <- ate_eif_genz$ATE
  v_hat_val <- ate_eif_genz$v_hat
  
  # Next, get error-prone ATE estimates (val data -> main)
  ate_eif_ep_val_genz <- generalizability_est(data,
                                              X,Y,Astar,val_idx,
                                              sl.lib,
                                              rho)
  tau_hat_ep_val <- ate_eif_ep_val_genz$ATE
  
  # Next, get error-prone ATE estimates over the main data
  ate_eif_ep_main <- get_tau_hat(data,X,Y,Astar,estimand,sl.lib)
  tau_hat_ep_main <- ate_eif_ep_main$estimates$RD['Estimate']
  
  # Now, construct CV estimator. First, get Gamma
  eif_val <- ate_eif_genz$EIF
  eif_val_ep <- ate_eif_ep_val_genz$EIF
  eif_main_ep <- ate_eif_ep_main$obs_est$aipw_eif1 -
    ate_eif_ep_main$obs_est$aipw_eif0
  gamma_hat <- 1/n * (cov(eif_val,eif_main_ep) -
                        cov(eif_val,eif_val_ep))
  
  # Next, get V
  V_hat <- 1/n * var(eif_val_ep + eif_main_ep)
  # print(gamma_hat/V_hat)
  
  # Finally, form the CV estimator
  tau_cv <- tau_hat_val - gamma_hat/V_hat * (tau_hat_ep_main - tau_hat_ep_val)
  
  # Get variance estimate for CV est variance
  var_hat <- v_hat_val - gamma_hat^2/V_hat
  
  return(list(tau_cv=tau_cv,
              var_hat=var_hat,
              gamma_hat=gamma_hat,
              V_hat=V_hat))
}
