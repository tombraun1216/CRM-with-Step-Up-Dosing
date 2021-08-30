#Unnormalized exposure function
get_exp <- function(dose_vect, t_curr, k_eff, k_e, admin)
{
  #This computes unnormalized AUC values as described in GWF, page 4
  #dose_vect = vector of three doses in schedule assigned to participant  
  #t_curr = current timepoint for which AUC is needed
  #k_eff = PK parameter used in Equation 7
  #k_e = PK parameter used in Equation 7
  #Time since each administration was given relative to current time
  t_vect <- pmax(0, t_curr-admin)
  
  #AUC for administrations 1, 2, and 3, respectively of doses in dose_vect
  exp1 <- dose_vect[1]*(exp(-k_e*t_vect[1]) - exp(-k_eff*t_vect[1]))
  exp2 <- dose_vect[2]*(exp(-k_e*t_vect[2]) - exp(-k_eff*t_vect[2]))
  exp3 <- dose_vect[3]*(exp(-k_e*t_vect[3]) - exp(-k_eff*t_vect[3]))
  #Total exposure
  k_eff/(k_eff - k_e)*(exp1+exp2+exp3)
}

#Unnormalized AUC function
get_auc <- function(dose_vect, t_curr, k_eff, k_e, admin)
{
  #This computes unnormalized AUC values as described in GWF, page 4
  #dose_vect = vector of three doses in schedule assigned to participant  
  #t_curr = current timepoint for which AUC is needed
  #k_eff = PK parameter used in Equation 7
  #k_e = PK parameter used in Equation 7
  #Time since each administration was given relative to current time
  t_vect <- pmax(0, t_curr-admin)
  
  #AUC for administrations 1, 2, and 3, respectively of doses in dose_vect
  auc1 <- dose_vect[1]*((1-exp(-k_e*t_vect[1]))/k_e - 
                        (1-exp(-k_eff*t_vect[1]))/k_eff)
  auc2 <- dose_vect[2]*((1-exp(-k_e*t_vect[2]))/k_e - 
                        (1-exp(-k_eff*t_vect[2]))/k_eff)
  auc3 <- dose_vect[3]*((1-exp(-k_e*t_vect[3]))/k_e - 
                        (1-exp(-k_eff*t_vect[3]))/k_eff)
#Total AUC
k_eff/(k_eff - k_e)*(auc1+auc2+auc3)
}

#Likelihood for all participants
L1 <- function(Tout, Yout, Dout, Ke, Keff, beta, admin)
{
  sapply(1:length(Tout), L1a, Tout=Tout, Yout=Yout, Dout=Dout,
         Ke=Ke, Keff=Keff, beta=beta, admin=admin)
}

#Likelihood for each participant
L1a <- function(i, Tout, Yout, Dout, Ke, Keff, beta, admin)
{
  exp0 <- get_exp(Dout[i,], Tout[i], Ke, Keff, admin)
  auc0 <- get_auc(Dout[i,], Tout[i], Ke, Keff, admin)
  H0 <- beta*auc0
  h0 <- beta*exp0
  S0 <- exp(-H0)
  f0 <- S0*h0
  f0*Yout[i] + S0*(1-Yout[i])
}

#Denominator for posterior means (normalizing constant)
post0 <- function(beta, Tout, Yout, Dout, admin, Ke, Keff, beta_mu, beta_sd)
{
  llh <- NULL
  for (j in 1:length(beta))
    llh <- c(llh, prod(L1(Tout, Yout, Dout, Ke, Keff, beta[j], admin)))
  prior_beta <- dlnorm(beta, beta_mu, beta_sd)
  llh*prior_beta
}

#Numerator of posterior mean for DLT probability for given schedule
post1 <- function(beta, Tout, Yout, Dout, admin, doseval, Ke, Keff, beta_mu, beta_sd)
{
  llh <- NULL
  for (j in 1:length(beta))
    llh <- c(llh, prod(L1(Tout, Yout, Dout, Ke, Keff, beta[j], admin)))
  prior_beta <- dlnorm(beta, beta_mu, beta_sd)
  prior_beta <- dlnorm(beta, beta_mu, beta_sd)
  pp <- 1-L1(21, 0, rbind(doseval), Ke, Keff, beta, admin)
  pp*llh*prior_beta
}

#Used recursively to update follow-up time of each participant
get_tm <- function(i)
{
  out_tm <- obstm[i]
  had_dlt <- (dlttm[i,]<obstm[i])*(dlt[i,]==1)
  if (max(had_dlt)==1) 
    out_tm <- min(dlttm[i,][had_dlt==1])
  out_tm
}   
