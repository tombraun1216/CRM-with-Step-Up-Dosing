denom <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout,
                  mu, sigma, a_t, b_t, a_g, b_g)
{
  prior <- dnorm(beta, mu, sigma)*dgamma(theta, a_t, b_t)*dgamma(gamma, a_g, b_g)
  exp(log_llh(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout))*prior        
}

numer_p1 <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout,
                     mu, sigma, a_t, b_t, a_g, b_g, doseval)
{
  prior <- dnorm(beta, mu, sigma)*dgamma(theta, a_t, b_t)*dgamma(gamma, a_g, b_g)
  prob <- doseval^(exp(beta))
  prob*exp(log_llh(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout))*prior        
}

numer_p2 <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout,
                     mu, sigma, a_t, b_t, a_g, b_g, doseval)
{
  prior <- dnorm(beta, mu, sigma)*dgamma(theta, a_t, b_t)*dgamma(gamma, a_g, b_g)
  prob <- doseval^(exp(beta-theta))
  prob*exp(log_llh(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout))*prior        
}

numer_p3 <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout,
                    mu, sigma, a_t, b_t, a_g, b_g, doseval)
{
  prior <- dnorm(beta, mu, sigma)*dgamma(theta, a_t, b_t)*dgamma(gamma, a_g, b_g)
  prob <- doseval^(exp(beta-theta-gamma))
  prob*exp(log_llh(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout))*prior        
}

numer_stp <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout,
                      mu, sigma, a_t, b_t, a_g, b_g)
{
  prior <- dnorm(beta+theta+gamma, mu, sigma)*dgamma(theta, a_t, b_t)*dgamma(gamma, a_g, b_g)
  exp(log_llh(beta+theta+gamma, theta, gamma, Yout, Dout, Wout, Zout, Nout))*prior        
}

log_llh <- function(beta, theta, gamma, Yout, Dout, Wout, Zout, Nout)
{
  e1 <- exp(beta)
  e2 <- exp(beta - theta)
  e3 <- exp(beta - theta - gamma)
  llh <- 0
  for (j in 1:Nout)
  {
    p11 <- p10 <- p21 <- p20 <- p31 <- p30 <- 1
    if(Zout[j]>=1)
    { 
      p11 <- Wout[j,1]*Dout[j,1]^e1
      p10 <- 1-p11
    }
    
    if (Zout[j]>=2)
    {
      p21 <- Wout[j,2]*(Dout[j,2]^e2 - Dout[j,1]^e1)
      p20 <- (1-Dout[j,1]^e1)-Wout[j,2]*(Dout[j,2]^e2-Dout[j,1]^e1)
    }
    
    if (Zout[j]==3)
    {
      p31 <- Wout[j,3]*(Dout[j,3]^e3 - Dout[j,2]^e2)
      p30 <- (1-Dout[j,2]^e2)-Wout[j,3]*(Dout[j,3]^e3-Dout[j,2]^e2)
    }
    
    if (Yout[j]==1) 
      llh_inc <- (Zout[j]==1)*log(p11) + (Zout[j]==2)*log(p21) + (Zout[j]==3)*log(p31)
    if (Yout[j]==0) 
      llh_inc <- (Zout[j]==1)*log(p10) + (Zout[j]==2)*log(p20) + (Zout[j]==3)*log(p30)
    llh <- llh + llh_inc
  }
  llh
}