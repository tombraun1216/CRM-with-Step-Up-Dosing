one.triangle <- function(i, t, theta1, theta2, TT, admin, doses)
{
  delta <- t-admin[i]
  i1 <- delta>TT
  i2 <- (delta>theta1 & delta<=TT)
  i3 <- (delta>0 & delta<=theta1)
  i4 <- (delta<=0)
  
  t2 <- doses[i]*theta2
  a1 <- TT*t2/2
  a2 <- theta1*t2/2+(TT-theta1)*t2/2 - (TT-delta)^2*t2/(TT-theta1)/2
  a3 <- t2/(2*theta1)*delta^2
  a4 <- 0
  area <- i1*a1+i2*a2+i3*a3+i4*a4

  h1 <- 0
  h2 <- t2/(TT-theta1)*(TT-delta)
  h3 <- t2/theta1*delta
  h4 <- 0
  hgt <- i1*h1+i2*h2+i3*h3+i4*h4
  c(area, hgt)
} 

all.triangle <- function(t, theta1, theta2, TT, admin, doses)
{
  sapply(1:3, one.triangle, t=t, theta1=theta1, theta2=theta2, TT=TT, admin=admin,
         doses=doses)
}

#CDF at a given timepoint
cdf <- function(t, theta1, theta2, TT, admin, doses)
{
  all <- all.triangle(t, theta1, theta2, TT, admin, doses)[1,]
  1-exp(-sum(all))
}

pdf <- function(t, the1, the2, TT, admin, doses)
{   
  all <- all.triangle(t, the1, the2, TT, admin, doses)[2,]
  sum(all)*(1-cdf(t, the1, the2, TT, admin, doses))
}

#Likelihood for all participants
L1 <- function(Tout, Yout, Dout, the1, the2, TT, admin)
{
  sapply(1:length(Tout), L1a, Tout=Tout, Yout=Yout, Dout=Dout,
         the1=the1, the2=the2, TT=TT, admin=admin)
}

#Likelihood (either hazard or survival function)
L1a <- function(i, Tout, Yout, Dout, the1, the2, TT, admin)
{
  a1 <- cdf(Tout[i], the1, the2, TT, admin, Dout[i,])
  a2 <- pdf(Tout[i], the1, the2, TT, admin, Dout[i,])
  a2*Yout[i] + (1-a1)*(1-Yout[i])
}

#Denominator for posterior probability of DLT (normalizing)
post0 <- function(the1, the2, Tout, Yout, Dout, TT, admin,
                  alpha, beta, lambda)
{
 a1 <- prod(L1(Tout, Yout, Dout, the1, the2, TT, admin))
 a2 <- exp(-the2/lambda)*(the1/TT)^(alpha-1)*(1-the1/TT)^(beta-1)
 a1*a2
}

#Numerator of posterior mean for probability of DLT by Day 21 for given schedule
post1 <- function(the1, the2, Tout, Yout, Dout, TT, admin, sched,
                  alpha, beta, lambda)
{
  a1 <- prod(L1(Tout, Yout, Dout, the1, the2, TT, admin))
  a2 <- exp(-the2/lambda)*(the1/TT)^(alpha-1)*(1-the1/TT)^(beta-1)
  a1*a2*cdf(21, the1, the2, TT, admin, skel[sched,])
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
