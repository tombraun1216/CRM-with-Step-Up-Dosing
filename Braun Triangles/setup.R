setwd("/Users/tombraun/Work Documents/Collaborations/Mercier (Roche)/Braun Triangles/")
source("tools.R")
library(pracma) #For doing integration

#Targeted DLT rate
target <- 0.25

#There are eight unique doses examined in the motivating example
doseval <- c(6, 18, 30, 90, 270, 800, 2400, 7200)/1000

#Create a table of the doses used in the six schedules in the motivating example
sched <- matrix(doseval[c(1,2,2,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8)], 
                ncol=3, nrow=6, byrow=T)
nsched <- nrow(sched)

#Skeleton dose values
orb <- 1.5
orw <- 1.5
skel1 <- 0.03
for (i in 2:nsched)
  skel1[i] <- 1/((1-skel1[i-1])/(orb*skel1[i-1])+1)
skel2 <- 1/((1-skel1)/(orw*skel1)+1)
skel <- cbind(skel1, skel2, skel2)

#Days at which administration is possible
d <- 0:2*7

#Number of participants
Npart <- 30

#Arrival times via Poisson process
iat_mean <- 7
arriv_tm <- c(0, cumsum(round(rpois(Npart-1, iat_mean))))

#Matrix of administration times (days) for all participants
admin_tm <- cbind(arriv_tm, arriv_tm+7, arriv_tm+14)

#Should we randomize assignments?
rand_assn <- F

######################################################
#Duration of possible toxicity from one administration
TT <- 10

#Days at which administration is possible
d <- 0:2*7

#Mean for theta2, which has exponential distribution
lambda <- 0.09

#Investigator supplies info for distribution of theta1
#m +/- d, where m is expected day for peak toxicity
#d is margin of error
toxpeak <- 6
toxerr <- 3
theta1.mean <- toxpeak/TT
theta1.var <- (toxerr/(2*TT))^2
r <- theta1.mean*(1-theta1.mean)/theta1.var - 1
alpha <- theta1.mean*r
beta <- (1-theta1.mean)*r

#Visualize prior distribution for DLT probabilities
Nsamp <- 10000
ll <- rexp(Nsamp, 1/lambda)
tt <- rbeta(Nsamp, alpha, beta)*TT

prior_mu2 <- prior_sd2 <- NULL
for (k in 1:nrow(skel))
{
  p1 <- p2 <- p3 <- NULL
  for (i in 1:Nsamp)
  {
    p1 <- c(p1, cdf(d[1]+7, tt[i], ll[i], TT, d, skel[k,]))
    p2 <- c(p2, cdf(d[2]+7, tt[i], ll[i], TT, d, skel[k,]))
    p3 <- c(p3, cdf(d[3]+7, tt[i], ll[i], TT, d, skel[k,]))
  }
  prior_mu2 <- rbind(prior_mu2, c(mean(p1), mean(p2), mean(p3)))
  prior_sd2 <- rbind(prior_sd2, c(sd(p1), sd(p2), sd(p3)))
}

