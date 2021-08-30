setwd("/Users/tombraun/Work Documents/Collaborations/Mercier (Roche)/Gunhan PK/")
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
#Define mean and standard deviation of beta
beta_mu <- -2
beta_sd <- 1

#These parameters are fixed
ke_fix <- log(2)/4*2
keff_fix <- exp(-2)

#Visualize prior distribution for DLT probabilities
Nsamp <- 10000
beta_samp <- rlnorm(Nsamp, beta_mu, beta_sd)

p1 <- p2 <- p3 <- NULL
for (j in 1:nsched)
{
  p1 <- cbind(p1, 1-exp(-beta_samp*get_auc(skel[j,], 7, keff_fix, ke_fix, d)))
  p2 <- cbind(p2, 1-exp(-beta_samp*get_auc(skel[j,], 14, keff_fix, ke_fix, d)))
  p3 <- cbind(p3, 1-exp(-beta_samp*get_auc(skel[j,], 21, keff_fix, ke_fix, d)))
}

prior_mu3 <- cbind(apply(p1, 2, mean), apply(p2, 2, mean), apply(p3, 2, mean))
prior_sd3 <- cbind(apply(p1, 2, sd), apply(p2, 2, sd), apply(p3, 2, sd))

