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
iat_mean <- 21
arriv_tm <- c(0, cumsum(round(rpois(Npart-1, iat_mean))))

#Matrix of administration times (days) for all participants
admin_tm <- cbind(arriv_tm, arriv_tm+7, arriv_tm+14)

#Should we randomize assignments?
rand_assn <- F

######################################################
#Parameter values for prior distributions
#For beta
#mu <- prod(log(-log(skel[,1])))^(1/nsched)
mu <- mean(log(-log(skel[,1])))
sigma <- 1#sqrt(mu)

#For theta and gamma
ff <- 1.6*2
mu_t <- -log((log(target)-log(ff))/(log(target)-2*log(ff)))
#ub_t <- -log((log(target+0.10)-log(ff))/(log(target+0.10)-2*log(ff)))
#sd_t <- (ub_t-mu_t)/2
a_t <- 1#mu_t^2/sd_t^2
b_t <- a_t/mu_t

mu_g <- -log(log(target)/(log(target)-log(ff)))
#ub_g <- -log(log(target+0.10)/(log(target+0.10)-log(ff)))
#sd_g <- (ub_g-mu_g)/2
a_g <- 1#mu_g^2/sd_g^2
b_g <- a_g/mu_g

#Visualize prior distribution for DLT probabilities
bb <- rnorm(10000, mu, sigma)
tt <- rgamma(10000, a_t, b_t)
gg <- rgamma(10000, a_g, b_g)

prior_mu1 <- prior_sd1 <- NULL
for (i in 1:6)
{
  p1 <- skel[i,1]^exp(bb)
  p2 <- skel[i,2]^exp(bb-tt)
  p3 <- skel[i,3]^exp(bb-tt-gg)
  prior_mu1 <- rbind(prior_mu1, c(mean(p1), mean(p2), mean(p3)))
  prior_sd1 <- rbind(prior_sd1, c(sd(p1), sd(p2), sd(p3)))
}