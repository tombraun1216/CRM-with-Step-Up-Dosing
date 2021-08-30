#Respective matrices of DLT indicators, dose indicators, weights 
#and IDs of all participants
ids <- dlt <- doses <- NULL

#Initialize simulation for first participant (who is assigned first schedule)
i <- 1
assn <- 1
max_comp_sched <- 0
stop_trial <- F

#Recursively evaluate MTD and enroll next subject
while (i < Npart)
{
   #Enroll first subject; record their toxicity time
   ids <- c(ids, i)
   xx <- as.numeric(U[i]<=ptox[assn[i],])
   obstm <- pmin(21, arriv_tm[i+1] - arriv_tm[1:i])
   dlt <- rbind(dlt, xx)
   doses <- rbind(doses, skel[assn[i],])
   
   Tout <- rbind(sapply(1:i, get_tm))
   Yout <- rbind(as.numeric(Tout!=obstm & Tout<21))
   Dout <- rbind(doses)

   #Determine which schedules have one participant with complete follow-up
   comp_sched <- assn[Yout==1 | Tout==21]
   if (length(comp_sched)>0) max_comp_sched <- max(comp_sched)
   
   #Do posterior computations via integral approximation
   den <- integrate(post0, 0, Inf, Tout=Tout, Dout=Dout, Yout=Yout, 
                    admin=d, Ke=ke_fix, Keff=keff_fix, beta_mu=beta_mu, 
                    beta_sd=beta_sd)$value
   num <- NULL
   for (uu in 1:nsched)
   {
      dv <- skel[uu,]
      num[uu] <- integrate(post1, 0, Inf, Tout=Tout, Dout=Dout, Yout=Yout, 
                           admin=d, doseval=dv, Ke=ke_fix, Keff=keff_fix, 
                           beta_mu=beta_mu, beta_sd=beta_sd)$value
   }
   ptox.mean <- num/den   
   
   #Determine the best schedule
   delta <- abs(ptox.mean-target)
   best <- (1:nsched)[delta==min(delta)]
   
   #Randomize assignments, if indicated
   if (rand_assn)
   {
      nhood <- best+(-1:1)
      nhood <- nhood[nhood %in% (1:nsched)]
      pwr <- 1/2
      p_rand <- 1/(delta[nhood]+0.001)^pwr/sum(1/(delta[nhood]+0.001)^pwr)
      best <- sample(nhood, 1, prob=p_rand)
   }
   
   next_assn <- min(best, max_comp_sched+1)
   
   #Stop trial if lowest dose looks too toxic
   p1_lb <- binom.test(sum(Yout[assn==1]), sum(assn==1), 
                       p=target, alt="g")$conf.int[1]
   if (p1_lb>target)
   {
      Nout <- i 
      i <- Npart
      stop_trial <- T
   }
   
   #Assign schedule to next participant and move forward in time
   if (p1_lb<=target)
   {
      assn <- c(assn, next_assn)
      i <- i+1
   }
}

if (stop_trial)
{
   assn <- c(assn, rep(NA, Npart-Nout))
   assn_tbl <- table(c(1:nsched, assn))-1
   obs_dlt <- c(Yout, rep(NA, Npart-Nout))
   best <- NA
}

#Now make final computation using remaining follow-up of all participants
if(!stop_trial)
{
 ids <- c(ids, i)
 xx <- as.numeric(U[i]<=ptox[assn[i],])
 obstm <- (arriv_tm[i]+21) - arriv_tm[1:i] 

 dlt <- rbind(dlt, xx)
 doses <- rbind(doses, skel[assn[i],])

 Tout <- rbind(pmin(21, sapply(1:i, get_tm)))
 Yout <- rbind(as.numeric(Tout!=obstm & Tout<21))
 Dout <- rbind(doses)
 
 #Compute posterior means
 den <- integrate(post0, 0, Inf, Tout=Tout, Dout=Dout, Yout=Yout, admin=d, 
                  Ke=ke_fix, Keff=keff_fix, beta_mu=beta_mu, beta_sd=beta_sd)$value
 num <- NULL
 for (uu in 1:nsched)
 {
    dv <- skel[uu,]
    num[uu] <- integrate(post1, 0, Inf, Tout=Tout, Dout=Dout, Yout=Yout,
                         admin=d, doseval=dv, Ke=ke_fix, Keff=keff_fix,
                         beta_mu=beta_mu, beta_sd=beta_sd)$value
 }
 ptox.mean <- num/den    

 #Determine the best schedule
 delta <- abs(ptox.mean-target)
 best <- (1:nsched)[delta==min(delta)]
 assn_tbl <- table(c(1:nsched, assn))-1
 obs_dlt <- Yout
}