#Respective matrices of DLT indicators, dose indicators, weights 
#and IDs of all participants
ids <- dlt <- doses <- NULL

#Initialize simulation for first participant (who is assigned first schedule)
i <- 1
assn <- 1
max_comp_sched <- 0
stop_trial <- F

while (i < Npart)
{
 #Collect DLT and dose information on newest participant
 ids <- c(ids, i)
 xx <- as.numeric(U[i]<=ptox[assn[i],])
 if (xx[1]==1) xx[2] <- xx[3] <- NA
 if (xx[1]==0 & xx[2]==1) xx[3] <- NA
 dlt <- rbind(dlt, xx)
 doses <- rbind(doses, skel[assn[i],])

 #Prepare accrued data for analysis
 ok <- ifelse(admin_tm<arriv_tm[i+1], 1, NA)
 Yout <- dlt[1:i,]*ok[1:i,]
 Dout <- doses[1:i,]*ok[1:i,]
 Wout <- pmin((arriv_tm[i+1]-admin_tm[1:i,])/7, 1)*ok[1:i,]
 Wout <- ifelse(Yout==1, 1, Wout)
 
 #Vectorize data when only one participant
 if (i==1)
 {
   Yout <- rbind(Yout)
   Dout <- rbind(Dout)
   Wout <- rbind(Wout)
 }

 Zout <- apply(!is.na(Yout), 1, sum)
 Yout <- apply(Yout, 1, max, na.rm=T)
 Nout <- i 
 
 #Determine which schedules have one participant with complete follow-up
 comp_sched <- assn[apply(Wout, 1, min, na.rm=T)==1]
 if (length(comp_sched)>0) max_comp_sched <- max(comp_sched)

 #Do posterior sampling/computation 
 source('run_integ.R')
 
 #Determine which schedule is current best schedule, subject to restrictions
 delta <- abs(p3_post-target)
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

if(stop_trial)
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
 if (xx[1]==1) xx[2] <- xx[3] <- NA
 if (xx[1]==0 & xx[2]==1) xx[3] <- NA
 dlt <- rbind(dlt, xx)
 doses <- rbind(doses, skel[assn[i],])
 Yout <- apply(dlt, 1, max, na.rm=T)
 Zout <- apply(!is.na(dlt), 1, sum)
 Dout <- doses
 Wout <- 1-is.na(Dout)
 Nout <- Npart
 source('run_integ.R')
 delta <- abs(p3_post-target)
 best <- (1:nsched)[delta==min(delta)]
 assn_tbl <- table(c(1:nsched, assn))-1
 obs_dlt <- Yout
}

