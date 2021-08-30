#Run simulations
all_best_sched <- all_assn_sched <- all_obs_dlt <- all_ind_assn <- NULL
rmse <- NULL

Nsim <- 1000

for (s in 1:Nsim)
{ 
  tryCatch(
    {
      if (s%%(Nsim/20)==0) cat("On simulation ",s,"...\n" ,sep="")
      set.seed(121667+s)
      #Vector of DLT probability thresholds for all participants
      U <- runif(Npart, 0, 1)
      dlttm <- cbind(runif(Npart,0,d[2]), 
                     runif(Npart,d[2],d[3]), 
                     runif(Npart,d[3],d[3]+7))
      source("onesim.R")
      if (stop_trial) 
        cat("Simulation", s, "stopped - all schedules too toxic!\n")
      all_best_sched[s] <- best
      all_assn_sched <- rbind(all_assn_sched, assn_tbl)
      all_obs_dlt <- rbind(all_obs_dlt, obs_dlt)
      all_ind_assn <- rbind(all_ind_assn, assn)
      rmse <- c(rmse, sqrt(mean(ptox.mean-ptox[,3])^2))
    }, 
    error=function(e){cat("WARNING :",conditionMessage(e), "\n")})
}

best_sched <- (table(c(1:nsched, all_best_sched))-1)/Nsim
assn_sched <- apply(all_assn_sched, 2, mean)/Npart
ind_assn <- apply(all_ind_assn, 2, mean)

simdir <- "/Users/tombraun/Desktop/Sim Results/"
if (!dir.exists(simdir)) dir.create(simdir)
save.image(paste(simdir, "scen", scen_num, iat_mean, "t.RData", sep=""))
