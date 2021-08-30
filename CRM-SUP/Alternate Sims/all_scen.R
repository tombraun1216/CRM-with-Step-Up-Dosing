setwd("/Users/tombraun/Work Documents/Collaborations/Mercier (Roche)/CRM-SUP/Alternate Sims/")

#Scenario 3 - best schedule is 3
#rm(list=ls())
#source("setup_sampsize1.R")
#scen_num <- 3
#x <- c(2, 6, 10, 14, 20, 28)/100
#ptox <- round(cbind(x/2, x/1.4, x)*2.4, 2)
#cat("SCENARIO ", scen_num, ":\n", sep="")
#source("one_scen.R")

#Scenario 4 - best schedule is 4
#rm(list=ls())
#source("setup_sampsize1.R")
#scen_num <- 4
#x <- c(2, 6, 10, 14, 20, 28)/100
#ptox <- round(cbind(x/3, x/1.6, x)*1.7, 2)
#ptox[3,] <- c(6, 10, 15)/100
#cat("SCENARIO ", scen_num, ":\n", sep="")
#source("one_scen.R")

#Scenario 5 - best schedule is 5
rm(list=ls())
source("setup_sampsize1.R")
scen_num <- 5
x <- c(2, 6, 10, 14, 20, 28)/100
ptox <- round(cbind(x/2, x/1.3, x)*1.3, 2)
ptox[,3] <- c(3, 8, 12, 15, 26, 36)/100
cat("SCENARIO ", scen_num, ":\n", sep="")
source("one_scen.R")
