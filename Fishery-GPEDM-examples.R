rm(list=ls())
library(TMB)
library(Matrix)
source("Fishery-GPEDM-Rfunctions.R")
############################################################################################################################################
# 1. Example for simulated data from Pella-Tomlinson production model (one-way trip increasing harvest rate scenario)
load("PellaTomlinson_increasingharvest_example.RData")
cpue <- B.all[,1] #simulated biomass/abundance from PT production dynamics
catch <- C.all[,1] #simulated catch from PT production dynamics
fgpedm <- fisheryGPEDMfit_example_fn(cpue=cpue,catch=catch,E=4,method="rprop",split=0.5,transformY="original")
plot(cpue,type="l",xlab="Time",ylab="Abundance index",main="Pella-Tomlinson model")
plot(catch,type="l",xlab="Time",ylab="Catch",main="Pella-Tomlinson model")
plot(fgpedm$obs,fgpedm$looPred,xlab="Observed abundance",ylab="LOO prediction",main="Fishery-GPEDM"); abline(0,1)
paste0("Out-of-sample LOO R2:",fgpedm$looR2)

############################################################################################################################################
# 2. Example for simulated data from prey-predator dynamics model (prey species harvested and one-way trip increasing harvest rate scenario)
load("PreyPredator_increasingharvest_example.RData")
cpue <- prey.all[,1] #simulated biomass/abundance from prey-predator dynamics
catch <- C.all[,1] #simulated catch from from prey-predator dynamics
fgpedm <- fisheryGPEDMfit_example_fn(cpue=cpue,catch=catch,E=5,method="rprop",split=0.5,transformY="log") #log transform of CPUE
plot(cpue,type="l",xlab="Time",ylab="Abundance index (prey)",main="Prey-predator dynamics with prey harvested")
plot(catch,type="l",xlab="Time",ylab="Catch",main="Prey-predator dynamics with prey harvested")
plot(fgpedm$obs,fgpedm$looPred,xlab="Observed abundance",ylab="LOO prediction",main="Fishery-GPEDM"); abline(0,1)
paste0("Out-of-sample LOO R2:",fgpedm$looR2)

############################################################################################################################################
# 3. Example for Gulf of Mexico Brown Shrimp empirical dynamics
# brownshrimp_cpue: Yearly Brown Shrimp SEAMAP CPUE averaged across the Gulf of Mexico from 1987 to 2019 (tails per tow)
# brownshrimp_catch: Yearly Brown Shrimp landings across the Gulf of Mexico from 1987 to 2019 (pounds)
# Note: The empirical data here are used only for demonstrating Fishery-centric GPEDM method and for reproducibility of the paper.
#       Further access to data and details please refer to the Southeast Area Monitoring and Assessment Program (SEAMAP) (1987-2019).
#       Data are available upon reasonable request Michell Masi and Molly Stevens and with permission of SEAMAP.
load("GulfMexicoBrownShrimp_meancpue_totalcatch.RData")
cpue <- brownshrimp_cpue # Brown Shrimp catch per unit effort (CPUE; tails per tow)
catch <- brownshrimp_catch*0.000001 # Brown Shrimp landings (million pounds)
fgpedm <- fisheryGPEDMfit_example_fn(cpue=cpue,catch=catch,E=5,method="rprop",split=0.9,transformY="logdifference") #log-difference transform of CPUE
plot(cpue,type="l",xlab="Time",ylab="CPUE (tails per tow)",main="Gulf of Mexico Brown Shrimp")
plot(catch,type="l",xlab="Time",ylab="Landings (million pounds)",main="Gulf of Mexico Brown Shrimp")
plot(fgpedm$obs,fgpedm$looPred,xlab="Observed CPUE",ylab="LOO prediction",main="Fishery-GPEDM"); abline(0,1)
paste0("Out-of-sample LOO R2:",fgpedm$looR2)


