## Fishery-GPEDM
R and TMB examples for fishery-centric Gaussian Process Empirical Dynamic Modelling (GP-EDM).

# Code: 
1. "Fishery-GPEDM-examples.R" demonstrates examples using simulated data and empirical Brown Shrimp data.
2. "Fishery-GPEDM-simulations.R" contains parameters and model structure for simulating Pella-Tomlinson and prey-predatory harvesting dyanmics.
3. "gpEDMtmb_v5_linearscaling.cpp" is the TMB package template written in C used for "Fishery-GPEDM-Rfunctions.R".
4.  "Fishery-GPEDM-Rfunctions.R" contains R functions for GPEDM.

# Data:
1. "PellaTomlinson_increasingharvest_example.RData" contains simulations of Pella-Tomlinson model.
2. "PreyPredator_increasingharvest_example.RData" contains simulations of Pella-Tomlinson model.
3. "GulfMexicoBrownShrimp_meancpue_totalcatch.RData" contains empirical Gulf of Mexico Brown Shrimp data. Note that, the empirical data here are used only for demo. Further access to data and details please refer to the Southeast Area Monitoring and Assessment Program (SEAMAP) (1987-2019). Data are available upon reasonable request Michell Masi and Molly Stevens and with permission of SEAMAP.
