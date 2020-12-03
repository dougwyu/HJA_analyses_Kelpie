### Run pipeline ####


## Only local: 
# setwd("J:/UEA/Oregon/Oproject/oregon_ada")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/oregon_ada")
dir()
rm(list = ls())


## 0. Create results folder to hold code and models #####
resFolder <- "res20201127_01"
modFolder <- file.path("results", resFolder, "models")


## 5 Compute model fit ####
nfolds = 5
fs <- list.files(modFolder, "^models_thin.*", full.names = TRUE)
source(file.path("results", resFolder, "code/S5_compute_model_fit.r"))


# i = 2