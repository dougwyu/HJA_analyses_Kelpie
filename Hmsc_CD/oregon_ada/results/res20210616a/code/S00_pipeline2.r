### Run pipeline ####


## Only local: 
# getwd()
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()

## copy code into code folder, edit the name of results folder appropriately and run this script
resFolder <- file.path("results", paste0("res", "20210616a"))

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
options(echo=TRUE) # if you want see commands in output file
# setwd("~/oregon_ada") # wd already set by starting in oregon_ada on ADA 
getwd()

library(Hmsc)
packageVersion("Hmsc")

## 0. Create results folder to hold code and models #####

print(resFolder)
dir.create(file.path(resFolder, "code"))
dir.create(file.path(resFolder, "data"))
dir.create(file.path(resFolder, "models"))

modFolder <- file.path(resFolder, "models")
modFolder

## 1. Read data #### 
# source(file.path(resFolder,"code/S1_read_data.r"))

# save data
# save(S.train, X.train, Y.train.pa, Y.train.qp, P, file = file.path(resFolder, "data", "allData.Rdata"))

## 2. define models ####
# source(file.path(resFolder,"code/S2_define_models.r"))
# save(models, modelnames, file = file.path(modFolder, "unfitted_models.rdata"))

## 3 Fit models #####
# uses models list from 2.
# samples_list = c(5,250,250,250,250,250)
# thin_list = c(1,1,10,100,1000,10000)

samples_list = c(50,200)
thin_list = c(1,10)

## samples_list = c(5,100)
## thin_list = c(1,5)
nChains <- 4

## TESTING
# samples_list = c(5)
# thin_list = c(1)
# nChains <- 2

# iterations per chain
(samples_list * thin_list) + ceiling(0.5*thin_list*samples_list)

# source(file.path(resFolder,"code/S3_fit_models.r"))
# models saves in S3
# rm(samples_list, thin_list, nChains)

### 4 Convergence ####
# fs is the list of model filenames
# fs <- list.files(modFolder, "models_thin.*", full.names = TRUE)
# source(file.path(resFolder,"code/S4_evaluate_convergence.r"))
#save results
# save(beta, file = file.path(resFolder, "beta_list.rdata"))

# vioP <- require(vioplot)
# pdf(file = file.path(resFolder, "MCMC_convergence%03d.pdf"))
# par(mfrow=c(2,1), mar = c(10,3,1,1))
# # remove NAs coming through from try errors first:
# if(vioP) sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2) else sapply(beta[!is.na(beta)], boxplot, las = 2)
# dev.off()

## 5 Compute model fit ####
nfolds = 5
# fs <- list.files(modFolder, "^models_thin.*", full.names = TRUE)
# print(fs)
# 1] "results/res20210616a/models/models_thin_1_samples_50_chains_4.Rdata"  
# [2] "results/res20210616a/models/models_thin_10_samples_200_chains_4.Rdata"

fs <- "results/res20210616a/models/models_thin_10_samples_200_chains_4.Rdata"
source(file.path(resFolder, "code/S5_compute_model_fit.r"))


