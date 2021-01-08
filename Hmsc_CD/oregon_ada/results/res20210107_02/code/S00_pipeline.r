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
resFolder <- paste0("res", format(Sys.Date(), "%Y%m%d"))
resFolder <- paste(resFolder, sprintf("%02d", sum(grepl(resFolder, list.dirs("results", recursive= F)))+1), sep="_")
print(resFolder)
dir.create(file.path("results", resFolder))
dir.create(file.path("results", resFolder, "code"))
dir.create(file.path("results", resFolder, "data"))
dir.create(file.path("results", resFolder, "models"))

modFolder <- file.path("results", resFolder, "models")
modFolder

# copy scripts, etc
rScripts <- list.files("code", pattern = "S[1-5].*\\.[rR]$|S00_.*", recursive = TRUE, full.names = TRUE)
rScripts
to <- file.path("results", resFolder, rScripts)
file.copy(rScripts, to)

rm(to, rScripts)

## 1. Read data #### 
source(file.path("results", resFolder,"code/S1_read_data_vif.r"))
# Save in specific data folder

# ## Set working director to current results
# setwd(paste0("~/oregon_ada/results/", resFolder))
# # setwd(paste0("./results/", resFolder)) # local only
# getwd()
# dir()

# save data
save(S.train, X.train, Y.train.pa, Y.train.qp, P, file = file.path("results", resFolder, "data", "allData.Rdata"))

## 2. define models ####
source(file.path("results", resFolder,"code/S2_define_models.r"))
save(models, modelnames, file = file.path(modFolder, "unfitted_models.rdata"))

## 3 Fit models #####
# uses models list from 2.
# samples_list = c(5,250,250,250,250,250)
# thin_list = c(1,1,10,100,1000,10000)
samples_list = c(5,50)
thin_list = c(1,5)
# samples_list = c(5)
# thin_list = c(1)
# iterations per chain
(samples_list * thin_list) + ceiling(0.5*thin_list*samples_list)
nChains <- 4
source(file.path("results", resFolder,"code/S3_fit_models.r"))
# models saves in S3
rm(samples_list, thin_list, nChains)

### 4 Convergence ####
# fs is the list of model filenames
fs <- list.files(modFolder, "models_thin.*", full.names = TRUE)
source(file.path("results", resFolder,"code/S4_evaluate_convergence.r"))
#save results
save(beta, file = file.path("results", resFolder, "beta_list.rdata"))

vioP <- require(vioplot)
pdf(file = file.path("results", resFolder, "MCMC_convergence%03d.pdf"))
par(mfrow=c(2,1), mar = c(10,3,1,1))
# remove NAs coming through from try errors first:
if(vioP) sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2) else sapply(beta[!is.na(beta)], boxplot, las = 2)
dev.off()

## 5 Compute model fit ####
nfolds = 5
fs <- list.files(modFolder, "^models_thin.*", full.names = TRUE)
source(file.path("results", resFolder, "code/S5_compute_model_fit.r"))


