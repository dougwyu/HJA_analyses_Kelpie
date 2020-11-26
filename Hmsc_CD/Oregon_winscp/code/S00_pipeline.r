### Run pipeline ####


## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/Oregon_winscp")
dir()
rm(list = ls())


## 0. Create results folder to hold code and models #####
resFolder <- paste0("res", format(Sys.Date(), "%d%m%Y"))
resFolder <- paste(resFolder, sprintf("%02d", sum(grepl(resFolder, list.dirs("results", recursive= F)))+1), sep="_")
dir.create(file.path("results", resFolder))
dir.create(file.path("results", resFolder, "code"))
dir.create(file.path("results", resFolder, "data"))
file.path("results", resFolder, "models")
dir.create(modFolder)

modFolder <- file.path(resFolder, "models")
# resFolder <- "res26112020_01"

# copy scripts, etc
rScripts <- list.files("code", pattern = "S[1-5].*\\.[rR]$", recursive = TRUE, full.names = TRUE)
to <- file.path("results", resFolder, rScripts)
file.copy(rScripts, to)

rm(to, rScripts)

## 1. Read data #### 
source("code/S1_read_data.r")
# Save in specific data folder

## Set working director to current results
setwd(paste0("~/Oregon_winscp/results/", resFolder))
# setwd(paste0("./results/", resFolder)) # local only
getwd()
dir()

# save data
save(SXY.train, S.train, X.train, Y.train, P, file = file.path("data", "allData.Rdata"))

## 2. define models ####
source("code/S2_define_models.r")
save(models, modelnames, file = file.path("models", "unfitted_models.rdata"))

## 3 Fit models #####
# samples_list = c(5,250,250,250,250,250)
# thin_list = c(1,1,10,100,1000,10000)
# samples_list = c(5,100,250)
# thin_list = c(1,5,10)
samples_list = c(2)
thin_list = c(1)
nChains <- 2
source("code/S3_fit_models.r")
# models saves in S3

### 4 Convergence ####
# fs is the list of model filenames
fs <- list.files("models", "models_thin.*", full.names = TRUE)
source("code/S4_evaluate_convergence.r")
#save results
save(beta, file = "beta_list.rdata")

pdf(file = "MCMC_convergence%03d.pdf")
par(mfrow=c(2,1), mar = c(10,3,1,1))
# remove NAs coming through from try errors first:
if(vioP) sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2) else sapply(beta[!is.na(beta)], boxplot, las = 2)
dev.off()

## 5 Compute model fit ####
fs <- list.files("models", "^models_thin.*", full.names = TRUE)
source("code/S5_compute_model_fit.r")

