### b3_sjSDM_prediction.r

## Run on ADA

## Predictions across landscape (on raster stack of predictors)
## 1. Create full model (from best tuning)
## 2. Predict across new data from raster stack (export data at regular grid of xy coords)
## convert to rasters and stack up
## 3. Do plots


#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada

library(dplyr)


resFolder <-"code_sjSDM/r20210406a/results"
abund <- "pa"

## load model data 
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
vars
varsName
head(env.vars)


## subset of vars for testing
vars <- c("be10", "Nss", "ndmi_stdDev", "ndvi_p5","ndmi_p50","tpi500","cov2_4.r1k","cut.r1k.pt","insideHJA")

## update for testing
# device <- "cpu"
# iter= 10
# sampling = 100

# formula.env = 'envDNN'
hidden <- list(c(50L,50L,10L), c(25L,25L,10L))

## get best tune run
res <- read.csv(file.path(resFolder,paste0("manual_tuning_sjsdm_", varsName, "_", k, "CV_M1S1_meanEVAL_", 
                                           abund, 
                                           "_min",
                                           minocc,
                                           "_nSteps",
                                           noSteps,
                                           ".csv")))

head(res)
res.best <- res[which.max(res$AUC.test_mean),,drop = T]
res.best
rm(res)


# Choose pa or qp reponse data and family
if(abund == "pa") {
  Y <- otu.pa.csv
  family <- stats::binomial('probit') } else {
    if(abund ==  "qp") {
      Y <- otu.qp.csv
      family <- stats::binomial('probit') # check other family?
    } else stop("check abund")
  } 


# select X data
scale.env <- env.vars %>%
  dplyr::select(all_of(vars))%>%
  dplyr::mutate(across(where(is.numeric), scale))

any(sapply(scale.env, function(x) any(is.na(x))))

# spatial data
scale.XY <- env.vars %>%
  select(c("UTM_E", "UTM_N")) %>%
  mutate(across(everything(), scale))

# do model - spatial and non spatial
modFull_sp <- sjSDM(
  
  Y = as.matrix(Y),
  
  env = DNN(data = scale.env, formula = ~.,
            hidden=hidden[[res.best$hidden.ind]],
            lambda = res.best$lambda.env,
            alpha = res.best$alpha.env,
            activation = res.best$acti.sp,
            dropout=res.best$drop,
            bias=T),
  
  biotic = bioticStruct(lambda=res.best$lambda.bio,
                        alpha=res.best$alpha.bio, 
                        on_diag=F, inverse = FALSE),
  
  spatial = linear(data=scale.XY, ~0+UTM_E*UTM_N,
                   lambda=res.best$lambda.sp,
                   alpha=res.best$alpha.sp),
  
  learning_rate = res.best$lr,
  iter = iter, 
  family = family, 
  sampling = sampling, # 150L, 5000L
  device = device
)

## Save model
saveRDS(modFull_sp, file.path(resFolder,paste0("s-jSDM_fullMod_sp_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rds")))


### Do prediction

# load newData
load(file.path(resFolder, paste0("predData.rdata")))

## Subset for vars

newData <- subset(allData, select = vars)

pred.sp = predict(modFull_sp, newdata = newData, SP = newSP)




