#### ecocopula analysis

library(mvabund)
library(ecoCopula)
library(doParallel)
library(foreach)

# only on local
# getwd()
# wd <- here::here()
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/oregon_ada")
dir()
rm(list = ls())

# get data
# source("L1_read_data.r") # local
source("code/S1_read_data.r")
# otu.pa.csv is Y.train.pa; otu.qp.csv is Y.train.qp
rm(Y.train.qp, Y.train.pa, X.train)

## Make list of formulae to loop through, 
colnames(env.vars)
# [1] "clearcut"         "insideHJA"        "oldGrowthIndex"   "elevation_m"      "canopyHeight_m"   "precipitation_mm"
# [7] "minT_annual"      "maxT_annual"      "mean.NDVI"        "mean.EVI"         "mean.green"       "mean.wet"        
# [13] "mean.bright"      "l_p25"            "l_rumple"         "B1_median"        "B2_median"        "B3_median"       
# [19] "B4_median"        "B5_median"        "B6_median"        "B7_median"        "B10_median"       "B11_median"      
# [25] "lg_DistStream"    "lg_DistRoad"      "lg_YrsDisturb"    "lg_cover2m_max"   "lg_cover2m_4m"    "lg_cover4m_16m"  
# [31] "dem500"           "tri.pt"

fm <- list(as.formula(otu.pa ~ dem500 + tri.pt),
           as.formula(otu.pa ~ dem500 + tri.pt + insideHJA + lg_YrsDisturb),
           as.formula(otu.pa ~ dem500 * oldGrowthIndex + tri.pt + insideHJA + lg_YrsDisturb),
           as.formula(otu.pa ~ dem500 * insideHJA + tri.pt + oldGrowthIndex + lg_YrsDisturb),
           as.formula(otu.pa ~ dem500 + tri.pt + insideHJA + oldGrowthIndex + lg_YrsDisturb + mean.NDVI + lg_DistRoad + canopyHeight_m),
           as.formula(otu.pa ~ dem500 * insideHJA + dem500 * oldGrowthIndex + tri.pt + lg_YrsDisturb + mean.NDVI + lg_DistRoad)
)

fm

# make a pa matrix as mvabund object
# otu.pa <- mvabund(otu.pa.csv)
# is.mvabund(otu.pa)

## ALternative data specification as list with btoh response and predictors as list elements 
# (foreach wasn't finding otu.pa)
dataN <- c(list(otu.pa = mvabund(otu.pa.csv)), as.list(env.vars))
str(dataN)

## set up parallel
nCores <- 16
# nCores <- 3

cl <- makeCluster(nCores)
registerDoParallel(cl)

nBoot <- 500


# i = 1

modList <- foreach(i = seq_along(fm),
                   .combine = list,
                   .multicombine = TRUE, # ,.export = c("otu.pa", "env.vars", "S.train")
                   .noexport = c("otu.pa.csv","otu.qp.csv", "P")
                   ) %dopar% {
                     
    # do glm model
    mod <- mvabund::manyglm(fm[[i]], family = "negative.binomial", data = dataN)
    
    # do ordination
    mod.ord <- ecoCopula::cord(mod)
    
    # Site scores: join ordination axes to data frame of predictors, coordinates
    site_res <- data.frame(mod.ord$scores, env.vars, S.train[,c("UTM_E", "UTM_N")])
    # Species scores
    sp_res <- data.frame(mod.ord$loadings, species = colnames(dataN$otu.pa))
  
    # do anova
    summ <- mvabund::summary.manyglm(mod, nBoot = nBoot, test = "LR")
    summ
    
    list(mod = mod, ord = mod.ord, site = site_res, sp=sp_res, summ = summ)
    
                   }

stopCluster(cl)

str(modList, max.level =2)
modList[[3]]$summ

save(modList, file = "results/ecocopula/ecocopula_modList_pilot.rdata")
