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

# env.vars already has topo data in it

## Make list of formulae to loop through, 
colnames(env.vars)
# [1] "clearcut"         "insideHJA"        "oldGrowthIndex"   "elevation_m"      "canopyHeight_m"   "precipitation_mm"
# [7] "minT_annual"      "maxT_annual"      "mean.NDVI"        "mean.EVI"         "mean.green"       "mean.wet"        
# [13] "mean.bright"      "l_p25"            "l_rumple"         "B1_median"        "B2_median"        "B3_median"       
# [19] "B4_median"        "B5_median"        "B6_median"        "B7_median"        "B10_median"       "B11_median"      
# [25] "lg_DistStream"    "lg_DistRoad"      "lg_YrsDisturb"    "lg_cover2m_max"   "lg_cover2m_4m"    "lg_cover4m_16m"  
# [31] "be10"             "tri"              "slope"            "aspect"           "Nss"              "Ess"             
# [37] "twi"              "ht"               "ht.r250"          "ht.r500"          "ht.r1k"           "cov2_4"          
# [43] "cov2_4.r250"      "cov2_4.r500"      "cov2_4.r1k"       "cov4_16"          "cov4_16.r250"     "cov4_16.r500"    
# [49] "cov4_16.r1k"      "be500"            "mTopo"            "cut.r1k.pt"

fm <- list(as.formula(otu.pa ~  be500 + tri + insideHJA + lg_YrsDisturb),
           as.formula(otu.pa ~ be500 * oldGrowthIndex + tri + insideHJA + lg_YrsDisturb),
           as.formula(otu.pa ~ be500 * insideHJA + tri + oldGrowthIndex + lg_YrsDisturb),
           as.formula(otu.pa ~ be500 + tri + insideHJA + oldGrowthIndex + lg_YrsDisturb + mean.NDVI + lg_DistRoad + canopyHeight_m),
           as.formula(otu.pa ~ be500 * insideHJA + be500 * oldGrowthIndex + tri + lg_YrsDisturb + mean.NDVI + lg_DistRoad),
           as.formula(otu.pa ~ be10 + Nss + Ess + ht + cov2_4 + cov4_16 + ht.r1k + cov2_4.r1k + cov4_16.r1k + mTopo),
           as.formula(otu.pa ~ be10 + Nss + Ess + ht + cov2_4 + cov4_16 + ht.r500 + cov2_4.r500 + cov4_16.r500 + mTopo),
           as.formula(otu.pa ~ be10 + Nss + Ess + ht + cov2_4 + cov4_16 + ht.r250 + cov2_4.r250 + cov4_16.r250 + mTopo),
           as.formula(otu.pa ~ be500 + slope + twi + ht + cov2_4 + cov4_16 + ht.r500 + cov2_4.r500 + cov4_16.r500 + mTopo)
)

fm

# make a pa matrix as mvabund object
# otu.pa <- mvabund(otu.pa.csv)
# is.mvabund(otu.pa)

## ALternative data specification as list with btoh response and predictors as list elements 
# (foreach wasn't finding otu.pa)
dataN <- c(list(otu.pa = mvabund(otu.pa.csv)), as.list(env.vars))
# str(dataN)

# check all terms are in data object
unique(unlist(sapply(fm, all.vars)))
all(unique(unlist(sapply(fm, all.vars))) %in% names(dataN))

## set up parallel
nCores <- 12
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

dir.exists("results/ecocopula")
