
#### ecocopula analysis - PLOTS

library(mvabund)
library(ecoCopula)
library(doParallel)
library(foreach)

getwd()
wd <- here::here()
setwd(wd)

# get data
source("Hmsc_CD/local/L1_read_data.r")
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
           as.formula(otu.pa ~ dem500 * insideHJA + dem500 * oldGrowthIndex + tri.pt + lg_YrsDisturb + mean.NDVI + lg_DistRoad),
           as.formula(otu.pa ~ dem500 * insideHJA + dem500 * oldGrowthIndex + tri.pt * oldGrowthIndex + lg_YrsDisturb + mean.NDVI + lg_DistRoad)
)

fm

# make a pa matrix as mvabund object
# otu.pa <- mvabund(otu.pa.csv)
# is.mvabund(otu.pa)

## ALternative data specification as list with btoh response and predictors as list elements 
# (foreach wasn't finding otu.pa)
dataN <- c(list(otu.pa = mvabund(otu.pa.csv)), as.list(env.vars))
str(dataN)


# check all terms are in data object
unique(unlist(sapply(fm, all.vars)))
all(unique(unlist(sapply(fm, all.vars))) %in% names(dataN))

## set up parallel
# nCores <- 16
nCores <- 2

cl <- makeCluster(nCores)
registerDoParallel(cl)

# nBoot <- 10

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
  
    # don't do this locally
    # do anova
    # summ <- mvabund::summary.manyglm(mod, nBoot = nBoot, test = "LR")
    # summ
    
    # list(mod = mod, ord = mod.ord, site = site_res, sp=sp_res, summ = summ)
    list(mod = mod, ord = mod.ord, site = site_res, sp=sp_res)
    
                   }

stopCluster(cl)

cor.preds <- colnames(env.vars)[sapply(env.vars, is.numeric)]

str(modList, max.level =2)
save(modList, fm, cor.preds, file = "Hmsc_CD/local/ecocopula_modList_pilot.rdata")

