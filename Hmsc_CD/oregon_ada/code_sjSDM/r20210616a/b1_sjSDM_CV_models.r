## Data set up for models on ADA - 5 fold cross validation

## Run on ADA

## Sets up multiple Cross validation on ADA. 

## 0. Contents #####
## 1. Gets data from github with criteria below
## 2. Creates k folds and sets up tuning grid, and results table for grid runs per k
## 3. Creates scaled data sets per k, runs, saves, evaluates model for each tune grid step
## 4. Writes results to csv. 
## 5. Averages results by k cv and writes to cv. 

## Individual species AUCs are not saved.

#  iter=170, nstep=800 (or 900. learning_rate = c(0.001, 0.002, 0.003) 


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
getwd() # always run sub from oregon_ada

library(dplyr)

resFolder <-"code_sjSDM/r20210616a/results"
if(!dir.exists(resFolder)) dir.create(resFolder, recursive = TRUE)

## Updated to new vars, also changes to elevation_m, canopy_height_m  to _f. 

# # model settings:
abund <- "qp"

spChoose <- "M1S1_v16"

device <- "gpu"
iter <- 170L
sampling <- 5000L

## Number of samples from tuning grid - random search
noSteps <- 1000

# no of CV folds
k <- 5

# timings... 
# models take approx to 0.28 mins run (Run time 00:36:11 for 125 models)
# noSteps * k * 0.28/60

noSteps * k

# Storage
# About ~600k per model .rds -- in GB
# noSteps * k * 600  / 1048576

# for testing on cpu
# device <- "cpu"
# iter <- 10L
# sampling <- 100L
# noSteps <- 2
# k <- 2

### 1. Get data from github #####

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",
                                       minimaprundate,"_",
                                       samtoolsfilter,"_", 
                                       samtoolsqual, 
                                       "_kelpie", 
                                       kelpierundate,
                                       "_",
                                       primer,
                                       "_vsearch97")

datFile <- paste0("sample_by_species_table_", 
                  samtoolsfilter, 
                  "_minimap2_",
                  minimaprundate,
                  "_kelpie",
                  kelpierundate,
                  "_FSL_qp.csv")

# file path:
fn <- file.path(gitHub, outputidxstatstabulatefolder, datFile)

# what file am i using?
basename(fn)

# when was it modified? - only if stored locally. 
file.mtime(fn)

# read complete data set
otuenv <- read.csv(fn, stringsAsFactors = FALSE, na.strings = "NA")

# keep OTUs with >= minocc incidences AND with presnece at both M1 or M2
minocc <- 6 # set to high number (e.g. 20) for testing

## get Species columns by M1 and M2, with minocc calculated per trap
## can choose below whether to include just species shared between M1 and M2
spM <- otuenv %>% 
  dplyr::filter(period == "S1") %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0) %>%
  filter(M1) %>% # CHOOOSE HERE FOR SINGLE. OR SHARED TRAP SPECIES GFROUP: filter(M1 & M2)
  tidyr::separate(col = OTU, into = c("ID", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID", "size"),
                  remove = FALSE, sep = "_") %>%
  #dplyr::filter(order == "Diptera")%>%
  dplyr::select(OTU)

nrow(spM)

# filter species here to those in sp.M1m2$OTU - already filtered for minocc
otu.qp.csv <- otuenv %>% 
  dplyr::filter(period == "S1" & trap == "M1") %>%
  dplyr::select(spM$OTU) ## 

# convert to presence/absence data
otu.pa.csv <- otu.qp.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
min(colSums(otu.pa.csv)) >= minocc # should be TRUE

# clean up
rm(datFile, gitHub, kelpierundate, minimaprundate, outputidxstatstabulatefolder, primer, samtoolsfilter, samtoolsqual, fn, spM)

# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates

## Load new version of vars
load("data/envVars.rdata")
# head(allVars)

env.vars <- otuenv %>% 
  dplyr::filter(period == "S1" & trap == "M1") %>%
  dplyr::select(trap, period, UTM_E, UTM_N, SiteName) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_")) %>%
  left_join(y = allVars, by = "SiteName") %>%
  mutate(lg_DistStream = log(DistStream + 0.001),
         lg_DistRoad = log(DistRoad + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001))

# str(env.vars)
# head(env.vars)
# cat(paste(colnames(env.vars), collapse = '", "'))

# old vars
# oldVars <- c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")
# 
# # new vars
# newvars <- c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")

# varsName <- "vars3"
# vars <- c("be10", "slope", "tri", "Nss", "Ess","ndmi_stdDev", "ndvi_p5", "ndvi_p50", "ndvi_p95", "ndmi_p5", "ndmi_p50", "ndmi_p95", "savi_p50", "LC08_045029_20180726_B1", "LC08_045029_20180726_B3", "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7", "LC08_045029_20180726_B10", "ndmi_stdDev_100m", "ndvi_p5_100m", "ndvi_p50_100m", "ndvi_p95_100m", "ndmi_p5_100m", "ndmi_p50_100m", "ndmi_p95_100m", "savi_p50_100m", "LC08_045029_20180726_B1_100m", "LC08_045029_20180726_B3_100m", "LC08_045029_20180726_B4_100m", "LC08_045029_20180726_B5_100m", "LC08_045029_20180726_B7_100m", "LC08_045029_20180726_B10_100m", "tpi250",  "tpi500", "tpi1k" , "ht", "ht.r250", "ht.r1k", "cov2_4.r250", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r1k", "mTopo", "cut.r1k.pt", "insideHJA", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "l_p25", "l_rumple")


# varsName <- "vars5"
# vars <- c("ht30", "gt4_r30", "gt4_500", "cut_r500", "cut_r250", "cut40_r500", "cut40_r250", "be30", "tri30", "slope30", "Nss30", "Ess30", "twi30", "tpi250", "tpi500", "tpi1k", "l_rumple", "insideHJA", "ndmi_stdDev_r100", "nbr_stdDev_r100", "ndvi_p5_r100", "ndvi_p50_r100", "ndvi_p95_r100", "ndmi_p5_r100", "ndmi_p50_r100", "ndmi_p95_r100", "LC08_045029_20180726_B1", "LC08_045029_20180726_B3", "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7", "LC08_045029_20180726_B10", "minT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m")

# varsName <- "vars6"
# vars <- c("ht30","gt4_r30","gt4_250",
#           "cut_r250","cut40_r1k","cut40_r250", 
#           "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k",
#           "l_Cover_2m_max", "l_Cover_4m_16m", "l_rumple",
#           "DistStream","DistRoad",
#           "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500", "ndvi_p5_r100","ndvi_p5_r500",
#           "ndvi_p50_r100","ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500",
#           "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10",
#           "minT_annual","precipitation_mm")

# load VIF vars and hard code below
# load(file.path(resFolder, "selected-vars.rdata"))
# cat(paste(vif.count$var, collapse = '", "'))

# varsName <- "vars7" # VIF5
# vars7 <- c("gt4_r30", "gt4_500", "cut_40msk", "cut_r1k", "cut_r250", "cut40_r1k", "cut40_r250", "tri30", "Nss30", "Ess30", "twi30", "tpi250", "tpi1k", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_rumple", "DistStream", "DistRoad", "insideHJA", "ndmi_stdDev_r100", "nbr_stdDev_r250", "ndvi_p5_r100", "ndvi_p95_r500", "ndmi_p95_r100", "ndmi_p95_r500", "LC08_045029_20180726_B3", "LC08_045029_20180726_B5", "LC08_045029_20180726_B10", "minT_annual")

varsName <- "vars8" # be30 replaces minT
vars <- c("Ess30", "DistRoad", "Nss30", "DistStream", "tri30", "LC08_045029_20180726_B5", "l_rumple", "cut_40msk", "insideHJA", "cut_msk", "cut40_r1k", "l_Cover_2m_4m", "l_Cover_4m_16m", "cut_r1k", "be30", "ndmi_p95_r100", "ndvi_p50_r500", "cut40_r250", "nbr_stdDev_r250", "precipitation_mm", "gt4_250", "tpi250", "cut_r250", "LC08_045029_20180726_B1", "ndvi_p50_r100", "gt4_r30", "twi30", "ndvi_p5_r100", "tpi1k", "ndvi_p5_r500", "l_p25", "LC08_045029_20180726_B10")


# varsName <- "vars6x" # selected from vars 6 - 25 highest performing variables
# vars <- c("minT_annual", "slope30", "cut_r250", "be30", "Ess30", "l_rumple", "l_Cover_4m_16m", "DistStream", "ndmi_stdDev_r500", "cut40_r250", "tpi1k", "DistRoad", "LC08_045029_20180726_B10", "Nss30", "tpi250", "l_Cover_2m_max", "ndvi_p5_r100", "ndmi_p50_r100", "precipitation_mm", "ht30", "cut40_r1k", "twi30", "LC08_045029_20180726_B5", "gt4_r30", "ndvi_p50_r500")

# setdiff(vars7, vars)
# intersect(vars7, vars)
# 
# check names 
all(vars %in% colnames(env.vars))

## Save model data
save(otu.pa.csv, otu.qp.csv, otuenv, env.vars,spChoose,
     k, minocc, noSteps, vars, varsName, abund, device, iter, sampling,
     file = file.path(resFolder, paste0("modelData_",abund,".rdata")))


### 2. Make testing and training k folds #####

# Make folds 
set.seed(100)

# make fold ids, length == to data set nrow, sample to randomise order
fold.id <- sample(rep_len(1:k, length.out = nrow(otu.pa.csv)))
table(fold.id)

### Create scaled data sets for each fold and run model
# i = 1

## Create tuning grid
# set variables
# formula.env = 'envDNN'
###CHECK THESE FROM OTHER SCRIPT
lambda.env = seq(0,1, length.out=6)	# .1
alpha.env = seq(0,1, length.out=6)		# .9
lambda.sp = seq(0,1, length.out=5)	# .1 ## Changed to 5
alpha.sp =  seq(0,1, length.out=5)	# .5 ## Changed to 5
hidden <- list(c(50L,50L,10L), c(25L,25L,10L))
hidden.ind = seq_along(hidden)
acti.sp = 'relu'
# drop = seq(0.1,0.4, length.out=4) # .3
drop = c(0.1,0.2) # .3
sample.bio = seq(0,1,length.out=11)
lr = c(0.001, 0.002)

## Make grid of priority tune parameters, choose, from these in sampling, then add lower priority parameters
tune.grid <- expand.grid(lambda.env = lambda.env, alpha.env= alpha.env, lambda.sp = lambda.sp,
                         alpha.sp = alpha.sp, lr= lr, hidden.ind = hidden.ind,
                         drop= drop, acti.sp = acti.sp, stringsAsFactors = FALSE)
head(tune.grid)

#  no models
# 4*4*7*7*2*3 *11
# data storage
# 4*7*7*2*3 * 642
# time
# 4*7*7*2*3*2 * .5 # 162 models - 1 hour

# get samples from tune grid
tune.rows <- sample(1:nrow(tune.grid), size = noSteps, replace = FALSE)


# add in lambda.bio and alpha.bio from random samples of sample.bio, and add results columns
tune.results <- data.frame(tr = 1:noSteps,
                           tune.grid[tune.rows,], 
                           lambda.bio = sample(sample.bio, size = noSteps, replace = TRUE),
                           alpha.bio = sample(sample.bio, size = noSteps, replace = TRUE),
                           loglike_sjSDM = NA, 
                           loss= NA, 
                           AUC.train = NA, 
                           AUC.test = NA,
                           ll.train = NA,
                           ll.test = NA,
                           nagel.train = NA,
                           nagel.test = NA,
                           plr.train = NA,
                           plr.test = NA,
                           tjur.train = NA,
                           tjur.test = NA,
                           cor.train = NA,
                           cor.test = NA,
                           auc.lt5.train = NA,
                           auc.lt5.test= NA)


str(tune.results)
# Add in k, to keep all cross validation runs. Average later.
tune.results <- cbind(tune.results[rep(seq(noSteps), k),], k = rep(1:k, each = noSteps))
rownames(tune.results) <- NULL
head(tune.results)


# clean up
rm(lambda.env, alpha.env, lambda.sp, alpha.sp, hidden.ind, drop, sample.bio, acti.sp)


# Choose pa or qp reponse data and family
if(abund == "pa") {
  Y <- otu.pa.csv
  family <- stats::binomial('probit') } else {
    if(abund ==  "qp") {
      Y <- otu.qp.csv
      family <- stats::binomial('probit') # check other family?
    } else stop("check abund")
  } 


### 3. Create scaled testing training data sets, run models and evaluate ####
# i= 1
for(i in 1:k){ # start fold loop 
  
  # select X data
  env.train <- env.vars[fold.id != i, vars]
  env.test <- env.vars[fold.id == i, vars]
  
  # any(sapply(env.train, function(x) any(is.na(x))))
  # any(sapply(env.test, function(x) any(is.na(x))))
  
  # select spatial data
  XY.train <- env.vars[fold.id != i, c("UTM_E", "UTM_N")]
  XY.test <- env.vars[fold.id == i, c("UTM_E", "UTM_N")]
  
  # select Y data
  s.otu.train <- as.matrix(Y[fold.id != i,])
  s.otu.test <- as.matrix(Y[fold.id == i,])
  
  
  factV <- colnames(env.train)[sapply(env.train, is.factor)]
  
  # scale X and spatial data for each fold
  # .. env data - without factors
  scale.env.train.all = dplyr::select(env.train, -any_of(factV)) %>% scale()
  # str(scale.env.train.all)
  
  # put factors back in 
  scale.env.train <- data.frame(scale.env.train.all, env.train[,factV, drop = F])
  
  # data frame of means and sd for scaling
  dd.env.scaler = data.frame(t(data.frame(env.mean = attr(scale.env.train.all, "scaled:center"), env.sd = attr(scale.env.train.all, "scaled:scale"))))
  
  scale.env.test = as.data.frame(do.call(rbind,apply(dplyr::select(env.test, -any_of(factV)), 1, function(x){(x-dd.env.scaler['env.mean',])/dd.env.scaler['env.sd',]} ) )) %>%
    tibble::add_column(env.test[,factV, drop = F])
  
  #dim(scale.env.train)
  #dim(scale.env.test)
  
  # .. spatial data
  XY.train.scale <- scale(XY.train)
  
  dd.xy.scaler = data.frame(t(data.frame(sp.mean = attr(XY.train.scale, "scaled:center"), 
                                         sp.sd = attr(XY.train.scale, "scaled:scale"))))
  
  XY.test.scale <- as.data.frame(do.call(rbind, 
                                         apply(XY.test, 1,function(x){(x-dd.xy.scaler['sp.mean',])/dd.xy.scaler['sp.sd',]})))
  
  rm(dd.xy.scaler, dd.env.scaler, scale.env.train.all)
  
  XY.train.scale <- data.frame(XY.train.scale)
  
  
  rm(env.test, env.train, XY.test, XY.train)
  # summary(scale.env.train)
  # summary(scale.env.test)

  ## Do tuning per k
  # j = 1 # j= 2
   
  
  for(j in seq_len(noSteps)){
    
    #lambda.bioN = sample(1:11,1)
    #alpha.bioN = sample(1:11,1)
  
    cat("\nk", i, "tune run", j)
    # print(tune.results[tune.results$k == i, ][j,])
    
    # subset this round of tuning parameters to make easier to insert in model specs
    tr <- subset(tune.results, k == i)[j,,drop = T]

    # do model    
    model.train <- sjSDM(
      
      Y = s.otu.train,
      
      env = DNN(data=scale.env.train, formula = ~.,
                hidden=hidden[[tr$hidden.ind]],
                lambda = tr$lambda.env,
                alpha = tr$alpha.env,
                activation = tr$acti.sp,
                dropout=tr$drop,
                bias=T),
      
      biotic = bioticStruct(lambda=tr$lambda.bio,
                            alpha=tr$alpha.bio, 
                            on_diag=F, inverse = FALSE),
      
      spatial = linear(data=XY.train.scale, ~0+UTM_E*UTM_N, 
                       lambda=tr$lambda.sp, 
                       alpha=tr$alpha.sp),
      
      learning_rate = tr$lr, # part of tuning grid
      iter = iter, 
      family = family, 
      sampling = sampling, # 150L, 5000L
      device = device
    )
    
    ## SAve each model
    saveRDS(model.train,
            file.path(resFolder,paste0("s-jSDM_tuning_model_", varsName, "_", k, "CV_", 
                                       i, "_tr_", j, "_", abund, ".rds")))
    
    # Do testing and save results in data frame
    tune.results$loglike_sjSDM[tune.results$k == i][j] <- logLik(model.train)
    tune.results$loss[tune.results$k == i][j] <-model.train$history[length(model.train$history)]
    
    
    ## Model evaluation - AUC, LL, r2 -- could put all this into same loop over species columns, then use colMeans, etc below.
    ## put togehter from different metric scripts. 
    
    for (pred in 1:2) {
      
      # pred = 1
      # 1 -> 'test'
      newdd = scale.env.test ; newsp = XY.test.scale; otudd = s.otu.test
      
      ## pred = 2 # training
      if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train}
      # Error in reticulate::py_is_null_xptr(fa) : 
      #   Cannot convert object to an environment: [type=double; target=ENVSXP].
      # if (pred==2) { newdd = scale.env.train; newsp = XY.train.scale; otudd = s.otu.train}
      
      # predict for all species = sites X columns
      pred.dd = apply(abind::abind(lapply(1:3, function(i) {
        predict(model.train, newdata=newdd, SP=newsp)}
        ),along = -1L), 2:3, mean)
      
      attr(pred.dd, 'dimnames') = NULL
      
      # convert observed to pa (if qp)
      otudd.pa = (otudd>0)*1
      # sum(colSums(otudd.pa)==0)
      
      # AUC for all species - if not present, then NA
      auc <- sapply(1:dim(otudd)[2], function(i) {
        
        tryCatch({
          as.numeric(pROC::roc(otudd.pa[,i], pred.dd[,i], direction = "<", quiet=T)$auc)},
          error = function(err){ return(NA)}
        )
      })
      
      ## Extra evaluation metrics
      # ll, nagel & plr for spp 
      rsq = data.frame(ll=rep(.1, length.out=ncol(pred.dd)), 
                       nagel=rep(.1, length.out=ncol(pred.dd)), 
                       plr=rep(.1, length.out=ncol(pred.dd)),
                       tjur = rep(NA, length.out=ncol(pred.dd)),
                       cor = rep(NA, length.out=ncol(pred.dd)))
      
      # m = 59
      for (m in 1:ncol(pred.dd)) { 
        
        p = pred.dd[,m]; y = otudd.pa[,m]
        
        loglikP = sum( log( p*y + (1-p)*(1-y) ) )
        loglikN = sum( log( mean(p)*y + (1-mean(p))*(1-y) ) )
        
        rsq$nagel[m] = (1-exp(2/length(p)*(loglikN-loglikP))) / (1-exp(2/length(p)*loglikN))
        rsq$ll[m] = loglikP
        
        tppp = sum(p*y)
        fppp = sum(p*(1-y))
        fapp = sum((1-p)*y) ### Weird behavious if fa object is here... check in sjSDM code... 
        tapp = sum((1-p)*(1-y))
        
        rsq$plr[m] = tppp/(tppp+fapp)/fppp*(fppp+tapp) # get NaN if species missing at all sites. OK> 
        
        tjur <- base::diff(tapply(p, y, mean, na.rm = T))
        rsq$tjur[m] <- ifelse(length(tjur) > 0, tjur, NA)
        rsq$cor[m] = suppressWarnings(cor(p, y)) # warning - with NaN when no presences in species in test set
        
      }
      
      
      # Tjur <- tryCatch(expr = mapply(function(x,y) {
      #   tj <- base::diff(tapply(x, y, mean, na.rm = T)) # difference of average predicted values at 1 and 0
      #   if(!length(tj) > 0) tj <- NA
      #   return(tj)
      # }, asplit(pred.dd,2), asplit(otudd.pa,2)), error = function(err){ return(NA)})
      # 
      
      ## add to data frame
      
      # if (pred==2) {tune.results$AUC.train[tune.results$k == i][j] = mean(auc, na.rm = TRUE)}
      # if (pred==1) {tune.results$AUC.test[tune.results$k == i][j] = mean(auc, na.rm = TRUE)}
      # 
      
      if (pred==2) {
        tune.results$AUC.train[tune.results$k == i][j] = mean(auc, na.rm = TRUE)
        tune.results$ll.train[tune.results$k == i][j] = mean(rsq$ll, na.rm = T)
        tune.results$nagel.train[tune.results$k == i][j] = mean(rsq$nagel, na.rm = T)
        tune.results$plr.train[tune.results$k == i][j]  = mean(rsq$plr, na.rm = T)
        tune.results$tjur.train[tune.results$k == i][j]  = mean(rsq$tjur, na.rm = T)
        tune.results$cor.train[tune.results$k == i][j]  = mean(rsq$cor, na.rm = T)
        tune.results$auc.lt5.train[tune.results$k == i][j]  = sum(auc < 0.5, na.rm = T)
        }
      
      if (pred==1) {
        tune.results$AUC.test[tune.results$k == i][j] = mean(auc, na.rm = TRUE)
        tune.results$ll.test[tune.results$k == i][j] = mean(rsq$ll, na.rm = T)
        tune.results$nagel.test[tune.results$k == i][j] = mean(rsq$nagel, na.rm = T)
        tune.results$plr.test[tune.results$k == i][j] = mean(rsq$plr, na.rm = T)
        tune.results$tjur.test[tune.results$k == i][j]  = mean(rsq$tjur, na.rm = T)
        tune.results$cor.test[tune.results$k == i][j]  = mean(rsq$cor, na.rm = T)
        tune.results$auc.lt5.test[tune.results$k == i][j]  = sum(auc < 0.5, na.rm = T)
        }
      
      
      
    } # end of evaluation loop
    
 rm(model.train)
    
  } # end of model loop
  
  

}

#### 4. WRite results to csv ####
write.table(tune.results,
            file = file.path(resFolder,paste0("manual_tuning_sjsdm_", varsName, "_", k, "CV_", spChoose, "_",
                                              abund, 
                                              "_min",
                                              minocc,
                                              "_nSteps",
                                              noSteps,
                                              ".csv")), row.names=F, sep=',')


head(tune.results)

### 5. Average AUC by tune runs ####
tune.mean <- tune.results %>%
  group_by(across(c(-loglike_sjSDM, -loss, -k, -contains(c("train", "test"))))) %>%
  #summarise(across(contains(c("train", "test"))), list(mean))
  summarise(across(contains(c("train", "test")), list(mean = mean)))%>%
  arrange(desc(AUC.test_mean))
   
head(data.frame(tune.mean))


write.table(tune.mean,
            file = file.path(resFolder,paste0("manual_tuning_sjsdm_", varsName, "_", k, "CV_", spChoose, 
                                              "_meanEVAL_", 
                                              abund, 
                                              "_min",
                                              minocc,
                                              "_nSteps",
                                              noSteps,
                                              ".csv")), row.names=F, sep=',')