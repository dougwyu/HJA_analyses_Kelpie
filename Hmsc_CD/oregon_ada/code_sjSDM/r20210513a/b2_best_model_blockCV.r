## Data set up for models on ADA - 5 fold cross validation

## Run on ADA

## Sets up multiple Cross validation on ADA. 

## 0. Contents #####
## 1. Gets data from github with criteria below
## 2. Creates k folds and best tuning parameters for final model
## 3. Creates scaled data sets per k, runs, saves, evaluates model 
## 4. Writes results to csv. 
## 5. Averages results by k cv and writes to cv. 

## Individual species AUCs are saved.


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
library(blockCV)
library(raster)
library(sf)

resFolder <-"code_sjSDM/r20210513a/results"
abund <- "pa"
gis <- "data/gis"

## load model data 
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
# if(!dir.exists(resFolder)) dir.create(resFolder, recursive = TRUE)

alpha <- 1
lambda <- 0.001
lr <- 0.003

# for testing on cpu
# device <- "cpu"
# iter <- 10L
# sampling <- 100L
# k <- 2

## Load spatial data for block CV, poitns and template raster

## load template rasters
load(file.path(gis, "templateRaster.rdata")) ## r, indNA
# maks NAs
r.aoi.pred <- raster::mask(r.aoi.pred, r,maskvalue = 0)

## Load sample points
## Load sample site points
load(file.path(gis, "sample_sites.rdata"))

## filter points to data set
sp.sf <- subset(xy.utm, SiteName %in% env.vars$SiteName)


### 2. Make testing and training k folds with block or buffer design #####

# Make folds 
set.seed(100)

## Random k fold ###
# make fold ids, length == to data set nrow, sample to randomise order
# fold.id <- sample(rep_len(1:k, length.out = nrow(otu.pa.csv)))
# table(fold.id)
# k <- length(unique(fold.id))

## Block design ####
# spatial blocking by specified range with random assignment
# sb <- blockCV::spatialBlock(speciesData = sp.sf,
#                             rasterLayer = r.aoi.pred,
#                             species= NULL,
#                             theRange = 3000, # size of the blocks
#                             k = 5,
#                             selection = "checkerboard", #random attempts an even spread of presences - per species
#                             iteration = 100, # find evenly dispersed folds
#                             numLimit = 0,
#                             showBlocks = TRUE,
#                             biomod2Format = FALSE,
#                             xOffset = 0, # shift the blocks horizontally
#                             yOffset = 0)
# fold.id <- sb$foldID
# k <- length(unique(fold.id))

## Buffered design
bf <- buffering(speciesData = sp.sf,
                 theRange = 2000,
                 species = NULL, # if presnece absnece - then both are used
                 spDataType = "PA", # presence-absence  data type
                 progress = FALSE)

# str(bf, max.level = 1)
# bf$folds[[1]]

folds.train <- lapply(bf$folds, function(x) 1:nrow(sp.sf) %in% x[[1]])
folds.test <- lapply(bf$folds, function(x) 1:nrow(sp.sf) %in% x[[2]])

lengths(folds.train)
lengths(folds.test)

k <- length(folds.train)
# testing: k <- 3

# add in lambda.bio and alpha.bio from random samples of sample.bio, and add results columns
eval.results <- data.frame(k = 1:k,
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

# clean up

# Choose pa or qp reponse data and family
if(abund == "pa") {
  Y <- otu.pa.csv
  family <- stats::binomial('probit') } else {
    if(abund ==  "qp") {
      Y <- otu.qp.csv
      family <- stats::poisson('log')
    } else stop("check abund")
  } 


# species results
sp.res.train <- vector(length = k, mode = "list")
sp.res.test <- vector(length = k, mode = "list")

### 3. Create scaled testing training data sets, run models and evaluate ####
# i= 1

for(i in 1:k){
  
  # select X data
  env.train <- env.vars[folds.train[[i]], vars]
  env.test <- env.vars[folds.test[[i]], vars]
  
  # select spatial data
  XY.train <- env.vars[folds.train[[i]], c("UTM_E", "UTM_N")]
  XY.test <- env.vars[folds.test[[i]], c("UTM_E", "UTM_N")]
  
  # select Y data
  s.otu.train <- as.matrix(Y[folds.train[[i]],])
  s.otu.test <- as.matrix(Y[folds.test[[i]],])
  
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
  
  
  cat(paste("\n run k", i))
 
  # do model    
  model.train <- sjSDM::sjSDM(
    
    Y = s.otu.train,
    
    env = sjSDM::linear(
      data=scale.env.train,
      formula = ~.,
      lambda = lambda,
      alpha = alpha
      ),
    
    biotic = bioticStruct(lambda=lambda,
                          alpha=alpha., 
                          on_diag=F, inverse = FALSE),
    
    spatial = linear(data=XY.train.scale, ~0+UTM_E*UTM_N, 
                     lambda=lambda, 
                     alpha=alpha),
    
    learning_rate = lr, # 0.003 recommended for high species number 
    iter = iter, 
    family = family, 
    sampling = sampling, # 150L, 5000L
    device = device
  )
  
  ## SAve each model
  saveRDS(model.train,
          file.path(resFolder,paste0("s-jSDM_blockCV_model_", varsName, "_", k, "CV_", 
                                     i, "_", abund, ".rds")))
  
  # Do testing and save results in data frame
  eval.results$loglike_sjSDM[eval.results$k == i] <- logLik(model.train)
  eval.results$loss[eval.results$k == i] <-model.train$history[length(model.train$history)]
  
  

  ## Model evaluation - AUC, LL, r2
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
      predict(model.train, newdata=newdd, SP=newsp, type = "link")}
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
      rsq$cor[m] = suppressWarnings(cor(p, y)) # just when speceis missing in test data
      
    }
    
    
    # Tjur <- tryCatch(expr = mapply(function(x,y) {
    #   tj <- base::diff(tapply(x, y, mean, na.rm = T)) # difference of average predicted values at 1 and 0
    #   if(!length(tj) > 0) tj <- NA
    #   return(tj)
    # }, asplit(pred.dd,2), asplit(otudd.pa,2)), error = function(err){ return(NA)})
    # 
    
    ## add to data frame

    if (pred==2) {
      eval.results$AUC.train[eval.results$k == i] = mean(auc, na.rm = TRUE)
      eval.results$ll.train[eval.results$k == i] = mean(rsq$ll, na.rm = T)
      eval.results$nagel.train[eval.results$k == i] = mean(rsq$nagel, na.rm = T)
      eval.results$plr.train[eval.results$k == i]  = mean(rsq$plr, na.rm = T)
      eval.results$tjur.train[eval.results$k == i]  = mean(rsq$tjur, na.rm = T)
      eval.results$cor.train[eval.results$k == i]  = mean(rsq$cor, na.rm = T)
      eval.results$auc.lt5.train[eval.results$k == i]  = sum(auc < 0.5, na.rm = T)
    }
    
    if (pred==1) {
      eval.results$AUC.test[eval.results$k == i] = mean(auc, na.rm = TRUE)
      eval.results$ll.test[eval.results$k == i] = mean(rsq$ll, na.rm = T)
      eval.results$nagel.test[eval.results$k == i] = mean(rsq$nagel, na.rm = T)
      eval.results$plr.test[eval.results$k == i] = mean(rsq$plr, na.rm = T)
      eval.results$tjur.test[eval.results$k == i]  = mean(rsq$tjur, na.rm = T)
      eval.results$cor.test[eval.results$k == i]  = mean(rsq$cor, na.rm = T)
      eval.results$auc.lt5.test[eval.results$k == i]  = sum(auc < 0.5, na.rm = T)
    }
    
    ## species results
    
    if (pred==1) {
      rsq$auc <- auc
      sp.res.test[[i]] <- rsq
    }
    
    if (pred==2) {
      rsq$auc <- auc
      sp.res.train[[i]] <- rsq
    }
    
    
  } # end of evaluation loop
  
  rm(model.train)
  
} # end of model loop

eval.results
apply(eval.results, 2, mean)

### 5. Average AUC by tune runs ####
# eval.mean <- eval.results %>%
#   group_by(across(c(-loglike_sjSDM, -loss, -k, -contains(c("train", "test"))))) %>%
#   #summarise(across(contains(c("train", "test"))), list(mean))
#   summarise(across(contains(c("train", "test")), list(mean = mean))) 


## species means
# str(sp.res.test, max.level = 2)

names <- colnames(sp.res.test[[1]])

sp.mn.test <- lapply(names, function(z){
  
  rowMeans(do.call(cbind, lapply(sp.res.test, function(x) x[[z]])), na.rm = T)
  
})

names(sp.mn.test) <- names
sp.mn.test


sp.mn.train <- lapply(names, function(z){
  
  rowMeans(do.call(cbind, lapply(sp.res.train, function(x) x[[z]])), na.rm = T)
  
})

names(sp.mn.train) <- names

sapply(sp.mn.test, mean)



save(eval.results, sp.mn.train, sp.mn.test, sp.res.test, sp.res.train, 
     file = file.path(resFolder, "sp_results_blockCV.rdata"))


### Plots 
# pdf(file.path(resFolder, "eval_metrics_pairs_test.pdf"))
# plot(sp.res.test[[1]])
# mtext("test")
# dev.off()

# pdf(file.path(resFolder,"eval_metrics_pairs_train.pdf"))
# plot(sp.res.train[[1]])
# mtext("test")
# dev.off()


pdf(file.path(resFolder, "eval_metrics_auc_test_train_blockCV.pdf"))
plot(sp.mn.train$auc, sp.mn.test$auc, xlim = c(0,1), ylim = c(0,1), xlab = "AUC train", ylab = "AUC test")
abline(0,1)
dev.off()

