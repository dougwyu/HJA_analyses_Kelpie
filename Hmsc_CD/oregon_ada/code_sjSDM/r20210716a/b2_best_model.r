## Data set up for models on ADA - 5 fold cross validation

## Run on ADA

## Sets up multiple Cross validation on ADA. 

## 0. Contents #####
## 1. Gets data from model tuning run (data saved for validation and testing beforehand)
## 2. Gets best tuning parameters for final model
## 3. Creates scaled data for full validation model and scales test set with these parameters
## 4. Evaluates test data with full validation model
## 4. Saves full model and test results

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
# [1] 0.1.3.9000
getwd() # always run sub from oregon_ada

library(dplyr)

resFolder <-"code_sjSDM/r202107016a/results"
abund <- "pa"

## load model data 
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
# if(!dir.exists(resFolder)) dir.create(resFolder, recursive = TRUE)

## Updated to new vars, also changes to elevation_m, canopy_height_m  to _f. 

### 2. Make validation full model and test on test data

set.seed(100)

# set variables from best tune
# formula.env = 'envDNN'
hidden <- list(c(50L,50L,10L), c(25L,25L,10L))

## get best tune run
# use these for now -- best so far.
# res <- read.csv(file.path(resFolder, "manual_tuning_sjsdm_5CV_M1S1_mean_AUC_pa_min_5_nSteps_1000.csv"))

res <- read.csv(file.path(resFolder,paste0("manual_tuning_sjsdm_", varsName, "_", k, "CV_", spChoose, 
                                           "_meanEVAL_", 
                                           abund, 
                                           "_min",
                                           minocc,
                                           "_nSteps",
                                           noSteps,
                                           ".csv")))


head(res)
res.best <- res[which.max(res$AUC.test_mean),,drop = T]
res.best


## Testing results table

# add in lambda.bio and alpha.bio from random samples of sample.bio, and add results columns
n <- 2
eval.results <- data.frame(type = rep(NA, n),
                           loglike_sjSDM = rep(NA, n),
                           loss=rep(NA, n),
                           AUC = rep(NA, n),
                           ll = rep(NA, n),
                           nagel= rep(NA, n),
                           plr = rep(NA, n),
                           tjur =rep(NA, n),
                           cor = rep(NA, n),
                           auc.lt5= rep(NA, n))

# clean up

# Choose pa or qp reponse data and family
if(abund == "pa") {
  Y <- otu.pa.csv
  Y.test <- otu.pa.csv.test
  family <- stats::binomial('probit') } else {
    if(abund ==  "qp") {
      Y <- otu.qp.csv
      Y.test <- otu.qp.csv.test
      family <- stats::poisson('log')
    } else stop("check abund")
  } 



### 3. Do full model with validation data and test on test data
# select X data
env.train <- env.vars[, vars]
env.test <- env.vars.test[, vars]
  
# select spatial data
XY.train <- env.vars[, c("UTM_E", "UTM_N")]
XY.test <- env.vars.test[, c("UTM_E", "UTM_N")]
  
# select Y data
s.otu.train <- as.matrix(Y)
s.otu.test <- as.matrix(Y.test)
  
 
factV <- colnames(env.train)[sapply(env.train, is.factor)]
  
# scale X and spatial data for full validation set  # .. env data - without factors
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
  
# do model    
model.train <- sjSDM(
    
    Y = s.otu.train,
    
    env = DNN(data=scale.env.train, formula = ~.,
              hidden=hidden[[res.best$hidden.ind]],
              lambda = res.best$lambda.env,
              alpha = res.best$alpha.env,
              activation = res.best$acti.sp,
              dropout=res.best$drop,
              bias=T),
    
    biotic = bioticStruct(lambda=res.best$lambda.bio,
                          alpha=res.best$alpha.bio, 
                          on_diag=F, inverse = FALSE),
    
    spatial = linear(data=XY.train.scale, ~0+UTM_E*UTM_N, 
                     lambda=res.best$lambda.sp, 
                     alpha=res.best$alpha.sp),
    
    learning_rate = res.best$lr, # 0.003 recommended for high species number 
    iter = iter, 
    family = family, 
    sampling = sampling, # 150L, 5000L
    device = device
)
  
## Save model
saveRDS(model.train,
          file.path(resFolder,paste0("s-jSDM_final_model_",varsName, "_",spChoose, "_",abund, ".rds")))
  
# Do testing and save results in data frame (2nd row is training)
eval.results[2,]$loglike_sjSDM <- logLik(model.train)
eval.results[2,]$loss <- model.train$history[length(model.train$history)]
  

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
    rsq$cor[m] = suppressWarnings(cor(p, y)) # just when speceis missing in test data
    
  }
  
  
  # Tjur <- tryCatch(expr = mapply(function(x,y) {
  #   tj <- base::diff(tapply(x, y, mean, na.rm = T)) # difference of average predicted values at 1 and 0
  #   if(!length(tj) > 0) tj <- NA
  #   return(tj)
  # }, asplit(pred.dd,2), asplit(otudd.pa,2)), error = function(err){ return(NA)})
  # 
  
  ## add to data frame
  eval.results[pred,]$type = c("test", "training")[pred]
  eval.results[pred,]$AUC = mean(auc, na.rm = TRUE)
  eval.results[pred,]$ll = mean(rsq$ll, na.rm = T)
  eval.results[pred,]$nagel = mean(rsq$nagel, na.rm = T)
  eval.results[pred,]$plr = mean(rsq$plr, na.rm = T)
  eval.results[pred,]$tjur = mean(rsq$tjur, na.rm = T)
  eval.results[pred,]$cor = mean(rsq$cor, na.rm = T)
  eval.results[pred,]$auc.lt5 = sum(auc < 0.5, na.rm = T)
  
  ## species results
  
  if (pred==1) {
    rsq$auc <- auc
    sp.res.test <- rsq
  }
  
  if (pred==2) {
    rsq$auc <- auc
    sp.res.train <- rsq
  }
  
  
} # end of evaluation loop
  
save(eval.results, sp.res.test, sp.res.train, file = file.path(resFolder, "sp_test_results.rdata"))

eval.results

# sp.res.test
# sp.res.train

colMeans(sp.res.test, na.rm = T)

# ## species means
# str(sp.res.test, max.level = 2)

# ### Plots 
pdf(file.path(resFolder, "eval_metrics_pairs_test.pdf"))
plot(sp.res.test)
mtext("test")
dev.off()
# 
# # pdf(file.path(resFolder,"eval_metrics_pairs_train.pdf"))
# # plot(sp.res.train[[1]])
# # mtext("test")
# # dev.off()
# 
# 
pdf(file.path(resFolder, "eval_metrics_auc_test_train.pdf"))
plot(sp.res.train$auc, sp.res.test$auc, xlim = c(0,1), ylim = c(0,1), xlab = "AUC train", ylab = "AUC test")
abline(0,1)
text(0.2,0.8, labels = sprintf("Mean test AUC = %0.2f", mean(sp.res.test$auc, na.rm = T)))
dev.off()
# 
