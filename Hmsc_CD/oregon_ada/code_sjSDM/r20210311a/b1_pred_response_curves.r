## Data set up for models on ADA - 5 fold cross validation

## Run on ADA

## Creates response curves for each variable (partial plots)
## 1. Create full model
## 2. Predict along a gradiente (min to max) of each variable in turn, with rest held at mean 
## or level 1 (factor) - could also do over multiple factor levels... 
##  3. Plot


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

resFolder <-"code_sjSDM/r20210311a/results"

## load model data 
load(file.path(resFolder, "modelData.rdata"))
# if(!dir.exists(resFolder)) dir.create(resFolder, recursive = TRUE)

## update for testing
# device <- "cpu"
# iter= 10
# sampling = 100

## get best tune parameters
hidden <- list(c(50L,50L,10L), c(25L,25L,10L))
# use these for now -- best so far.
res <- read.csv(file.path(resFolder, "manual_tuning_sjsdm_10CV_M1S1_mean_AUC_pa_min_6_nSteps_1000.csv"))
res.best <- res[which.max(res$AUC.test_mean),,drop = T]
res.best


## Make full model

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
  select(all_of(vars))%>%
  mutate(across(where(is.numeric), scale))

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

    # spatial = NULL,
    
    learning_rate = 0.003, # 0.003 recommended for high species number 
    iter = iter, 
    family = family, 
    sampling = sampling, # 150L, 5000L
    device = device
)

## Non spatial
modFull_nsp <- sjSDM(
  
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
  
  # spatial = linear(data=scale.XY, ~0+UTM_E*UTM_N,
  #                  lambda=res.best$lambda.sp,
  #                  alpha=res.best$alpha.sp),
  # 
  spatial = NULL,
  
  learning_rate = 0.003, # 0.003 recommended for high species number 
  iter = iter, 
  family = family, 
  sampling = sampling, # 150L, 5000L
  device = device
)


## Save model
saveRDS(modFull_sp, file.path(resFolder,paste0("s-jSDM_fullMod_sp_", varsName, "_", spChoose, "_", abund, ".rds")))
saveRDS(modFull_nsp, file.path(resFolder,paste0("s-jSDM_fullMod_nsp_", varsName, "_", spChoose, "_", abund, ".rds")))

## Make newData for each var

varsN <- colnames(env.vars[,vars])[sapply(env.vars[,vars], is.numeric)]
varsF <- colnames(env.vars[,vars])[sapply(env.vars[,vars], is.factor)]

# v <- varsN[1]

n <- 250

predList_nsp <- vector(mode = "list", length = length(varsN))
predList_sp <- vector(mode = "list", length = length(varsN))

for(v in varsN) {
# for(v in varsN[1:3]) {
  
  ## Make gradient of each var
  vRange <- seq(min(env.vars[,v]), max(env.vars[,v]),length.out = n)
  
  # scale the vRange
  vRange.scale <- (vRange - mean(env.vars[,v]))/sd(env.vars[,v])
  
  # rest are at 0 -- scaled
  newdd <- data.frame(matrix(0, ncol = length(vars), nrow = n, dimnames = list(NULL, vars)))
  # head(newdd)
  
  # update focal var
  newdd[,v] <- vRange.scale
  
  # update factors
  level1 <- sapply(env.vars[,varsF, drop = F], levels)[1,]
  
  newdd[,varsF] <- level1
  
  # convert back to same factors
  for(f in varsF){
    
    newdd[, f] <- factor(newdd[,f], levels = levels(env.vars[,f]))
    
  }
  
  # head(newdd)
  # tail(newdd)
  # str(newdd)
  
  # do centroid of coordinates
  newSP <- matrix(c(0,0), ncol = 2, nrow = n, byrow = T, 
                  dimnames = list(NULL, colnames(scale.XY)))
  
  
  ## predict gives matrix of sites x species predictions
  pred.nsp = predict(modFull_nsp, newdata = newdd, SP = NULL)
  pred.sp = predict(modFull_sp, newdata = newdd, SP = newSP)
  
  # mean of 3 predictions
  # pred.dd = apply(abind::abind(lapply(1:3, function(i) 
  #   {predict(modFull, newdata=newdd, SP=newsp)}),along = -1L), 2:3, mean)
  # 
  
  ## Add vRange to pred matrix - for plotting
  pred.nsp <- cbind(vRange, pred.nsp)
  pred.sp <- cbind(vRange, pred.sp)
  
  colnames(pred.nsp) <- c(v, colnames(Y))
  colnames(pred.sp) <- c(v, colnames(Y))
  
  # head(pred.nsp)
  
  # attr(pred.dd, 'dimnames') = NULL

  predList_sp[[which(v == varsN)]] <- pred.sp
  predList_nsp[[which(v == varsN)]] <- pred.nsp
  
  names(predList_sp)[which(v == varsN)] <- v
  names(predList_nsp)[which(v == varsN)] <- v
  
}


save(modFull_nsp, modFull_sp, predList_sp, predList_nsp, file= file.path(resFolder, "responseData.rdata"))




### Plots 
# pdf(file.path(resFolder, "eval_metrics_pairs_test.pdf"))
# plot(sp.res.test[[1]])
# mtext("test")
# dev.off()
# 
# # pdf(file.path(resFolder,"eval_metrics_pairs_train.pdf"))
# # plot(sp.res.train[[1]])
# # mtext("test")
# # dev.off()
# 
# 
# pdf(file.path(resFolder, "eval_metrics_auc_test_train.pdf"))
# plot(sp.mn.train$auc, sp.mn.test$auc, xlim = c(0,1), ylim = c(0,1), xlab = "AUC train", ylab = "AUC test")
# abline(0,1)
# dev.off()
# 
# 
