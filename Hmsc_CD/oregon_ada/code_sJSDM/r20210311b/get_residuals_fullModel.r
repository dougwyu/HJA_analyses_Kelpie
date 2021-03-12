## Data set up for models on ADA - 5 fold cross validation

## Run on ADA

## Get residuals from full models - with spatial and without spatial components. 
## 1. Create full model
## 2. Predicted values - observered. 


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

resFolder <-"code_sjSDM/r20210311b/results"

## load model data 
load(file.path(resFolder, "modelData.rdata"))
# if(!dir.exists(resFolder)) dir.create(resFolder, recursive = TRUE)

## update for testing
device <- "cpu"
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


## predict gives matrix of sites x species predictions
pred.nsp = predict(modFull_nsp, newdata = NULL, SP = NULL)
pred.sp = predict(modFull_sp, newdata = NULL, SP = NULL)

# mean of 3 predictions
  # pred.dd = apply(abind::abind(lapply(1:3, function(i) 
  #   {predict(modFull, newdata=newdd, SP=newsp)}),along = -1L), 2:3, mean)
  # 

dim(as.matrix(Y))

resids.nsp <- scale(as.matrix(Y) - pred.nsp)
resids.sp <- scale(as.matrix(Y) - pred.sp)

## Mean and sd of residual across all species 
resids.nsp.mean <- rowMeans(resids.nsp)
resids.sp.mean <- rowMeans(resids.sp)

resids.nsp.sd <- apply(resids.nsp, 1, sd)
resids.sp.sd <- apply(resids.sp, 1, sd)



## get coordinates 

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269

library(sf)
pts.sf <- otuenv %>%
  select(c(UTM_N, UTM_E, SiteName, period, trap)) %>%
  dplyr::filter(period == "S1" & trap == "M1") %>%
  st_as_sf(coords = c("UTM_E", "UTM_N"), crs = nadutm10)

resids.sf <- cbind(pts.sf, mean.nsp = resids.nsp.mean, 
                   mean.sp = resids.sp.mean, 
                   sd.nsp = resids.nsp.sd,
                   sd.sp = resids.sp.sd)

resids.sf

save(pred.sp, pred.nsp, resids.sp, resids.nsp, resids.sf, file= file.path(resFolder, "residualData.rdata"))


## PLOTS

pred.nsp <- data.frame(pred.nsp)
colnames(pred.nsp) <- colnames(resids.nsp)

library(ggplot2)

pred_resid.nsp <- data.frame(resids.nsp) %>%
  tibble::add_column(SiteName = env.vars$SiteName) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_to = "residual")%>%
  left_join(
    pred.nsp %>%
      tibble::add_column(SiteName = env.vars$SiteName) %>%
      tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_to = "predicted")
  )


ggplot(pred_resid.nsp, aes(x = predicted, y = residual))+
  geom_point(size = 0.1)+
  facet_wrap(~OTU)


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
