##### Univariate models ####


## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()

setwd("~/Oregon_winscp")

# previous workspace loadded by default... so 
rm(list = ls())


load(file.path("data", "allData.Rdata"))

## Look at change in deviance with univariate and all-bar-one models

library(Hmsc)
# vioP <- require(vioplot)

preds <- colnames(X.train)
length(preds)

thin <- 10
samples <- 250
nChains <- 2
(samples * thin) + ceiling(0.5*thin*samples)

# for testing
# thin <- 1
# samples <- 2
# nChains <- 2
# 

# i = 1
# i = 2

## set up parallel
library(foreach)
library(doParallel)


nCores <- 21

cl <- makeCluster(nCores)
registerDoParallel(cl)


prRes <- foreach(i = seq_along(preds),
                 
                 .combine=list,
                 .multicombine=TRUE,
                 .packages = c("Hmsc", "coda")) %dopar% {
                   
                   ## univariate
                   m2 <- Hmsc::Hmsc(
                     Y = Y.train, 
                     XData = X.train, XFormula = as.formula(paste("~", preds[i])),
                     phyloTree = P,
                     distr = "probit"
                   )
                   
                   m2 <- try(Hmsc::sampleMcmc(m2, samples = samples, thin=thin,
                                              transient = ceiling(0.5*thin*samples),
                                              verbose = 0, # about 1/3 of samples
                                              nChains = nChains, nParallel = 1))
                   
                   if(!inherits(m2, "try-error")){
                     mpost2 = Hmsc::convertToCodaObject(m2, spNamesNumbers = c(TRUE,FALSE), 
                                                        covNamesNumbers = c(TRUE,FALSE))
                     psrf.beta2 = coda::gelman.diag(mpost2$Beta,multivariate=FALSE)$psrf[,1]
                   } else psrf.beta2 <- NA
                   
                   
                   ## model fit
                   preds2 = Hmsc::computePredictedValues(m2)
                   MF2 <- Hmsc::evaluateModelFit(hM=m2, predY=preds2)
                   
                   # cross validation - predictive  evaluation metrics
                   partition2 <- createPartition(m2, nfolds = 2)
                   predsCV = try(computePredictedValues(m2, partition=partition2, nParallel = 1))
                   if(inherits(predsCV, "try-error")) MFCV2 <- NA else {
                     MFCV2 <- evaluateModelFit(hM=m2, predY=predsCV)
                     }
                   
                   
                   res <- list(psrf.beta2, MF2, MFCV2)
                   names(res) <- c(paste0("beta_", preds[i], "_uni"), 
                                   paste0("MF_",preds[i], "_uni"),
                                   paste0("MFCV_",preds[i], "_uni"))
                   
                   #return
                   res
                   
                   }
                   

# 1 minute for 5 samples, thin1, 

# str(prRes, max.level = 3)
# lapply(prRes, names)


save(prRes, file = paste0("results/predSelection/prRes_univariate_thin", thin, "_samp", samples, ".rdata"))

# load(paste0("results/predSelection/prRes_univariate_thin", thin, "_samp", samples, ".rdata"))
# 
# all <- unlist(prRes, recursive = F)
# str(all, max.level = 1)
# names(all)
# MF <- all[grepl("MF", names(all))]
# beta <- all[grepl("beta", names(all))]
# 
# mfRes <- sapply(MF, function(x) sapply(x, mean))

# sort by AUC
# barplot(mfRes[,order(mfRes["AUC",])], beside = T, las = 2)
# sort(mfRes[2,], decreasing = TRUE)

# sort by RMSE
# plot(sort(mfRes[1,]))
# sort(mfRes[1,])


# par(mfrow=c(2,1), mar = c(10,2,1,1))
# if(vioP) vioplot::vioplot(beta[!is.na(beta)], las = 2) else boxplot(beta[!is.na(beta)], las = 2)
