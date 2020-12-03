##### Univariate models ####


## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
dir()

setwd("~/Oregon_winscp")

load(file.path("data", "allData.Rdata"))

## Look at change in deviance with univariate and all-bar-one models

library(Hmsc)

preds <- colnames(X.train)

thin <- 5
samples <- 50
nChains <- 2
(samples * thin) + ceiling(0.5*thin*samples)

# thin <- 1
# samples <- 5
# nChains <- 2
# 

# i = 1
# i = 2

betaRes <- vector(mode = "list", length = length(preds))
MFRes <- vector(mode = "list", length = length(preds))

for(i in 1){
#for(i in seq_along(preds)){
  
  print(paste("no", i, preds[i]))
  
  # model formula
  XFormula <- as.formula(paste0("~ ", paste0(preds[-i], collapse = " + ")))
  
  ## all but
  m1 <- Hmsc(
    Y = Y.train, 
    XData = X.train, XFormula = XFormula,
    phyloTree = P,
    distr = "probit"
  )
  
  ## univariate
  m2 <- Hmsc(
    Y = Y.train, 
    XData = X.train, XFormula = as.formula(paste("~", preds[i])),
    phyloTree = P,
    distr = "probit"
  )
  
  
  m1 <- try(sampleMcmc(m1, samples = samples, thin=thin,
                      transient = ceiling(0.5*thin*samples), 
                      verbose = 0, # about 1/3 of samples
                      nChains = nChains, nParallel = nChains))
  
  
  m2 <- try(sampleMcmc(m2, samples = samples, thin=thin,
                       transient = ceiling(0.5*thin*samples), 
                       verbose = 0, # about 1/3 of samples
                       nChains = nChains, nParallel = nChains))
  
  
  # convergence
  mpost1 = convertToCodaObject(m1, spNamesNumbers = c(TRUE,FALSE), covNamesNumbers = c(TRUE,FALSE))
  psrf.beta1 = gelman.diag(mpost1$Beta,multivariate=FALSE)$psrf[,1]
  # str(psrf.beta1)
  
  if(!inherits(m2, "try-error")){
    mpost2 = convertToCodaObject(m2, spNamesNumbers = c(TRUE,FALSE), covNamesNumbers = c(TRUE,FALSE))
    psrf.beta2 = gelman.diag(mpost2$Beta,multivariate=FALSE)$psrf[,1]
  } else psrf.beta2 <- NA
  
  
  beta <- list(psrf.beta1, psrf.beta2)
  names(beta) <- c(paste0("beta_", preds[i], "_all_but"), 
                   paste0("beta_", preds[i], "_uni"))
  betaRes[[i]] <- beta
  
  ## model fit
  preds1 = computePredictedValues(m1)
  MF1 <- evaluateModelFit(hM=m1, predY=preds1)
  
  preds2 = computePredictedValues(m2)
  MF2 <- evaluateModelFit(hM=m2, predY=preds2)
  
  # cross validation - predictive  evaluation metrics
  # partition1 <- createPartition(m1, nfolds = 2)
  # preds1 = try(computePredictedValues(m1, partition=partition1, nParallel = nChains))
  # if(inherits(preds1, "try-error")) MFCV1 <- NA else {
  #   MFCV1 <- evaluateModelFit(hM=m1, predY=preds1)
  
   mf <- list(MF1, MF2)
   names(mf) <- c(paste0("MF_",preds[i], "_all_but"),
                   paste0("MF_",preds[i], "_uni"))
   MFRes[[i]] <- mf


}


save(betaRes, MFRes, file = "results/prRes_loop.rdata")

# str(betaRes, max.level = 2)
# x <- betaRes[[1]]

# lapply(betaRes, vioplot::vioplot, las = 1)

# x <- MFRes[[1]]
# lapply(MFRes, function(x) sapply(unlist(x, recursive = F), mean))


# sapply(beta[!is.na(beta)], boxplot, las = 2)

# # results
# ev1 <- sapply(MF1, mean)
# ev1 <- data.frame(pred = preds[i], type = "allB", RMSE = ev1[])
