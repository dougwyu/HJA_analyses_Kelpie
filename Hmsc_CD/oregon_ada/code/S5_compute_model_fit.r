### Evaluate Fit #####


## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("~/Oregon_winscp")
# dir()
# rm(list = ls())


library(Hmsc)

# get models in models
# fs <- list.files("models", "^models_thin.*", full.names = TRUE)
#fs
#i = 1; i = 3

for (i in seq_along(fs)) {
  
  filename = fs[i]
  print(filename)
  nChains <- as.numeric(sub(".*_chains_([[:digit:]]{1}).*", "\\1", filename))
  #nChains = 2
  
  load(filename)
  
  MF = list()
  MFCV = list()
  WAIC = list()
  
  for(model in seq_along(models)){
    
    m <- models[[model]]
    
    if(inherits(m, "try-error")) {
      
      MF[[model]] <- NA
      MFCV[[model]] <- NA
      WAIC[[model]] <- NA
      
    } else {
      
      # explanatory - training evaluation metrics
      preds = computePredictedValues(m)
      MF[[model]] <- evaluateModelFit(hM=m, predY=preds)
      
      # cross validation - predictive  evaluation metrics
      partition <- createPartition(m, nfolds = nfolds)
      
      preds = try(computePredictedValues(m,partition=partition, nParallel = nChains))
      if(inherits(preds, "try-error")) MFCV[[model]] <- NA else {
        MFCV[[model]] <- evaluateModelFit(hM=m, predY=preds)
      }
      
      # WAIC..   Error in X %*% Beta : non-conformable arguments -- is there a problem in the X matrix.. just use
      # variables in the RRR??? 
      waic <- try(computeWAIC(m))
      if(inherits(waic, "try-error")) WAIC[[model]] <- NA else {
        WAIC[[model]] <- waic
      }
    }
    
  }
  
  m.names <- paste0(modelnames, sub(".*_thin_([[:digit:]]*)_samples_([[:digit:]]*)_.*", "_\\1_\\2", filename))
  names(WAIC) <- names(MFCV) <- names(MF) <- m.names
  
  
  
  filename_out <- file.path(dirname(filename), paste0("MF_", "CV", nfolds, "_", sub("^models_", "", basename(filename))))
  save(MF,MFCV,WAIC,modelnames, partition, file = filename_out)
  
  print(sapply(MF[!is.na(MF)], function(x) sapply(x, mean)))
  print(sapply(MFCV[!is.na(MFCV)], function(x) sapply(x, mean)))
  
  
}

ls()
