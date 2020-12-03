### Evaluate Fit #####


## Only local: 
# setwd("J:/UEA/r_transfer/code_r/Oregon")
# dir()


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/Oregon") # relative to home director... not working directory... 
dir()



library(Hmsc)

#You may wish to loop over samples_list and thinning_list as done in Script S3 
nChains = 4
thin = 1
samples = 250 #250
print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))

filename_in = paste("models/models_thin_", as.character(thin),
                    "_samples_", as.character(samples),
                    "_chains_",as.character(nChains),
                    ".Rdata",sep = "")
filename_in

load(file = filename_in) #models, modelnames
# models[[1]]
nm = length(models)

MF = list()
MFCV = list()
WAIC = list()

# model = 1

for(model in 1:nm){
  
  print(paste0("model = ",as.character(model)))
  m = models[[model]]
  
  preds = computePredictedValues(m)
  MF[[model]] = evaluateModelFit(hM=m, predY=preds)
  
  partition = createPartition(m, nfolds = 2)
  preds = computePredictedValues(m,partition=partition, nParallel = nChains)
  
  MFCV[[model]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[model]] = computeWAIC(m)
}


filename_out = paste("models/MF_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                      ".Rdata",sep = "")

filename_out

save(MF,MFCV,WAIC,modelnames, file = filename_out)

sapply(MF, function(x) sapply(x, mean))
sapply(MFCV, function(x) sapply(x, mean))


