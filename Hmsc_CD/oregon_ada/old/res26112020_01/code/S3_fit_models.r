#### Model set up #### 

## Christian D
## 12/11/2020

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

# load(file = "models/unfitted_models.rdata") #models, modelnames
#models

# samples_list = c(5,250,250,250,250,250)
# thin_list = c(1,1,10,100,1000,10000)

# trials
# samples_list = c(5,100,250)
# thin_list = c(1,5,10)

# samples_list = c(5,50)
# thin_list = c(1,5)


#number of iterations for each element above
(samples_list * thin_list) + ceiling(0.5*thin_list*samples_list)
# But will end up with samples_list  no of samples

# nChains = 2; Lst = 1; model = 1
nChains <- 2

# for(Lst in 1:2){ # for first trials only 

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  
  for (model in 1:nm) {
    print(paste0("model = ", modelnames[model]))
    m <- models[[model]]
    m <- try(sampleMcmc(m, samples = samples, thin=thin,
                   # adaptNf=rep(ceiling(0.4*samples*thin), m$nr),# ...  almost same as default
                   transient = ceiling(0.5*thin*samples), verbose = 0, # about 1/3 of samples
                   nChains = nChains, nParallel = nChains))
    models[[model]] <- m
  }
  
  filename <- paste(modFolder, "/models_thin_", as.character(thin), # modFolder is from S00_pipeline
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata", sep = "")
  
  save(models, modelnames, file=filename)
}

