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

preds <- colnames(X.train)

thin <- 5
samples <- 50
nChains <- 2
(samples * thin) + ceiling(0.5*thin*samples)

# thin <- 1
# samples <- 2
# nChains <- 2
# 

# i = 1
# i = 2

## set up parallel
library(foreach)
library(doParallel)


nCores <- 2

cl <- makeCluster(nCores)
registerDoParallel(cl)


#prRes <- foreach(i = seq_along(preds), 
prRes <- foreach(i = 1, 
              .combine=list,
              .multicombine=TRUE,
              .packages = c("Hmsc", "coda")) %dopar% {
                
                # model formula
                XFormula <- as.formula(paste0("~ ", paste0(preds[-i], collapse = " + ")))
                
                ## all but
                m1 <- Hmsc::Hmsc(
                  Y = Y.train, 
                  XData = X.train, XFormula = XFormula,
                  phyloTree = P,
                  distr = "probit"
                )
                
                m1 <- try(Hmsc::sampleMcmc(m1, samples = samples, thin=thin,
                                     transient = ceiling(0.5*thin*samples), 
                                     verbose = 0, # about 1/3 of samples
                                     nChains = nChains, nParallel = 1))
                
                
                # convergence
                if(!inherits(m1, "try-error")){
                mpost1 = Hmsc::convertToCodaObject(m1, spNamesNumbers = c(TRUE,FALSE), 
                                                   covNamesNumbers = c(TRUE,FALSE))
                psrf.beta1 = coda::gelman.diag(mpost1$Beta,multivariate=FALSE)$psrf[,1]
                # str(psrf.beta1)
                } else psrf.beta1 <- NA
                
                
                ## model fit
                preds1 = Hmsc::computePredictedValues(m1)
                MF1 <- Hmsc::evaluateModelFit(hM=m1, predY=preds1)
                
                res <- list(psrf.beta1, MF1)
                names(beta) <- c(paste0("beta_", preds[i], "_all_but"), paste0("MF_",preds[i], "_all_but"))
                  
                # return
                res
                
              }


# about 15 mins on laptop, 2 cores, 2 predictors, 2 chains, 2 samples, 1 thin
# # 6,794,128K  AveRSS on ADA for same job. 22 mins (/2 for total time)  6.7Gb for 2 cores

save(prRes_allB, file = paste0("results/predSelection/prRes_allB_thin", thin, "_samp", samples, ".rdata"))