
## Christian D
## 12/11/2020


## On ADA server
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/oregon_ada/results/RRR_test3")
# dir()
rm(list = ls())

library(Hmsc)
#library(colorspace)
vioP <- require(vioplot) # make sure it's installed on server
 

# filenames of fitted models
dir("models")
fs <- list.files("models", "models_thin.*", full.names = TRUE)
# fs

beta <- list()

# i = 1;

for (i in seq_along(fs)) {
  
  filename = fs[i]
    
  load(filename)
  
  # str(models, max.level = 1)
  #str(models[[1]]$X)
  #str(models[[1]]$Y)
  # make a list of vectors - then it doesn't matter if they are different lengths, if using differnt predictors
  
  ma <- lapply(models, function(m){
    
    if(inherits(m, "try-error")) return(NA) else{
      mpost = convertToCodaObject(m, spNamesNumbers = c(TRUE,FALSE), covNamesNumbers = c(TRUE,FALSE))
      psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
      #str(psrf.beta); rownames(psrf.beta)
      
      # tmp = summary(psrf.beta)
      psrf.beta[,1]
    }
    
  })

  names(ma) <- paste0(modelnames, sub(".*_thin_([[:digit:]]*)_samples_([[:digit:]]*)_.*", "_\\1_\\2", filename))
  
  # str(ma, max.level= 1)
  beta[[i]] <- ma
  
  

}

save(beta, file = "results/beta.rdata")

# Scale reduction factor - beta parameter, mean +/- sd
lapply(beta, function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))

pdf(file = file.path("results", "MCMC_convergence%03d.pdf"))
par(mfrow=c(2,1), mar = c(10,3,1,1))
# remove NAs coming through from try errors first:
if(vioP) sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2) else sapply(beta[!is.na(beta)], boxplot, las = 2)
dev.off()

