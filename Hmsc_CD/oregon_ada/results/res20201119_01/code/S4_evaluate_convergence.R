
## Christian D
## 12/11/2020

## Only local: 
# setwd("J:/UEA/Oregon")
# dir()


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/Oregon_winscp")
dir()


library(Hmsc)
#library(colorspace)
vioP <- require(vioplot)
 

# get models in models
fs <- list.files("models", "models_thin.*", full.names = TRUE)
# fs

beta <- list()

# i = 1;
# i = 2

for (i in seq_along(fs)) {
  
  filename = fs[i]
    
  load(filename)
  
  # str(models, max.level = 1)
  #str(models[[1]]$X)
  #str(models[[1]]$Y)
  # ma <- matrix(nrow = models[[1]]$ns * models[[1]]$nc, ncol = nm)
  # m <- models[[1]]
  
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

# str(beta, max.level = 2)
save(beta, file = "results/beta_list.rdata")
# load("results/beta_list.rdata")

pdf(file = "results/MCMC_convergence%03d.pdf")
par(mfrow=c(2,1), mar = c(10,2,1,1))
# takes a matrix, or individual vectors... or a list... 
#vioplot::vioplot(ma, main="psrf(beta)")

# remove NAs coming through from try errors first:
if(vioP) sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2) else sapply(beta[!is.na(beta)], boxplot, las = 2)

dev.off()

