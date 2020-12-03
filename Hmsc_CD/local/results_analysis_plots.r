
library(Hmsc)

here::i_am("local/results_analysis_plots.r")
wd <- here::here()
setwd(wd)

## load results
resFolder <- "res20201127_01"
fp <- file.path(wd, "oregon_ada/results", resFolder)
dir(fp, recursive = T)

# define thin, samples
thin <- 10
samples <- 250


## get model definitions
list.files(fp, paste0("^models_.*thin_", thin, "_samples_", samples, ".*\\.Rdata$"), full.names = T, recursive = T)
load(list.files(fp, paste0("^models_.*thin_", thin, "_samples_", samples, ".*\\.Rdata$"), full.names = T, recursive = T))

models

m <- models[[2]]
str(m, max.level=1)
m$XFormula
m$XRRRFormula
nrow(m$XData)
m$ranLevels
# no of RRR factors
m$ncRRR
m$phyloTree

names <- names(models)

tryE <- which(sapply(models, inherits, what = "try-error"))

sapply(models[-tryE], function(m) {

  cbind(
    predictors = as.character(m$XFormula),
    obs = nrow(m$XData),
    phylo = !is.null(m$phyloTree),
    rL = m$ranLevelsUsed,
    RRR_predictors = as.character(m$XRRRFormula),
    ncRRR = m$ncRRR
    
  )
  
  
})


## 1. Do model convergence checks ####
load(list.files(fp, "beta.*\\.rdata", full.names = T))
str(beta, max.level = 1)
sapply(beta, names)

# Scale reduction factor - beta parameter, mean +/- sd
lapply(beta, function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))

# 5 - 100 OK for non RRR, 
# 10 250 ok for RRR
# spatial RRR not working still.. CHECK

op <- par(mfrow=c(3,1), mar = c(10,3,1,1))
# remove NAs coming through from try errors first:
sapply(beta[!is.na(beta)], vioplot::vioplot, las = 2)
par(op)

## 2.Do model evalauation ###
list.files(fp, "MF.*10_samples_250.*", full.names = T, recursive = T)
load(list.files(fp, "MF.*10_samples_250.*", full.names = T, recursive = T))

print(sapply(MF[!is.na(MF)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y)))))
print(sapply(MFCV[!is.na(MFCV)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y)))))

## Species associations 
# residual associations (taking into account covariates)
m <- models[[2]]
name <- names(models[2])

omegaCor <- computeAssociations(m)
sLevel <- 0.95

toPlot <- ((omegaCor[[1]]$support > sLevel)
           + (omegaCor[[1]]$support < (1-sLevel)) > 0) * omegaCor[[1]]$mean

str(toPlot)
labels <-attr(toPlot, "dimnames")[[1]]

family <- sub(".*__[[:alpha:]]*_([[:alpha:]]*)_.*", "\\1", labels)
length(unique(family))
family <- as.factor(family)
levels(family)

# assign new labels to cor matrix
attr(toPlot, "dimnames") <- list(substr(family, 1,3), substr(family, 1,3))

pdf(file.path(fp,"spCorPlot.pdf"))
corrplot::corrplot(toPlot, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(1,1,1,1), tl.col = rainbow(length(unique(family)))[family], 
                   tl.cex = 0.4)
dev.off()

## Use latent variables to show species associations

etaPost = getPostEstimate(m, "Eta") # sites, circles, latent factor
lambdaPost = getPostEstimate(m, "Lambda") # species, triangles, latent factors

# how many latent factors are there? decided in model..
dim(etaPost$mean) # 88, 5 factors.

# as with ordination, first factors will have most information... 

m$XFormula
Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
       colVar = "lg_YrsDisturb", main = name, spNames = NA) 
# blue lowest values, red, highest values.  colorRampPalette(c("blue","white","red"))

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes
str(m$XData)

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = "oldGrowthIndex", main = name, spNames = NA) # numeric, low to high(red)



# 
Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,3),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(3, 4),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes

# modified biplot function to allow species colors... eg by taxonomy or trait.. 
biPlot_hmsc(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,2),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = "black", spCols = rainbow(length(unique(family)))[family], asp = 1) # levels' no, yes
