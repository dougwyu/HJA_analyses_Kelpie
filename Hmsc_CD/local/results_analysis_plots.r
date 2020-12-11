
library(Hmsc)

here::i_am("local/results_analysis_plots.r")
wd <- here::here()
setwd(wd)

## load results
resFolder <- "RRR_test3"
fp <- file.path(wd, "oregon_ada/results", resFolder)
dir(fp, recursive = T)

# define thin, samples
thin <- 12
samples <- 100


## get model definitions
list.files(fp, paste0("^models_.*thin_", thin, "_samples_", samples, ".*\\.Rdata$"), full.names = T, recursive = T)
load(list.files(fp, paste0("^models_.*thin_", thin, "_samples_", samples, ".*\\.Rdata$"), full.names = T, recursive = T))

models

# m <- models[[2]]
# str(m, max.level=1)
# m$XFormula
# m$XRRRFormula
# nrow(m$XData)
# m$ranLevels
# # no of RRR factors
# m$ncRRR
# m$phyloTree

names <- names(models)

tryE <- !sapply(models, inherits, what = "try-error")

# print models, predictors, no of Obs, random levels, etc
sapply(models[tryE], function(m) {

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


## Relationship of AUC with prevalence

m <- models[[2]]
name <- names(models)[2]
m
name

# Get order of species prevalence - no of sites present at
sp.prev <- colSums(m$Y)/nrow(m$Y) # pa data

str(MFCV, max.level = 1)
mfcv <- MFCV[[paste(name, thin, samples, sep = "_")]]

plot(mfcv$AUC ~ sp.prev)
## lower limit of AUC little bit explained by prevalence.... ???

plot(mfcv$TjurR2 ~ sp.prev)
plot(mfcv$RMSE ~ sp.prev)

## Top predictor coefficients for highest AUC,... highest prevalence???

## Species associations 
# residual associations (taking into account covariates)

# sort species names by prevalence
which.max(sp.prev)
which.min(sp.prev)

sort(sp.prev)[c(1, length(sp.prev))]
prev.sort <- sort(sp.prev)

omegaCor <- computeAssociations(m)
str(omegaCor, max.level = 2)

# support level to show positive/negative interactions
sLevel <- 0.95
sLevel <- 0.99

toPlot <- ((omegaCor[[1]]$support > sLevel)
           + (omegaCor[[1]]$support < (1-sLevel)) > 0) * omegaCor[[1]]$mean

# str(toPlot)
# rearrange toPlot in order of species prevalence, rarest to most common
dimnames(toPlot)

head(dimnames(toPlot)[[1]])
prev.names[1]

toPlot_prev <- toPlot[match(names(prev.sort), dimnames(toPlot)[[1]]), match(names(prev.sort), dimnames(toPlot)[[2]])]
toPlot_prev[1:2, 1:2]
tmp <- (length(sp.prev)-1):length(sp.prev)
toPlot_prev[tmp, tmp]

labels <-attr(toPlot_prev, "dimnames")[[1]]
tx.order <- sub(".*__[[:alpha:]]*_([[:alpha:]]*)_.*", "\\1", labels)
tx.order <- as.factor(tx.order)
levels(tx.order)


# assign new labels to cor matrix
attr(toPlot_prev, "dimnames") <- list(substr(tx.order, 1,3), substr(tx.order, 1,3))

# No labels
#pdf(file.path(fp,"spCorPlot.pdf"))
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "n", mar = c(1,1,1,1),
                   tl.cex = 0.4)
# dev.off()

# tax order labels
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(1,1,1,1), tl.col = rainbow(length(unique(tx.order)))[tx.order], 
                   tl.cex = 0.4)

# Show quantiles of rarity/commonness with labels
quantile(prev.sort)
qLabs <- as.numeric(cut(prev.sort, quantile(prev.sort), include.lowest = T, labels = 1:4))
table(qLabs, useNA = "always")

print("\U25A0")
attr(toPlot_prev, "dimnames") <- list(rep("\U25A0", dim(toPlot_prev)[1]), rep("\U25A0", dim(toPlot_prev)[1]))

# show rarity quantiles
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(1,1,1,1), tl.col = c("grey", "black", "grey", "black")[qLabs], 
                   tl.cex = 0.45)


pdf(file.path(fp,"spCorPlot.pdf"))
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(2,1,2,1), tl.col = c("grey", "black", "grey", "black")[qLabs], 
                   tl.cex = 0.45)
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
             colVar = , main = name, spNames = NA) 


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

source(file.path(wd, "local","biplot_hmsc.r"))
# modified biplot function to allow species colors... eg by taxonomy or trait.. and removed par setting asp... see
# if species patterns can be seen better....

biPlot_hmsc(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,2),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = "black", spCols = rainbow(length(unique(tx.order)))[tx.order], asp = 1) # levels' no, yes

legend("top", ncol = 2, bty = "n", cex = 0.7, legend = unique(tx.order), pch = 17, col = rainbow(length(unique(tx.order))))
