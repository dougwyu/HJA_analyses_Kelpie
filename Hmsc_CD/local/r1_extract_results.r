
## Extract and display results by model

library(Hmsc)
wd <- here::here()
setwd(file.path(wd, "Hmsc_CD"))

## load results
# get results folders 
resF <- list.files("oregon_ada/results", pattern = "res\\d*_\\d{2}$", include.dirs = TRUE, full.names = T)
resF <- c(resF, list.files("oregon_ada/results", pattern = "RRR_test2", include.dirs = TRUE, full.names = T))

## Loop through results folders, extract model definition, convergence results and evaluation

# rf <- resF[3]
modRes <- list()
meRes <- list()

for(rf in resF){
  
  # get model paths
  ms <- list.files(rf, pattern = "^models_.", full.names = TRUE, recursive = TRUE)
  
  # get model evaluation paths
  mf <- list.files(rf, pattern = "MF.*", full.names = TRUE, recursive = TRUE)
  
  # get highest thin and samples
  thin_vec <- sub(".*_thin_(\\d*)_.*", "\\1", ms)
  samples_vec <- sub(".*_samples_(\\d*)_.*", "\\1", ms)
  thin <- thin_vec[which.max(as.numeric(thin_vec))] # generally higher thin has highest samples.. for now.
  samples <- samples_vec[which.max(as.numeric(thin_vec))]
  rm(thin_vec, samples_vec)
  
  # load model 
  load(ms[grepl(paste0(".*_thin_", thin, "_samples_", samples, ".*\\.Rdata$"), ms)])

  # check for errors
  tryE <- !sapply(models, inherits, what = "try-error")
  
  # print models, predictors, no of Obs, random levels, etc
  #m <- models[[1]]; str(m)
  mod.specs <- lapply(models[tryE], function(m) {
    
    if(is.null(m$ranLevelsUsed)) rL <- NA else rL <- m$ranLevelsUsed
    if(is.null(m$XRRRFormula)) RRR_predictors = NA else RRR_predictors = as.character(m$XRRRFormula)[2]
    
    data.frame(
      predictors = as.character(m$XFormula)[2],
      obs = m$ny,
      nSp = m$ns,
      phylo = !is.null(m$phyloTree),
      rL = rL,
      RRR_predictors = RRR_predictors,
      ncRRR = m$ncRRR
      )
    
    })
  
  mod.df <- do.call(rbind, mod.specs)
  mod.df$name <- names(models)[tryE]
  rownames(mod.df) <- NULL  
  mod.df
  rm(mod.specs)
  
  ## Do convergence 
  load(list.files(rf, "beta.*\\.[rR]data", full.names = TRUE, recursive = TRUE))
  
  # Scale reduction factor - beta parameter, mean +/- sd
  all.psrf <- unlist(lapply(beta, function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y)))))
  mod.df$psrf <- all.psrf[match(paste(mod.df$name, thin, samples, sep = "_"), names(all.psrf))]
  rm(beta, all.psrf)
  
  ## do model evaluation
  load(mf[grepl(paste0(".*_thin_", thin, "_samples_", samples, ".*\\.Rdata$"), mf)]) # WAIC, MF, MFCV
  
  MFres <- t(do.call(data.frame,
                   lapply(MF[!is.na(MF)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))))
  MFCVres <- t(do.call(data.frame, 
                       lapply(MFCV[!is.na(MFCV)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))))
  
  colnames(MFres) <- paste0(colnames(MFres), "_expl")
  colnames(MFCVres) <- paste0(colnames(MFCVres), "_pred")
  
  mfs <- merge(MFres,MFCVres, by = "row.names")
  mod.df <- cbind(mod.df,  mfs[match(paste(mod.df$name, thin, samples, sep = "_"), mfs$Row.names),2:7])
  
  mod.df$thin <- as.numeric(thin)
  mod.df$samples <- as.numeric(samples)
  
  ## Collect all AUC data with species prevalence
  prevList <- lapply(models[tryE], function(m) {
    mev <- data.frame(prev = colSums(m$Y)/nrow(m$Y))
  })
  
  ind <- is.na(MF)
  
  ME_stats <- mapply(function(x,y,z,k) {
    
    mfdf <- do.call(data.frame, x)
    colnames(mfdf) <- paste0(colnames(mfdf), "_expl")
    
    mfcvdf <- do.call(data.frame, y)
    colnames(mfcvdf) <- paste0(colnames(mfcvdf), "_pred")
    
    
    cbind(rf = basename(rf),
          name = k, 
          mfdf,
          mfcvdf, 
          prev = z)}, MF[!ind], MFCV[!ind], prevList[!ind], names(prevList[!ind]), SIMPLIFY = FALSE)

  meRes[[which(rf == resF)]] <- do.call(rbind,ME_stats)
  
  rm(MF, MFCV, MFres, MFCVres, mfs, WAIC, thin, samples, ms, mf,modelnames)
  
  modRes[[which(rf == resF)]] <- mod.df
  
  }

all.df <- do.call(rbind, modRes)
all.me <- do.call(rbind, meRes)

rm(mod.df, modRes, resF, rf)

all.df <- all.df[order(all.df$AUC_pred, decreasing = TRUE),]
head(all.df)

head(all.me)

getwd()
save(all.df, all.me, file = "local/modelResults.rdata")

## Relationship of AUC with prevalence
library(tidyr)
library(ggplot2)

me.l <- pivot_longer(all.me, cols = colnames(all.me)[grepl("_expl|_pred", colnames(all.me))], 
                            names_to = c("metric", "type"), names_sep ="_")
# ends_with("_expl", "_pred")
head(me.l)
tail(me.l)
unique(me.l$type)

unique(me.l$rf)

res <- c("res20201202_01", "res20201204_01", "res20201209_01", "RRR_test2")

me.l %>%
  filter(!grepl("RRR|adj_evi|q", name) & rf %in% res & type == "pred") %>%
  ggplot(aes(x = prev, y = value))+
  geom_point()+
  facet_grid(rows = vars(metric), cols = vars(name), scales = "free_y")

getwd()
ggsave("local/prevPlot.png")  

ggplot(subset(me.l, rf %in% res & type == "pred"), aes(x = prev, y = value))+
  geom_point()+
  facet_wrap(~ metric + name)


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
