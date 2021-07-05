#### ecocopula analysis

## copied initially from oregon_ada/code/P1_ecocopula_analysis_v3.r

options(echo=TRUE) # if you want see commands in output file

library(mvabund)
library(ecoCopula)
library(raster)

# only on local
# getwd()
# wd <- here::here()
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

## on ADA
gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"

baseFolder <- "code_sjSDM/r20210610a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_results.rdata")) # sp.mn.test
rm(eval.results, sp.mn.train, sp.res.test, sp.res.train)

## Mean AUC per species (and other eval metrics)
str(sp.mn.test, max.level = 1)
head(sp.mn.test$auc)

## Filter species by auc
auc.filt <- 0.65
# threshold for presence absence data
tr <- 0.5

# how many species after AUC filter?
sum(sp.mn.test$auc > auc.filt)

# incidence 
incidence <- colSums(otu.pa.csv)/nrow(otu.pa.csv)


load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) # pred.mn, pred.sd, 

## local
# load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 

dim(pred.mn)

## filter for species performance
pred.in <- pred.mn[,sp.mn.test$auc > auc.filt]
dim(pred.in)

# reduce data set for testing
rFact <- 1000
smp <- seq(1, nrow(pred.in), rFact)

pred.mod <- pred.in[smp, ]
dim(pred.mod)

# make a pa matrix as mvabund object
pred.prob <- mvabund::mvabund(pred.mod)
pred.pa <- mvabund::mvabund((pred.mod >= tr)*1)

pred.prob[1:10, 1:10]
pred.pa[1:10, 1:10]

# hist(log(colSums(pred.pa)))


# do glm model
mod.pa <- mvabund::manyglm(pred.pa~1, family = binomial(link="cloglog"))
# mod.prob <- mvabund::manyany(pred.prob~1, family = binomial(link="cloglog"))

# do ordination
mod.pa.ord <- ecoCopula::cord(mod.pa)

# mod.prob.ord <- ecoCopula::cord(mod.prob)

# head(mod.prob.ord$scores)
# length(mod.prob.ord$scores[,"Factor1"])


## load raster templates
load(file.path(gis_out, "templateRaster.rdata")) ## r, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting


# ## make site scores into raster
rSites.pa <- raster(r.msk)
rSites.pa[] <- NA
rSites.pa[indNA][smp] <- mod.pa.ord$scores[,"Factor1"]
rSites.pa

# 
# rSites.prob <- r.msk
# #rSites.prob[indNA][smp] <- mod.prob.ord$scores[,"Factor1"]

plot(rSites.pa)

# # make into df for ggplot
# coords <- xyFromCell(rSites, seq_len(ncell(rSites)))
# df1 <- as.data.frame(values(rSites))
# df1 <- cbind(coords, df1)

# Species scores
# sp_res <- data.frame(mod.ord$loadings, species = colnames(dataN$otu.pa))

save(mod.prob, mod.prob.ord, file = file.path(resFolder, "ecocop_res_prob.rdata"))
# save(mod.pa, mod.pa.ord, file = file.path(resFolder, "ecocop_res.rdata"))


