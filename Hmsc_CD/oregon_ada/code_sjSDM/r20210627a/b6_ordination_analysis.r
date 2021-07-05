#### ecocopula analysis

## copied initially from oregon_ada/code/P1_ecocopula_analysis_v3.r

options(echo=TRUE) # if you want see commands in output file

# library(mvabund)
# library(ecoCopula)
library(raster)
library(Rtsne)
library(vegan)
# library(dplyr)

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

# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210627a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_test_results.rdata")) # # eval.results, sp.res.test, sp.res.train


## Mean AUC per species (and other eval metrics)
str(sp.res.test, max.level = 1)
head(sp.res.test$auc)

## Filter species by auc
auc.filt <- 0.65
# threshold for presence absence data
# tr <- 0.5

# how many species after AUC filter?
sum(sp.res.test$auc > auc.filt, na.rm = T)

# incidence 
incidence <- colSums(otu.pa.csv)/nrow(otu.pa.csv)


load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) # pred.mn, pred.sd, 

## local
# load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 

dim(pred.mn)

## filter for species performance
pred.in <- pred.mn[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
dim(pred.in)


## load raster templates
load(file.path(gis_out, "templateRaster.rdata")) ## r, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting

rList <- lapply(data.frame(pred.in), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})

# plot(tmp)
rStack <- stack(rList)
#names(rStack) <- spp.in$best.name
rStack

# f <- 3
f <- 5
# f <- 10
# f <- 50

rStack.agg <- raster::aggregate(rStack, f)
rStack.agg
# plot(rStack.agg, 1)

pred.mod <- values(rStack.agg)
NAs <- complete.cases(pred.mod)
sum(NAs)

## reduced data set
Xmat <- pred.mod[NAs, ]
r <- raster(rStack.agg)

## Full data set
# Xmat <- pred.in
# r <- raster(rStack)
# NAs <- indNA

# pa version
# Xmat <- (pred.mod[indNA2, ] >= tr)*1


dim(Xmat)
Xmat[1:10, 1:10]

## T SNE version
# library(Rtsne)
# ?Rtsne

perplexity <- 30

3 * perplexity < nrow(Xmat) - 1

system.time(
  tsne <- Rtsne(Xmat, dims = 2, theta = 0.5, partial_pca = F, num_threads = 0) # don't think I'm using openMP
)

# system.time(
#   tsne <- Rtsne(Xmat, dims = 2, theta = 0.5, partial_pca = T, num_threads = 0)
# )

# plot(tsne$Y, asp = 1, pch = ".")
# str(tsne, max.level =1)
# plot(tsne$Y, asp = 1)


# library(vegan)

system.time(
  pca <- rda(X = Xmat, scale = T)
)

str(pca, max.level = 1)
# biplot(pca, pch = ".")
# screeplot(pca)
# round(cumsum(100*pca$CA$eig/sum(pca$CA$eig)),2)[1:10]
# pca$CA$eig[1:15]

# plot(pca$Ybar[,1:2], asp = 1, pch = ".")

## put site scores into raster
makeR <- function(r, siteScores, NAs) {

  rSites <- raster(r)
  rSites[] <- NA
  rSites[NAs] <- siteScores
  rSites
  
}

rSites1 <- makeR(r, tsne$Y[,1], NAs)
rSites2 <- makeR(r, tsne$Y[,2], NAs)

# plot(stack(rSites1, rSites2))

pcaR1 <- makeR(r, scores(pca, 1, "sites"), NAs)
pcaR2 <- makeR(r, scores(pca, 2, "sites"), NAs)
# plot(stack(pcaR1, pcaR2))

# plot(stack(rSites1, rSites2,pcaR1, pcaR2))


ord.stck <- stack(rSites1, rSites2,pcaR1, pcaR2)
names(ord.stck) <- c("tsne1", "tsne2", "pca1", "pca2")


pdf(file.path(plotsFolder, "ord_plots.pdf"))
plot(ord.stck)
dev.off()


save(ord.stck, pca, tsne, r, f, NAs, file = file.path(resFolder, "ord_res.rdata"))

# 
# save(r, f, indNA2, file = file.path(resFolder, "rast_template_data.rdata"))



# 
# 
# # make a pa matrix as mvabund object
# pred.prob <- mvabund::mvabund(pred.mod[indNA2, ])
# pred.pa <- mvabund::mvabund((pred.mod[indNA2, ] >= tr)*1)
# 
# pred.prob[1:10, 1:10]
# pred.pa[1:10, 1:10]
# 
# dim(pred.pa)
# # hist(log(colSums(pred.pa)))
# 
# 
# #######
# 
# # do glm model
# mod.pa <- mvabund::manyglm(pred.pa~1, family = binomial(link="cloglog"))
# # mod.prob <- mvabund::manyany(pred.prob~1, family = binomial(link="cloglog"))
# 
# # do ordination
# mod.pa.ord <- ecoCopula::cord(mod.pa)
# 
# # mod.prob.ord <- ecoCopula::cord(mod.prob)
# 
# # head(mod.prob.ord$scores)
# # length(mod.prob.ord$scores[,"Factor1"])
# 
# 
# # ## make site scores into raster
# rSites.pa <- raster(rStack.agg)
# rSites.pa[] <- NA
# rSites.pa[indNA2] <- mod.pa.ord$scores[,"Factor1"]
# rSites.pa
# 
# # plot(rSites.pa)
# 
# # # make into df for ggplot
# # coords <- xyFromCell(rSites, seq_len(ncell(rSites)))
# # df1 <- as.data.frame(values(rSites))
# # df1 <- cbind(coords, df1)
# 
# # Species scores
# # sp_res <- data.frame(mod.ord$loadings, species = colnames(dataN$otu.pa))
# 
# r <- raster(rStack.agg)
# 
# 
# # save(mod.prob, mod.prob.ord, rSites.pa, file = file.path(resFolder, "ecocop_res_prob.rdata"))
# save(mod.pa, mod.pa.ord, rSites.pa, r, f, indNA2, file = file.path(resFolder, "ecocop_res.rdata"))
# 
# save(r, f, indNA2, file = file.path(resFolder, "rast_template_data.rdata"))
