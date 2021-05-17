
#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada

library(dplyr)
library(ggplot2)

resFolder <-"code_sjSDM/r20210513a/results"
abund <- "pa"

# Thus, a second model:
#   
#   Causal inference model to gain insight into species niches:  linear + light regularization (alpha = 1, lambda = 0.001 for env, space, and bio)
# 
# There is no tuning for this, so it should be fast to run them, but we probably need to do some spatial blocking (blockCV package in R) to remove the effect of space per se and isolate the effects of the environmental covariates. Even so, this will be quick, and we are interested in the covariates that appear most often. 

library(blockCV)
library(automap)
library(raster)
library(sf)


## load example model data
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
# vars, varsName, abund, minocc, otu.pa.csv, otu.qp.csv, env.vars
rm(otuenv, k, noSteps, device, iter, sampling)
head(env.vars) # can use SiteName here to filter points to this data set

gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis <- "data/gis"

## load template rasters
load(file.path(gis, "templateRaster.rdata")) ## r, indNA
r
plot(stack(r, r.aoi.pred), colNA = "black")

r.aoi.pred <- raster::mask(r.aoi.pred, r,maskvalue = 0)
plot(stack(r, r.aoi.pred), colNA = "black")


## Load sample points
## Load sample site points
load(file.path(gis, "sample_sites.rdata"))
xy.utm

## filter by points in model


# spatial blocking by specified range with random assignment
sb <- blockCV::spatialBlock(speciesData = xy.utm,
                   rasterLayer = r.aoi.pred,
                   theRange = 2500, # size of the blocks
                   k = 5,
                   selection = "systematic",
                   iteration = 100, # find evenly dispersed folds
                   numLimit = 0,
                   showBlocks = TRUE,
                   biomod2Format = FALSE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)

sb

str(sb, max.level = 1)

foldExplorer(sb,r.aoi.pred, xy.utm)
# str(sb$plots, max.level = 1)

sb$folds

## do spatial autocorrelation of predictors..... 
allBrck <- brick(file.path(gis_out, "r_utm/allStack_aoi.tif"))
allBrck
names(allBrck)
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- c(allNames, "lg_DistStream", "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m")

nCores <- parallel::detectCores()-1

autoRange <- spatialAutoRange(
  rasterLayer = allBrck[[c("ht30", "tpi250", "DistStream", "DistRoad", "ndvi_p50_r100")]],
  sampleNumber = 500,
  border = aoi.pred.sf,
  doParallel = TRUE,
  nCores = nCores,
  showPlots = TRUE)

plot(autoRange)
summary(autoRange)


## filter by only those on scale < 500m
names(allBrck)
rLayers <- names(allBrck)[!grepl("r500|1k", names(allBrck))]
# remove insideHJA and the cut binary
rLayers <- rLayers[!rLayers %in% c("cut_r", "insideHJA")]
rLayers

aRange <- spatialAutoRange(
  rasterLayer = allBrck[[rLayers]],
  sampleNumber = 5000,
  border = aoi.pred.sf,
  doParallel = TRUE,
  nCores = nCores,
  showPlots = TRUE)

save(aRange, file = "Hmsc_CD/oregon_ada/data/gis/aRange.rdata")

summary(aRange)
# plot(aRange)

str(aRange, max.level = 1)
plot(aRange$variograms[[1]])
# OJO corresponds to position 9 in summary table

plot(aRange$variograms[[9]])

## Check autocorrelation in presence/absence data

## eg one species

s.var <- autofitVariogram(otu.pa.csv[,7]~1, 
                 input_data = as(subset(xy.utm, SiteName %in% env.vars$SiteName), "Spatial"))

s.var
plot(s.var)
class(s.var)

dist <- st_distance(subset(xy.utm, SiteName %in% env.vars$SiteName))
table(cut(dist[lower.tri(dist)], breaks = s.var$exp_var$dist))
range(dist)

# check for all species
# semi variance is variance across binned points (wihtin distance band) , ie distance from mean of values within bin
# values to left of range are correlated, 

s.varL <- lapply(1:ncol(otu.pa.csv), function(i) {
  
  
  dat <- as(subset(xy.utm, SiteName %in% env.vars$SiteName), "Spatial")
  s.var <- autofitVariogram(otu.pa.csv[,i]~1, 
                            input_data = dat)
  s.var
})

str(s.varL[2])
s.varL[[1]]$var_model$range[2]

sapply(s.varL, function(x) x$var_model$range[2])

## get blocking CV ids

# spatial blocking by specified range with random assignment
sb <- blockCV::spatialBlock(speciesData = xy.utm,
                            rasterLayer = r.aoi.pred,
                            theRange = 3000, # size of the blocks
                            k = 5,
                            selection = "random",
                            iteration = 100, # find evenly dispersed folds
                            numLimit = 0,
                            showBlocks = TRUE,
                            biomod2Format = FALSE,
                            xOffset = 0, # shift the blocks horizontally
                            yOffset = 0)

sb
foldExplorer(sb,r.aoi.pred, xy.utm)

xy.utm$foldID <- sb$foldID

sb$plots+
  geom_sf(data = xy.utm, alpha = 0.5, aes(col = as.factor(foldID)))+
  scale_color_viridis_d(name = "fold id", option = "A")

sb.chk <- blockCV::spatialBlock(speciesData = xy.utm,
                            rasterLayer = r.aoi.pred,
                            theRange = 3000, # size of the blocks
                            k = 5,
                            selection = "checkerboard",
                            iteration = 100, # find evenly dispersed folds
                            numLimit = 0,
                            showBlocks = TRUE,
                            biomod2Format = FALSE,
                            xOffset = 0, # shift the blocks horizontally
                            yOffset = 0)

xy.utm$foldID <- sb.chk$foldID

sb.chk$plots+
  geom_sf(data = xy.utm, alpha = 0.5, aes(col = as.factor(foldID)))+
  scale_color_manual(name = "fold id", values = c("darkred", "darkblue"))




rangeExplorer(r.aoi.pred)

## with buffers

# buffering with presence-absence data
bf1 <- buffering(speciesData = xy.utm,
                 theRange = 2000,
                 species = NULL, # to count the number of presences and absences/backgrounds
                 spDataType = "PA", # presence-absence  data type
                 progress = TRUE)

bf1
str(bf1, max.level = 1)

bf1$folds[[1]]

foldExplorer(bf1,r.aoi.pred, xy.utm)

sp.sf <- cbind(subset(xy.utm, SiteName %in% env.vars$SiteName), otu.pa.csv)

bf2 <- buffering(speciesData = sp.sf,
                 theRange = 2000,
                 species = colnames(sp.sf)[2], # to count the number of presences and absences/backgrounds
                 spDataType = "PA", # presence-absence  data type
                 progress = TRUE)

bf3 <- buffering(speciesData = sp.sf,
                 theRange = 2000,
                 species = NULL, # to count the number of presences and absences/backgrounds
                 spDataType = "PA", # presence-absence  data type
                 progress = TRUE)
bf2
str(bf2, max.level = 1)
str(bf3, max.level = 1)

foldExplorer(bf2,r.aoi.pred, xy.utm)

# model trained on data outside of buffer distance around each point in turn
sapply(bf2$folds, lengths)
## eg some 2-12 poitns left out each model
hist(nrow(sp.sf)-sapply(bf2$folds, lengths)[1,])

sp.sf$fold1 <- 1:nrow(sp.sf) %in% bf2$folds[[1]][[1]]
sp.sf$fold2 <- 1:nrow(sp.sf) %in% bf2$folds[[2]][[1]]
sp.sf$fold3 <- 1:nrow(sp.sf) %in% bf2$folds[[3]][[1]]

ggplot(sp.sf, aes(col = fold1))+
  geom_sf()+
  geom_sf(data = sp.sf[1,], col = "black", shape = 17, size = 2)

ggplot(sp.sf, aes(col = fold2))+
  geom_sf()+
  geom_sf(data = sp.sf[2,], col = "black", shape = 17, size = 2)

ggplot(sp.sf, aes(col = fold3))+
    geom_sf()+
  geom_sf(data = sp.sf[3,], col = "black", shape = 17, size = 2)
