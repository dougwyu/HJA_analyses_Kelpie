
## Local
wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"
setwd(wd)
# dir()
getwd()


# Thus, a second model:
#   
#   Causal inference model to gain insight into species niches:  linear + light regularization (alpha = 1, lambda = 0.001 for env, space, and bio)
# 
# There is no tuning for this, so it should be fast to run them, but we probably need to do some spatial blocking (blockCV package in R) to remove the effect of space per se and isolate the effects of the environmental covariates. Even so, this will be quick, and we are interested in the covariates that appear most often. 

library(blockCV)
library(raster)
library(sf)

gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"

## load template rasters
load("Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata") ## r, indNA
r
plot(stack(r, r.aoi.pred), colNA = "black")

r.aoi.pred <- raster::mask(r.aoi.pred, r,maskvalue = 0)
plot(stack(r, r.aoi.pred), colNA = "black")


## Load sample points
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

## filter by points in model


# spatial blocking by specified range with random assignment
sb <- blockCV::spatialBlock(speciesData = xy.utm,
                   rasterLayer = r.aoi.pred,
                   theRange = 2000, # size of the blocks
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

aRange <- spatialAutoRange(
  rasterLayer = allBrck[[c("ht30", "cut_r1k", "tpi250", "DistStream", "DistRoad", "ndvi_p50_r100")]],
  sampleNumber = 500,
  border = aoi.pred.sf,
  doParallel = TRUE,
  nCores = nCores,
  showPlots = TRUE)

aRange <- spatialAutoRange(
  rasterLayer = allBrck,
  sampleNumber = 5000,
  border = aoi.pred.sf,
  doParallel = TRUE,
  nCores = nCores,
  showPlots = TRUE)

save(aRange, file = "Hmsc_CD/oregon_ada/data/gis/aRange.rdata")

summary(aRange)
plot(aRange)
