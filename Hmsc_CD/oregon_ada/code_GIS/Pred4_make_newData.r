## Create new data for prediction across study area ####

## Local
wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"
setwd(wd)
# dir()
getwd()

library(raster)
library(sf)
library(dplyr)


gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"


## Load rasters
load("Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata") ## r, indNA

# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))
# get names and name groups
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck
names(allBrck)

plot(allBrck$insideHJA)
plot(allBrck$cut_msk)
plot(allBrck$cut_40msk)

## Make template raster - mask
# m.na <- values(is.na(allBrck))
# colSums(m.na)

# indNA <- complete.cases(values(dropLayer(allBrck, "cut_r")))
# r <- raster(allBrck)
# r[] <- indNA
# 
# plot(r, colNA = "black")
# save(r, indNA, file = "Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata")

## Add transformed raster here
allBrck$lg_DistStream <- log(allBrck$DistStream + 0.001)
allBrck$lg_DistRoad = log(allBrck$DistRoad + 0.001)
allBrck$lg_cover2m_max = log(allBrck$l_Cover_2m_max + 0.001)
allBrck$lg_cover2m_4m = log(allBrck$l_Cover_2m_4m + 0.001)
allBrck$lg_cover4m_16m = log(allBrck$l_Cover_4m_16m + 0.001)

plot(allBrck[[c("lg_DistStream", "DistStream")]])
plot(allBrck[[c("lg_DistRoad", "DistRoad")]])

plot(allBrck[[c("lg_cover2m_max", "l_Cover_2m_max", "lg_cover2m_4m", "l_Cover_2m_4m")]])

## Scale whole data set - apart from categorical predictors
allBrck.sc <- scale(dropLayer(allBrck, c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")))
# stores scale parameters in the @data slot
allBrck.sc # in memory
str(allBrck.sc)

## add back categorical - but bring into memory first
catRasters <- readAll(allBrck[[c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")]])
catRasters[[1]]

allBrck.sc <- addLayer(allBrck.sc, catRasters)
names(allBrck.sc)

## save scaled rasters
save(r, indNA, allBrck.sc, file = "Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")
# load("Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")

## data frame of coordinates
newXY <- coordinates(r)

newXY.sc <- scale(newXY)
str(newXY.sc)

colnames(newXY.sc) <- c("UTM_E", "UTM_N")
## filter for NAs in predictors

newXY.sc <- newXY.sc[indNA,]
dim(newXY.sc)

## Load sample points
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

xy.sites <- st_coordinates(xy.utm)
head(xy.sites)

# scale sample site coords with same parameters as complete data set
attr(newXY.sc, "scaled:center")
attr(newXY.sc, "scaled:scale")

xy.sites.sc <- cbind(
  X = (xy.sites[,"X"] - attr(newXY.sc, "scaled:center")["x"])/attr(newXY.sc, "scaled:scale")["x"],
  Y = (xy.sites[,"Y"] - attr(newXY.sc, "scaled:center")["y"])/attr(newXY.sc, "scaled:scale")["y"]
)

# Join sitenames and change names
xy.sites.sc <- data.frame(xy.sites.sc, SiteName = xy.utm$SiteName)
colnames(xy.sites.sc) <- c("UTM_E", "UTM_N", "SiteName")
head(xy.sites.sc)

## extract site env vars from scaled data set
allVars.sc <- data.frame(SiteName = xy.utm$SiteName, raster::extract(allBrck.sc, xy.utm))
head(allVars.sc)

## Get new data as data.frame:
newData.sc <- data.frame(values(dropLayer(allBrck.sc, "cut_r")))
head(newData.sc)
str(newData.sc)

length(indNA)

### remove NAs, for faster prediction and less storage - add back in for raster creation
newData.sc <- newData.sc[indNA, ] # only complete cases

## change categorical to predictor values to match model
newData.sc[, "insideHJA"] <- ifelse(newData.sc[, "insideHJA"] == 0, "no", "yes")
table(newData.sc[,"insideHJA"])

summary(newData.sc)

save(newData.sc, xy.sites.sc, newXY.sc, allVars.sc, file = "Hmsc_CD/oregon_ada/data/newData_scaled.rdata")
