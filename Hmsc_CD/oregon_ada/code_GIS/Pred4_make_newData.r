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

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")

utm10N <- 32610

gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"


## create a reduced prediction area - convex hull around (all points + HJA) + buffer
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

## convex hull 
aoi.pred.sf <- st_buffer(st_union(st_convex_hull(st_union(xy.utm)), st_convex_hull(hja.utm)), dist = 1500)

st_write(aoi.pred.sf, "J:/UEA/Oregon/gis/s_utm/aoi_pred_sf.shp")

plot(st_geometry(xy.utm))
plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, col = NA)

# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))
# get names and name groups
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck
names(allBrck)

plot(allBrck$insideHJA)
# plot(allBrck$cut_msk)
# plot(allBrck$cut_40msk)

plot(hja.utm, add = T)
plot(aoi.pred.sf, add =T)
plot(xy.utm, add = T, pch = 16, col = "black")

## Add transformed raster here
allBrck$lg_DistStream <- log(allBrck$DistStream + 0.001)
allBrck$lg_DistRoad = log(allBrck$DistRoad + 0.001)
allBrck$lg_cover2m_max = log(allBrck$l_Cover_2m_max + 0.001)
allBrck$lg_cover2m_4m = log(allBrck$l_Cover_2m_4m + 0.001)
allBrck$lg_cover4m_16m = log(allBrck$l_Cover_4m_16m + 0.001)

plot(allBrck[[c("lg_DistStream", "DistStream")]])
plot(allBrck[[c("lg_DistRoad", "DistRoad")]])

plot(allBrck[[c("lg_cover2m_max", "l_Cover_2m_max", "lg_cover2m_4m", "l_Cover_2m_4m")]])

# get names for excel description
# write.csv(names(allBrck), "clipboard", row.names = F)

origin(allBrck)

##### CROP #####
## crop to reduced aoi.pred - get extent of the convex hull around sample points and HJA
aoi.pred.ext <- adjExt(st_bbox(aoi.pred.sf), outF = "Extent", d = 30) ####
# make raster template from above extent
r.aoi.pred <- raster(aoi.pred.ext, crs = utm10N, res = 30)

origin(r.aoi.pred)

# mask it with aoi (points and HJA)
r.aoi.pred[] <- 1
r.aoi.pred <- raster::mask(r.aoi.pred, st_as_sf(aoi.pred.sf))
plot(r.aoi.pred, colNA = "grey")

plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, border = "blue")
plot(xy.utm, add = T, pch = 16, col = "black")

## crop raster brick
allBrck <- crop(allBrck, r.aoi.pred)

origin(allBrck)
allBrck
r.aoi.pred

names(allBrck)
brNames <- names(allBrck)
# names(allBrck) <- c(allNames, "lg_DistStream", "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m")

plot(allBrck$gt4_r30)

indNA <- complete.cases(values(dropLayer(allBrck, "cut_r")))
sum(!indNA)

r.msk <- r.aoi.pred
r.msk[!indNA] <- NA

plot(stack(r.aoi.pred, r.msk), colNA = "black")

## Mask allBrick
allBrck <- mask(allBrck, r.msk)
plot(allBrck$gt4_r30)

# save raster as single tif
writeRaster(allBrck, filename = file.path(gis_out, "r_utm/allStack_aoi.tif"), overwrite = TRUE)

# plot(r, colNA = "black")
save(r.msk, indNA, r.aoi.pred, aoi.pred.sf, file = "Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata")
save(r.msk, indNA, r.aoi.pred, aoi.pred.sf, file = "J:/UEA/Oregon/gis/processed_gis_data/templateRaster.rdata")

## Scale whole data set - apart from categorical predictors
allBrck.sc <- scale(dropLayer(allBrck, c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")))
# stores scale parameters in the @data slot
allBrck.sc # in memory
# str(allBrck.sc)
inMemory(allBrck.sc[[1]])

sapply(1:nlayers(allBrck.sc), function(x) inMemory(allBrck.sc[[x]]))

## add back categorical - but bring into memory first
catRasters <- readAll(allBrck[[c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")]])
catRasters[[1]]

allBrck.sc <- addLayer(allBrck.sc, catRasters)
names(allBrck.sc)

## save scaled rasters
# save(r.msk, r.aoi.pred, indNA, allBrck.sc, file = "Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")
# load("Hmsc_CD/oregon_ada/data/gis/predRaster_sc.rdata")
save(r.msk, r.aoi.pred, indNA, allBrck.sc, file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/predRaster_sc.rdata")

# save scaled raster
writeRaster(allBrck.sc, filename = file.path(gis_out, "r_utm/AllStack_aoi_sc.tif"), overwrite = TRUE)

#### MAKE NEW DATA #####

## data frame of coordinates
newXY <- coordinates(r.msk)

newXY.sc <- scale(newXY)
str(newXY.sc)

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

# fix colnames and NAs
colnames(newXY.sc) <- c("UTM_E", "UTM_N")
## filter for NAs in predictors

newXY.sc <- newXY.sc[indNA,]
dim(newXY.sc)
sum(!complete.cases(newXY.sc))

## extract site env vars from scaled data set
allVars.sc <- data.frame(SiteName = xy.utm$SiteName, raster::extract(allBrck.sc, xy.utm))
head(allVars.sc)

## Get new data as data.frame:
allBrck.sc
newData.sc <- data.frame(values(dropLayer(allBrck.sc, "cut_r")))
head(newData.sc)
str(newData.sc)

length(indNA) == ncell(allBrck.sc)

### remove NAs, for faster prediction and less storage - add back in for raster creation
newData.sc <- newData.sc[indNA, ] # only complete cases

nrow(newData.sc) + sum(!indNA) == ncell(allBrck)

sum(!complete.cases(newData.sc))

## change categorical to predictor values to match model
newData.sc[, "insideHJA"] <- ifelse(newData.sc[, "insideHJA"] == 0, "no", "yes")
table(newData.sc[,"insideHJA"], useNA = "always")

summary(newData.sc)

#save(newData.sc, xy.sites.sc, newXY.sc, allVars.sc, file = "Hmsc_CD/oregon_ada/data/newData_scaled.rdata")
save(newData.sc, xy.sites.sc, newXY.sc, allVars.sc, 
     file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_scaled.rdata")



