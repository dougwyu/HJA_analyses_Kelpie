## Create new data for prediction across study area ####

## Local
wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"
setwd(wd)
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")
# dir()
getwd()

library(raster)
library(sf)
library(dplyr)

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")

utm10N <- 32610

gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"

# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"

## create a reduced prediction area - convex hull around (all points + HJA) + buffer
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm; xy.all.utm

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

## convex hull 
aoi.pred.sf <- st_buffer(st_union(st_convex_hull(st_union(xy.utm)), st_convex_hull(hja.utm)), dist = 500)

st_write(aoi.pred.sf, "J:/UEA/Oregon/gis/s_utm/aoi_pred_sf_500.shp", delete_layer = T)

plot(aoi.pred.sf, add =F, col = NA)
plot(st_geometry(xy.utm),add =T)
plot(hja.utm, add = T, col = NA)


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

plot(allBrck$be30)
plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, col = NA)
plot(xy.utm, add = T, pch = 16, col = "black")

## Try difference buffer distances
# aoi.pred.lst <- lapply(c(500, 1000, 1500), function(i){
#   st_buffer(st_union(st_convex_hull(st_union(xy.utm)), st_convex_hull(hja.utm)), dist = i)
# })
# 
# ## add to map
# plot(allBrck$be30)
# plot(st_geometry(hja.utm), add = T, col = NA)
# plot(st_geometry(xy.utm), add = T, pch = 16, col = "black")
# sapply(seq_along(aoi.pred.lst), function(x) plot(aoi.pred.lst[[x]], add = T, col = NA, border = x))

### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))


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
aoi.pred.ext <- adjExt(st_bbox(aoi.pred.sf_edit), outF = "Extent", d = 30) ####
# make raster template from above extent
r.aoi.pred <- raster(aoi.pred.ext, crs = utm10N, res = 30)

origin(r.aoi.pred)

# mask it with aoi (points and HJA)
r.aoi.pred[] <- 1
r.aoi.pred <- raster::mask(r.aoi.pred, st_as_sf(aoi.pred.sf_edit))
plot(r.aoi.pred, colNA = "grey")

plot(hja.utm, add = T, col = NA)
plot(aoi.pred.sf, add =T, border = "blue")
plot(aoi.pred.sf_edit, add =T, col = NA)
plot(xy.utm, add = T, pch = 16, col = "black")

## crop raster brick
allBrck <- crop(allBrck, r.aoi.pred)

origin(allBrck)
allBrck
r.aoi.pred

names(allBrck)
brNames <- names(allBrck)
# names(allBrck) <- c(allNames, "lg_DistStream", "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m")
# save allBrck names
save(brNames, file = file.path(gis_out, "brNames.rdata"))

plot(allBrck$gt4_r30)

indNA <- complete.cases(values(dropLayer(allBrck, "cut_r")))
sum(!indNA)

## add NAs
r.msk <- r.aoi.pred
r.msk[!indNA] <- NA

plot(stack(r.aoi.pred, r.msk), colNA = "black")

## UPdate index of complete cases including those from new aoi edit
indNA <- !is.na(values(r.msk))

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

# Get new data prior from unscaled version
## extract site env vars from original data set
allVars <- data.frame(xy.all.utm[, c("SiteName", "trap", "period")], raster::extract(allBrck, xy.all.utm))
head(allVars)

## Get new data as data.frame:
allBrck
newData <- data.frame(values(dropLayer(allBrck, "cut_r")))
# remove NAs
newData <- newData[indNA, ] # only complete cases
## change categorical to predictor values to match model
newData[, "insideHJA"] <- ifelse(newData[, "insideHJA"] == 0, "no", "yes")
table(newData[,"insideHJA"], useNA = "always")

save(allVars, newData, indNA, file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_unscaled.rdata")


gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
load("J:/UEA/Oregon/gis/processed_gis_data/r_oversize/predRaster_sc.rdata") # r.msk, r.aoi.pred, indNA, allBrck.sc

## data frame of coordinates
newXY <- coordinates(r.msk)

newXY.sc <- scale(newXY)
str(newXY.sc)

## Load sample points
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm # and xy.all.utm - all sites
xy.all.utm

xy.sites <- st_coordinates(xy.all.utm) ## OJO some of the coordinates are repeated (different site names)
head(xy.sites)

# scale sample site coords with same parameters as complete data set
attr(newXY.sc, "scaled:center")
attr(newXY.sc, "scaled:scale")

xy.sites.sc <- cbind(
  X = (xy.sites[,"X"] - attr(newXY.sc, "scaled:center")["x"])/attr(newXY.sc, "scaled:scale")["x"],
  Y = (xy.sites[,"Y"] - attr(newXY.sc, "scaled:center")["y"])/attr(newXY.sc, "scaled:scale")["y"]
)

# Join sitenames and change names
xy.sites.sc <- data.frame(xy.sites.sc, xy.all.utm[, c("SiteName", "trap", "period")])
colnames(xy.sites.sc) <- c("UTM_E", "UTM_N", "SiteName", "trap", "period")
head(xy.sites.sc)

# fix colnames and NAs
colnames(newXY.sc) <- c("UTM_E", "UTM_N")
## filter for NAs in predictors

newXY.sc <- newXY.sc[indNA,]
dim(newXY.sc)
sum(!complete.cases(newXY.sc))

## extract site env vars from scaled data set
allVars.sc <- data.frame(xy.all.utm[, c("SiteName", "trap", "period")], raster::extract(allBrck.sc, xy.all.utm))
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



