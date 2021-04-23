
## Do canopy height locally (avoids uploading extra data to ADA)

getwd()
wd <- here::here()
setwd(wd)
dir()

library(sf)
library(raster)

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269

# gis <- "J:/UEA/Oregon/gis"
gis <- file.path(wd, "HJA_scripts/10_eo_data/raw_gis_data")


# ### Do canopy cover
hh <- raster(file.path(gis, "marie/latlong_highesthit.tif"))
hh
be <- raster(file.path(gis, "marie/latlong_bare_earth.tif"))
# crs: +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs 
plot(hh)
plot(be)

plot(stack(hh, be))

## canopy height
ht <- hh - be

# set values below 0 to 0
ht[ht < 0] <- 0
ht
plot(ht)

## create gap layer
ht_gt4m <- ht > 13.1234
ht_gt4m
plot(ht_gt4m)
writeRaster(ht_gt4m, filename = file.path(gis, "r_lcc/ht_gt4m.tif"), datatype = "INT1U")
writeRaster(ht, filename = file.path(gis, "r_lcc/ht_ft.tif"), datatype = "FLT4S")

## crs in raster package: 
# # crs: +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs

# warning after saving
# Warning message:
#   In showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj = prefer_proj) :
#   Discarded datum Unknown based on GRS80 ellipsoid in CRS definition

# CRS afer loading saved raster - same crs, same warning
ht_gt4m <- raster(file.path(gis, "r_lcc/ht_gt4m.tif"))
ht_gt4m
# +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs 

# separate script on ADA
# # do focal stats over gap layer
# ft <- 3.280839895 # ft to m
# cp.r500 <- (focalWeight(ht_gt4m, d = 500*ft, type = "circle")>0)*1
# cp.r250 <- (focalWeight(ht_gt4m, d = 250*ft, type = "circle")>0)*1
# 
# capCov.r250 <- focal(ht_gt4m, w = cp.r250, fun = mean, na.rm = T)
# capCov.r500 <- focal(ht_gt4m, w = cp.r500, fun = mean, na.rm = T)
# 
# writeRaster(capCov.r250, filename = file.path(gis, "r_lcc/capCov_250.tif"), datatype = "FLT4S")
# writeRaster(capCov.r500, filename = file.path(gis, "r_lcc/capCov_500.tif"), datatype = "FLT4S")

