## Create new predictors from cover, elevation

getwd()
wd <- here::here()
wd
setwd(wd)
dir()
rm(wd)

library(sf)
library(raster)

utm10N <- 32610
gis <- "J:/UEA/Oregon/gis"
## 
dir(gis)

# get points
# get data
source("Hmsc_CD/local/L1_read_data.r")
rm(Y.train.pa, Y.train.qp, P, otu.pa.csv, otu.qp.csv, X.train)

head(S.train)
xy.sf <- st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = utm10N)
rm(S.train)

# write
st_write(xy.sf, file.path(gis, "s_utm/m1s1_utm10.shp"))

# distance between points
dist <- st_distance(xy.sf)
diag(dist)<- NA
nnd <- apply(dist, 1, min, na.rm = T)
hist(nnd)
summary(nnd)

## bring in cut areas

cut <- st_read(file.path(gis, "marie/disturbance.shp"))
cut.utm <- st_transform(cut, crs = utm10N)

# get elevation raster
# bareearth, projected to UTM 10 N, bilinear, 10m res, converted to m (/0.30480061) in arcgis
be <- raster(file.path(gis, "r_utm/bareEarth_m_10m.tif"))
be

# bareearth, projected to UTM 10 N, bilinear, 10m res, converted to m (/0.30480061) in arcgis
be30 <- raster(file.path(gis, "r_utm/bareEarth_m_30m.tif"))
be30

# get cover data as brick
list.files(file.path(gis,"r_utm"), "lidar")
# "lidar_metric_mosaic_rumple.tif" 
cov2_4 <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_Cover_2m_4m.tif"))
cov4_16 <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_Cover_4m_16m.tif"))
ht <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_p95.tif"))

cov2_4
cov4_16
ht

hist(ht)
hist(ht[ht<100])

# some extreme values...  6. remove these
sum(values(ht)>100, na.rm = T)

ht[ht>100] <- NA
hist(ht)
hist(cov2_4)

hist(be)

## do focal stats, at 250m, 500m, 1000 / 2
# make focal weigth matrix for radius of 125 250, 500
r125 <- (focalWeight(ht, d = 125, type = "circle")>0)*1
r125
r250 <- (focalWeight(ht, d = 250, type = "circle")>0)*1
r500 <- (focalWeight(ht, d = 500, type = "circle")>0)*1

ht.r125 <- focal(ht, w = r125, fun = mean, na.rm = T)
ht.r250 <- focal(ht, w = r250, fun = mean, na.rm = T)
ht.r500 <- focal(ht, w = r500, fun = mean, na.rm = T)

# same for cover
identical((focalWeight(cov, d = 125, type = "circle")>0)*1,r125)

cov2_4.r125 <- focal(cov2_4, w = r125, fun = mean, na.rm = T)
cov2_4.r250 <- focal(cov2_4, w = r250, fun = mean, na.rm = T)
cov2_4.r500 <- focal(cov2_4, w = r500, fun = mean, na.rm = T)

cov4_16.r125 <- focal(cov4_16, w = r125, fun = mean, na.rm = T)
cov4_16.r250 <- focal(cov4_16, w = r250, fun = mean, na.rm = T)
cov4_16.r500 <- focal(cov4_16, w = r500, fun = mean, na.rm = T) 

covStack <- stack(ht.r125, ht.r250, ht.r500, 
                  cov2_4.r125, cov2_4.r250, cov2_4.r500, 
                  cov4_16.r125, cov4_16.r250, cov4_16.r500)

names(covStack) <- c("ht.r125", "ht.r250", "ht.r500", 
                    "cov2_4.r125", "cov2_4.r250", "cov2_4.r500",
                    "cov4_16.r125", "cov4_16.r250", "cov4_16.r500")

terr <- terrain(be, opt = c("slope", "aspect", "TRI"), units= "degrees")
terr

## Make eastness and northness
Nss <- cos(terr$aspect* pi / 180) # "Northness (aspect)"
Ess <- sin(terr$aspect* pi / 180)  # eastness


## topographic index

## average elevation

## TWI
# dynatopmodel::upslope.area
system.time(
  twi <- dynatopmodel::upslope.area(be, atb = T)
)
twi
# $atb is the twi

# user  system elapsed  # on laptop
# 1299.66    4.90 1424.00


terr <- addLayer(terr, Nss, Ess, twi$atb)
names(terr) <- c( "tri","slope","aspect","Nss","Ess","twi")

writeRaster(terr, bylayer = T, filename = file.path(gis, "r_utm/terr.tif"), datatype = "FLT4S", suffix = "names")
writeRaster(covStack, bylayer = T, filename = file.path(gis, "r_utm/l.tif"), datatype = "FLT4S", suffix = "names")

save(terr, covStack, file = file.path(gis, "r_utm/elev_cov.rdata"))

terr30 <- terrain(be30, opt = c("slope", "aspect", "TRI"), units= "degrees")
terr30

## Make eastness and northness
Nss30 <- cos(terr30$aspect* pi / 180) # "Northness (aspect)"
Ess30 <- sin(terr30$aspect* pi / 180)  # eastness

system.time(
  twi30 <- dynatopmodel::upslope.area(be30, atb = T)
)
twi30

terr30 <- addLayer(terr30, Nss30, Ess30, twi30$atb)
names(terr30) <- c( "tri","slope","aspect","Nss","Ess","twi")
writeRaster(terr30, bylayer = T, filename = file.path(gis, "r_utm/terr30.tif"), datatype = "FLT4S", suffix = "names")

save(terr30, covStack, file = file.path(gis, "r_utm/elev_cov30.rdata"))
