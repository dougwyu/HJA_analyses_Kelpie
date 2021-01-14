## Create new predictors from cover, elevation

getwd()
wd <- here::here()
setwd(wd)
dir()
rm(wd)

library(sf)
library(raster)

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269

gis <- "J:/UEA/Oregon/gis"
## 
dir(gis)

# get points
# get data
source("Hmsc_CD/local/L1_read_data.r")
rm(Y.train.pa, Y.train.qp, P, otu.pa.csv, otu.qp.csv, X.train)

head(S.train)
xy.sf <- st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = nadutm10)
rm(S.train)

# transform to wgs utm to match rasters
xy.utm <- st_transform(xy.sf, crs = utm10N)
rm(xy.sf)

# write
# st_write(xy.utm, file.path(gis, "s_utm/m1s1_utm10.shp"), delete_layer = T)

# distance between points
dist <- st_distance(xy.utm)
diag(dist)<- NA
nnd <- apply(dist, 1, min, na.rm = T)
hist(nnd)
summary(nnd)

## bring in cut areas
cut <- st_read(file.path(gis, "marie/disturbance.shp"))
cut
cut[order(cut$YEAR, na.last = F),] # 0 shown in arcgis come in as NAs...
sum(is.na(cut$YEAR))
cut.utm <- st_transform(cut, crs = utm10N)
cut.utm
plot(cut.utm[, c("YEAR")])
plot(st_geometry(xy.utm), add = T, pch = 2, col = "black", cex = 1.5)
rm(cut)

# get elevation raster
# bareearth, projected to UTM 10 N, bilinear, 10m res, converted to m (/0.30480061) in arcgis
be <- raster(file.path(gis, "r_utm/bareEarth_m_10m.tif"))
be

# bareearth, projected to UTM 10 N, bilinear, 10m res, converted to m (/0.30480061) in arcgis
# be30 <- raster(file.path(gis, "r_utm/bareEarth_m_30m.tif"))
# be30

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
r250 <- (focalWeight(ht, d = 250, type = "circle")>0)*1
dim(r250)
r500 <- (focalWeight(ht, d = 500, type = "circle")>0)*1
r1k <- (focalWeight(ht, d = 1000, type = "circle")>0)*1

ht.r250 <- focal(ht, w = r250, fun = mean, na.rm = T)
ht.r500 <- focal(ht, w = r500, fun = mean, na.rm = T)
ht.r1k <- focal(ht, w = r1k, fun = mean, na.rm = T)

# same for cover
identical((focalWeight(cov2_4, d = 250, type = "circle")>0)*1,r250)

cov2_4.r250 <- focal(cov2_4, w = r250, fun = mean, na.rm = T)
cov2_4.r500 <- focal(cov2_4, w = r500, fun = mean, na.rm = T)
cov2_4.r1k <- focal(cov2_4, w = r1k, fun = mean, na.rm = T)

cov4_16.r250 <- focal(cov4_16, w = r250, fun = mean, na.rm = T)
cov4_16.r500 <- focal(cov4_16, w = r500, fun = mean, na.rm = T) 
cov4_16.r1k <- focal(cov4_16, w = r1k, fun = mean, na.rm = T)

covStack <- stack(ht, ht.r250, ht.r500, ht.r1k, 
                  cov2_4,cov2_4.r250, cov2_4.r500,cov2_4.r1k, 
                  cov4_16, cov4_16.r250, cov4_16.r500,cov4_16.r1k)

names(covStack) <- c("ht", "ht.r250", "ht.r500", "ht.r1k", 
                    "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k",
                    "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k")

terr <- terrain(be, opt = c("slope", "aspect", "TRI"), unit= "degrees")
terr

## Make eastness and northness
Nss <- cos(terr$aspect* pi / 180) # "Northness (aspect)" # convert to radians for cos/sin functions
Ess <- sin(terr$aspect* pi / 180)  # eastness

hist(terr$aspect)
hist(Ess)

plot(terr$aspect)
plot(Nss)

plot(Ess)

## TWI
# dynatopmodel::upslope.area
system.time(
  twi <- dynatopmodel::upslope.area(be, atb = T)
)
twi
# $atb is the twi
save(twi$atb, file = file.path(gis, "r_utm/twi.rdata"))

# user  system elapsed  # on laptop
# 1299.66    4.90 1424.00

terr <- stack(be, terr, Nss, Ess, twi$atb)
names(terr) <- c("be10","tri","slope","aspect","Nss","Ess","twi")


writeRaster(terr, bylayer = T, filename = file.path(gis, "r_utm/terr10.tif"), datatype = "FLT4S", suffix = "names", overwrite = TRUE)
writeRaster(covStack, bylayer = T, filename = file.path(gis, "r_utm/l.tif"), datatype = "FLT4S", suffix = "names")
save(terr, covStack, file = file.path(gis, "r_utm/elev_cov.rdata"))

# load(file.path(gis, "r_utm/elev_cov.rdata"))

# terr30 <- terrain(be30, opt = c("slope", "aspect", "TRI"), units= "degrees")
# terr30
# 
# ## Make eastness and northness
# Nss30 <- cos(terr30$aspect* pi / 180) # "Northness (aspect)"
# Ess30 <- sin(terr30$aspect* pi / 180)  # eastness
# 
# system.time(
#   twi30 <- dynatopmodel::upslope.area(be30, atb = T)
# )
# twi30
# 
# terr30 <- addLayer(terr30, Nss30, Ess30, twi30$atb)
# names(terr30) <- c( "tri","slope","aspect","Nss","Ess","twi")
# writeRaster(terr30, bylayer = T, filename = file.path(gis, "r_utm/terr30.tif"), datatype = "FLT4S", suffix = "names")
# 
# save(terr30, covStack, file = file.path(gis, "r_utm/elev_cov30.rdata"))

## do cut within 1 km
hist(cut.utm$YEAR)
str(cut)
# 0 values are NA, and outside study area
# convert to raster and do focal
extent(cut.utm)
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")
ext <- adjExt(extent(cut.utm))
r.temp <- raster(ext, crs = utm10N, res = 30)
r.temp
origin(r.temp)
cut.r <- rasterize(cut.utm, r.temp, field = "YEAR", fun = "max")
cut.r
plot(cut.r)

# make a binary cut layer
cut.msk <- cut.r
cut.msk
cut.msk[cut.msk > 1] <- 1
plot(cut.msk)

# convert NA to 0, as this is 0 disturbance, only important where no disturbance in whole focal window
cut.msk[is.na(cut.msk)] <- 0

r1k <- (focalWeight(cut.msk, d = 1000, type = "circle")>0)*1
cut.r1k <- focal(cut.msk, w = r1k, fun = sum, na.rm = T) # number of cells with cut from any year within 1k circle
# convert to proportion of 1k circle area
cut.r1k <- cut.r1k/sum(r1k) # number of cells
plot(cut.r1k)

# roughly equal... 
sum(r1k) * 30*30
pi * 1000^2

plot(cut.r1k)

cutStack <- stack(cut.r, cut.msk, cut.r1k)
names(cutStack) <- c("cut_r", "cut_msk", "cut_r1k")
save(cutStack, file = file.path(gis, "r_utm/cut_stack.rdata"))
writeRaster(cutStack, bylayer = T, filename = file.path(gis, "r_utm/disturb.tif"), suffix = "names", overwrite = TRUE)

# load(file.path(gis, "r_utm/cut_stack.rdata"))

## Extract values and add to data frame
cut.r1k.pt <- extract(cutStack$cut_r1k, xy.utm)

## average elevation
dem500 <- raster::extract(terr$be10, xy.utm, buffer = 500, fun=mean, na.rm = T)
# point elevation
dem.pt <- raster::extract(terr$be10, xy.utm)

## topographic index
mTopo <- dem.pt - dem500

## All terr
terr.pt <- raster::extract(terr, xy.utm)

# lidar cover
cov.pt <- raster::extract(covStack, xy.utm)

## make data frame
topo.df <- data.frame(uniqueID = xy.utm$uniqueID, terr.pt, cov.pt, be500 = dem500, mTopo, cut.r1k.pt)
head(topo.df)
str(topo.df)
summary(topo.df)

save(topo.df, file = "Hmsc_CD/oregon_ada/data/topo_data.rdata")

cat(paste(colnames(topo.df), collapse = '", "'))

pairs(topo.df[,c("be10", "slope", "Nss", "Ess", "twi", "ht", "ht.r1k", "cov2_4", "cov2_4.r1k", "cov4_16", "cov4_16.r1k", "mTopo", "cut.r1k.pt")])


c("be10", "slope", "Ess", "twi", "ht", "cov2_4", "cov4_16", "mTopo") # , "cut.r1k.pt"
c("be10", "slope", "Ess", "twi","ht.r1k", "cov2_4.r1k", "cov4_16.r1k", "mTopo") # , "cut.r1k.pt"





