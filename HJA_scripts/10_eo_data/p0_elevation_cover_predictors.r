## Create new predictors from cover, elevation

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
 
dir(gis)
# raster files are here:
dir(file.path(gis, "r_utm"))

# get points
samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",minimaprundate,"_",samtoolsfilter,"_", 
                                       samtoolsqual, "_kelpie", kelpierundate,"_", primer,"_vsearch97")

datFile <- paste0("sample_by_species_table_", samtoolsfilter, "_minimap2_", minimaprundate,"_kelpie",
                  kelpierundate,"_uncorr.csv")

otuenv <- read.csv(file.path(gitHub, outputidxstatstabulatefolder, datFile))

otuenv[1:6,1:10]
coords <- unique(otuenv[,c("SiteName", "UTM_E", "UTM_N")])
xy.sf <- st_as_sf(coords, coords = c("UTM_E", "UTM_N"), crs = nadutm10)

rm(gitHub, otuenv, outputidxstatstabulatefolder, datFile, primer, 
   kelpierundate, minimaprundate, samtoolsfilter, samtoolsqual)

# transform to wgs utm to match rasters
xy.utm <- st_transform(xy.sf, crs = utm10N)
rm(xy.sf)

length(unique(xy.utm$SiteName))

# write
# st_write(xy.utm, file.path(gis, "s_utm/m1s1_utm10.shp"), delete_layer = T)

# distance between points
dist <- st_distance(xy.utm)
diag(dist)<- NA
nnd <- apply(dist, 1, min, na.rm = T)
hist(nnd)
summary(nnd)

## bring in cut areas
cut <- st_read(file.path(gis, "s_nad_utm/disturbance.shp"))
cut
cut[order(cut$YEAR, na.last = F),] # 0 shown in arcgis come in as NAs...
sum(is.na(cut$YEAR))
cut.utm <- st_transform(cut, crs = utm10N)
cut.utm
plot(cut.utm[, c("YEAR")])
plot(st_geometry(xy.utm), add = T, pch = 2, col = "black", cex = 1.5)
rm(cut)

# get elevation raster
# bareearth, projected to UTM 10N wgs84, bilinear, at 10m res, converted to m (/0.30480061) in arcgis
be10 <- raster(file.path(gis, "r_utm/bareEarth_m_10m.tif"))
be10

# latlong_bare_earth.tif (2.8GB) projected to UTM 10 N, bilinear, 10m res, converted to m (/0.30480061) in arcgis
# saved as "r_utm/bareEarth_m_30m.tif"

# get cover data as brick
list.files(file.path(gis,"r_utm"), "lidar")
# "lidar_metric_mosaic_rumple.tif" 
cov2_4 <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_Cover_2m_4m.tif"))
cov4_16 <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_Cover_4m_16m.tif"))
ht <- raster(file.path(gis, "r_utm", "lidar_metric_mosaic_p95.tif"))

cov2_4
cov4_16
ht

plot(ht)

hist(ht)
hist(ht[ht<100])

# some extreme values...  6. remove these
sum(values(ht)>100, na.rm = T)

ht[ht>100] <- NA
hist(ht)
plot(ht)

hist(cov2_4)

hist(be10)


## do focal stats, at 250m, 500m, 1000 / 2
# make focal weigth matrix for radius of 250, 500, 1000
r250 <- (focalWeight(ht, d = 250, type = "circle")>0)*1
dim(r250)
r500 <- (focalWeight(ht, d = 500, type = "circle")>0)*1
r1k <- (focalWeight(ht, d = 1000, type = "circle")>0)*1

ht.r250 <- focal(ht, w = r250, fun = mean, na.rm = T)
ht.r500 <- focal(ht, w = r500, fun = mean, na.rm = T)
ht.r1k <- focal(ht, w = r1k, fun = mean, na.rm = T)

d250 <- (focalWeight(be10, d = 250, type = "circle")>0)*1
dim(d250)
prod(dim(d250))/2 # get central cell
d500 <- (focalWeight(be10, d = 500, type = "circle")>0)*1
d1k <- (focalWeight(be10, d = 1000, type = "circle")>0)*1

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

terr <- terrain(be10, opt = c("slope", "aspect", "TRI", "TPI"), unit= "degrees")
terr

# f <- matrix(1, nrow=3, ncol=3)
# TPI <- focal(be10, w=f, fun=function(x, ...) x[5] - mean(x[-5]), pad=TRUE, padValue=NA)
# 
# plot(terr$tpi)
# plot(TPI)

## Do TPI at three different scales
# d9 <- matrix(1,3,3)
# d250 <- (focalWeight(be10, d = 250, type = "circle")>0)*1
# dim(d250)
# d500 <- (focalWeight(be10, d = 500, type = "circle")>0)*1
# d1k <- (focalWeight(be10, d = 1000, type = "circle")>0)*1
# 
# ## Do TPI at three different scales
# c250 <- ceiling(prod(dim(d250))/2)
# c500 <- ceiling(prod(dim(d500))/2) # get central cell
# c1k <- ceiling(prod(dim(d1k))/2) # get central cell
# 
# 
# TPI250 <- focal(be10, w=d250, fun=function(x, ...) x[c250] - mean(x[-c250]), pad=TRUE, padValue=NA)
# writeRaster(TPI250, filename = "data/tpi250.tif", datatype = "FLT4S")
# 
# TPI500 <- focal(be10, w=d500, fun=function(x, ...) x[c500] - mean(x[-c500]), pad=TRUE, padValue=NA)
# writeRaster(TPI500, filename = "data/tpi500.tif", datatype = "FLT4S")
# 
# TPI1k <- focal(be10, w=d1k, fun=function(x, ...) x[c1k] - mean(x[-c1k]), pad=TRUE, padValue=NA)
# writeRaster(TPI1k, filename = "data/tpi1k.tif", datatype = "FLT4S")

## load above from ADA
TPI250 <- raster(file.path(gis, "r_utm", "tpi250.tif"))
TPI500 <- raster(file.path(gis, "r_utm", "tpi500.tif"))
TPI1k <- raster(file.path(gis, "r_utm", "tpi1k.tif"))

TPI250

plot(TPI250)
plot(TPI500)
plot(TPI1k)

plot(be10)

## make SD categories from TPI 500
sd <- cellStats(TPI500, "sd")
mn <- cellStats(TPI500, "mean")

# sd around the mean
breaks.1sd <- c(-Inf, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, Inf)*sd + mn
sd500.1sd <- cut(TPI500, breaks = breaks.1sd)

# weiss method + slope
# classes:       1   2     3   4   5
breaks <- c(-Inf, -1, -0.5, 0.5, 1, Inf)*sd + mn
breaks

sd500 <- cut(TPI500, breaks = breaks)
sd500
table(values(sd500))

plot(sd500, colNA = "black", col = rainbow(5))

sd500_6c <- raster::overlay(sd500, terr$slope, 
                             fun = function(x,z) {
                               ifelse((x == 5) & !is.na(z), 1, # ridge > 1 sd
                                      ifelse((x == 4) & !is.na(z), 2, # upper slope > 0.5 sd  =<1 sd
                                             ifelse((x == 3) & (z > 5), 3, # middle slope
                                                    ifelse((x == 3) & (z <= 5), 4, # flat slope
                                                           ifelse((x == 2) & !is.na(z), 5, # lower slope
                                                                  ifelse((x == 1) & !is.na(z), 6, NA)))))) # valleys
     
                                                       })
sd500_6c
plot(sd500_6c)
plot(stack(sd500.1sd, sd500, sd500_6c), col = rainbow(6))

# SiteName %in% c("076361", "159408") # two points not coverd

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
  twi <- dynatopmodel::upslope.area(be10, atb = T)
)
twi
# $atb is the twi
save(twi, file = file.path(gis, "r_utm/twi.rdata"))

# user  system elapsed  # on laptop
# 1299.66    4.90 1424.00

terr <- stack(be10, terr, Nss, Ess, twi$atb)
names(terr) <- c("be10","tri","slope","aspect","Nss","Ess","twi")


writeRaster(terr, bylayer = T, filename = file.path(gis, "r_utm/terr10.tif"), datatype = "FLT4S", suffix = "names", overwrite = TRUE)
writeRaster(covStack, bylayer = T, filename = file.path(gis, "r_utm/l.tif"), datatype = "FLT4S", suffix = "names")
save(terr, covStack, file = file.path(gis, "r_utm/elev_cov.rdata"))

# load(file.path(gis, "r_utm/elev_cov.rdata"))

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

### Extract values  ####

## Extract values and add to data frame
cut.r1k.pt <- extract(cutStack$cut_r1k, xy.utm)

## average elevation
be500 <- raster::extract(terr$be10, xy.utm, buffer = 500, fun=mean, na.rm = T)
# point elevation
dem.pt <- raster::extract(terr$be10, xy.utm)

## topographic index
mTopo <- dem.pt - be500

## All terr
terr.pt <- raster::extract(terr, xy.utm)

# lidar cover
cov.pt <- raster::extract(covStack, xy.utm)

## make data frame
topo.df <- data.frame(siteName = xy.utm$SiteName, terr.pt, cov.pt, be500 = be500, mTopo, cut.r1k.pt)
head(topo.df)
str(topo.df)
summary(topo.df)

save(topo.df, file = "Hmsc_CD/oregon_ada/data/topo_data.rdata")
save(topo.df, file = "HJA_scripts/10_eo_data/topo_data.rdata")

cat(paste(colnames(topo.df), collapse = '", "'))

pairs(topo.df[,c("be10", "slope", "Nss", "Ess", "twi", "ht", "ht.r1k", "cov2_4", "cov2_4.r1k", "cov4_16", "cov4_16.r1k", "mTopo", "cut.r1k.pt")])


# c("be10", "slope", "Ess", "twi", "ht", "cov2_4", "cov4_16", "mTopo") # , "cut.r1k.pt"
# c("be10", "slope", "Ess", "twi","ht.r1k", "cov2_4.r1k", "cov4_16.r1k", "mTopo") # , "cut.r1k.pt"


### ADd annual metrics 

# # stdDev bandNames
# stdN <- c("B1_stdDev", "B2_stdDev", "B3_stdDev", "B4_stdDev", "B5_stdDev", "B6_stdDev", "B7_stdDev", "B10_stdDev", "B11_stdDev", "ndvi_stdDev", "ndmi_stdDev", "nbr_stdDev")
# 
# ## q bandNames
# qN <- c("B1_p0", "B1_p5", "B1_p50", "B1_p95", "B1_p100", "B2_p0", "B2_p5", "B2_p50", "B2_p95", "B2_p100", "B3_p0", "B3_p5", "B3_p50", "B3_p95", "B3_p100", "B4_p0", "B4_p5", "B4_p50", "B4_p95", "B4_p100", "B5_p0", "B5_p5", "B5_p50", "B5_p95", "B5_p100", "B6_p0", "B6_p5", "B6_p50", "B6_p95", "B6_p100", "B7_p0", "B7_p5", "B7_p50", "B7_p95", "B7_p100", "B10_p0", "B10_p5", "B10_p50", "B10_p95", "B10_p100", "B11_p0", "B11_p5", "B11_p50", "B11_p95", "B11_p100", "ndvi_p0", "ndvi_p5", "ndvi_p50", "ndvi_p95", "ndvi_p100", "ndmi_p0", "ndmi_p5", "ndmi_p50", "ndmi_p95", "ndmi_p100", "nbr_p0", "nbr_p5", "nbr_p50", "nbr_p95", "nbr_p100")
# 
# # least cloudy
# cldN <- c("LC08_045029_20180726_B1","LC08_045029_20180726_B2","LC08_045029_20180726_B3","LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B6","LC08_045029_20180726_B7","LC08_045029_20180726_B10","LC08_045029_20180726_B11","LC08_045029_20180726_ndvi","LC08_045029_20180726_ndmi","LC08_045029_20180726_nbr","LC08_046029_20180818_B1","LC08_046029_20180818_B2","LC08_046029_20180818_B3","LC08_046029_20180818_B4","LC08_046029_20180818_B5","LC08_046029_20180818_B6","LC08_046029_20180818_B7","LC08_046029_20180818_B10","LC08_046029_20180818_B11","LC08_046029_20180818_ndvi","LC08_046029_20180818_ndmi","LC08_046029_20180818_nbr","LC08_045030_20180726_B1","LC08_045030_20180726_B2","LC08_045030_20180726_B3","LC08_045030_20180726_B4","LC08_045030_20180726_B5","LC08_045030_20180726_B6","LC08_045030_20180726_B7","LC08_045030_20180726_B10","LC08_045030_20180726_B11","LC08_045030_20180726_ndvi","LC08_045030_20180726_ndmi","LC08_045030_20180726_nbr")

## Get rasters

std <- brick(file.path(gis, "r_utm/gee/stdDev.tif"))
qnt <- brick(file.path(gis, "r_utm/gee/quantiles.tif"))
cld <- brick(file.path(gis, "r_utm/gee/leastCloud2.tif"))

# set NA value
NAvalue(std) <- -9999
NAvalue(qnt) <- -9999
NAvalue(cld) <- -9999

std
names(std)
names(qnt)
names(cld)
cellStats(std$ndvi_stdDev, range)
NAvalue(cld)

cellStats(cld$LC08_045029_20180726_B2, range)

plotRGB(cld, r = 4, g = 3, b = 2, stretch = "lin", colNA = "black")
plot(cld$LC08_045029_20180726_nbr, colNA = "black")

plot(std[[c("ndvi_stdDev", "savi_stdDev")]])
plot(qnt[[c("ndvi_p50", "savi_p50")]])

plot(qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
       "ndmi_p5", "ndmi_p50", "ndmi_p95",
       "savi_p5", "savi_p50", "savi_p95")]])



## Make subset stack

annualStack <- raster::stack(std[[c("ndmi_stdDev")]],
                     qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
                           "ndmi_p5", "ndmi_p50", "ndmi_p95",
                            "savi_p50")]],
                     cld[[c("LC08_045029_20180726_B1",
                            "LC08_045029_20180726_B3","LC08_045029_20180726_B4",
                            "LC08_045029_20180726_B5",
                            "LC08_045029_20180726_B7","LC08_045029_20180726_B10")]])

annualStack
xy.utm

annual.pts <- extract(annualStack, xy.utm)
pairs(annual.pts)

annual.100m <- extract(annualStack, xy.utm, buffer = 100, fun=mean, na.rm = T)

head(annual.pts)
head(annual.100m)
colnames(annual.100m) <- paste0(colnames(annual.100m), "_100m")

## Add TPI to data set
# tpi <- extract(sd500_6c, xy.utm)
tpi_stck <- stack(TPI250, TPI500, TPI1k, sd500_6c)
names(tpi_stck) <- c("tpi250", "tpi500", "tpi1k", "tpi")
tpi <- extract(tpi_stck, xy.utm)

sum(is.na(tpi))
which(is.na(tpi))

plot(sd500_6c)
plot(xy.utm[which(is.na(tpi)),], pch = 16, add = T)

plot(be10)
plot(xy.utm[which(is.na(tpi)),], pch = 16, add = T)

plot(TPI500)
plot(xy.utm[which(is.na(tpi)),], pch = 16, add = T)


annual.df <- cbind(annual.pts, annual.100m, tpi)
colnames(annual.df)
head(annual.df)

# Join with previous
load("Hmsc_CD/oregon_ada/data/topo_data.rdata") # topo.df
head(topo.df)

ann.topo <- data.frame(SiteName = topo.df$siteName, annual.df)
str(ann.topo)
head(ann.topo)

save(ann.topo, file = "Hmsc_CD/oregon_ada/data/ann_topo.df")

which(is.na(ann.topo$tpi))


ann.topo[1:10, 1:10]
pairs(ann.topo[, c(2,4,7:10, 24:26, 30:32)])

pairs(ann.topo[, c(18:30, 52)])

## export some plots

##
library(ggplot2)

load(file.path(gis, "r_utm/cut_stack.rdata")) #cutStack
load(file.path(gis, "r_utm/elev_cov.rdata")) # terr, covStack

png("Hmsc_CD/local/plots/cut_pred.png", width = 300, height = 300, units = "mm", res = 200)
plot(cutStack)
dev.off()
png("Hmsc_CD/local/plots/terr_pred.png", width = 300, height = 300, units = "mm", res = 200)
plot(terr)
dev.off()
png("Hmsc_CD/local/plots/cov_pred.png", width = 300, height = 300, units = "mm", res = 200)
plot(covStack)
dev.off()

## reduce resolution for plotting
covStack <- aggregate(covStack, 10) # 30 x 30 to ... 
plot(covStack)

terr <- aggregate(terr, 30)

# prepare for ggplot geom_raster
coords <- xyFromCell(terr, seq_len(ncell(terr)))
terr.df <- utils::stack(as.data.frame(getValues(terr)))
terr.df <- cbind(coords, terr.df)
head(terr.df)

coords <- xyFromCell(cutStack, seq_len(ncell(cutStack)))
cut.df <- cbind(coords, utils::stack(data.frame(values(cutStack))))
head(cut.df)

coords <- xyFromCell(covStack, seq_len(ncell(covStack)))
cov.df <- cbind(coords, utils::stack(data.frame(values(covStack))))
head(cov.df)

## do the points
xy <- data.frame(st_coordinates(xy.utm))
xy$ind <- unique(terr.df$ind)[1]
head(xy)

ggplot(terr.df) + 
  geom_tile(aes(x, y, fill = values)) +
  geom_point(aes(X,Y, size = 2), data = xy)+
  facet_wrap(~ ind) +
  scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal()


xy$ind <- unique(cov.df$ind)[1]
head(xy)

cov.df %>%
  filter(grepl("ht", ind)) %>%
  ggplot()+ 
  geom_tile(aes(x, y, fill = values)) +
  geom_point(aes(X,Y), size = 1, data = xy)+
  facet_wrap(~ ind) +
  scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal()
ggsave("Hmsc_CD/local/plots/ht_pred.png", width = 300, height = 300, units = "mm")


xy$ind <- unique(cut.df$ind)[1]
ggplot(cut.df) + 
  geom_tile(aes(x, y, fill = values)) +
  geom_point(aes(X,Y), data = xy)+
  facet_wrap(~ ind) +
  scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal()


#### Predictor plots

hllshd <- raster(file.path(gis, "r_utm/bE_30m_hlshd.tif"))

# prepare for ggplot geom_raster
coords <- xyFromCell(hllshd, seq_len(ncell(hllshd)))
hllshd.df <- utils::stack(as.data.frame(getValues(hllshd)))
hllshd.df <- cbind(coords, hllshd.df)
head(hllshd.df)


stck.sub <- annualStack[[c("ndmi_stdDev", "ndvi_p5", "ndvi_p50", "ndvi_p95")]]
# ## stack sa long form, for facet
coords <- xyFromCell(stck.sub, seq_len(ncell(stck.sub)))
stck.df <- utils::stack(as.data.frame(getValues(stck.sub)))
stck.df <- cbind(coords, stck.df)
head(stck.df)

ggplot(data = stck.df, aes(x = x, y = y, fill = values)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "B", name = "NDMI annual std dev") +
  facet_wrap(~ind)+
  coord_equal()
ggsave("local/plots/predictors_annual_facet.png", width = 300, height = 300, units = "mm")


stck.sub <- annualStack[[c("ndmi_stdDev", "ndvi_p5", "ndvi_p50", "ndvi_p95")]]
# ## stack wide - not for facet
coords <- xyFromCell(stck.sub, seq_len(ncell(stck.sub)))
stck.df <- data.frame(values(stck.sub))
stck.df <- cbind(coords, stck.df)
head(stck.df)

sapply(stck.df, range, na.rm=T)


p1 <- ggplot(data = stck.df, aes(x = x, y = y, fill = ndmi_stdDev)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "B", name = "NDMI annual std dev") +
  #geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill= "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  coord_equal()
  
p2 <- ggplot(data = stck.df, aes(x = x, y = y, fill = ndvi_p5)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "C", name = "5 percentile NDVI") +
  #geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill= "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  coord_equal()

p3 <- ggplot(data = stck.df, aes(x = x, y = y, fill = ndvi_p50)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "B", name = "Median NDVI") +
  #geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill= "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  coord_equal()

p4 <- ggplot(data = stck.df, aes(x = x, y = y, fill = ndvi_p95)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "C", name = "95 percentile NDVI") +
  #geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill= "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  coord_equal()
  
plot_grid(p1, p2, p3, p4, ncol = 2)
ggsave("local/plots/predictors_annual.png", width = 300, height = 300, units = "mm")
