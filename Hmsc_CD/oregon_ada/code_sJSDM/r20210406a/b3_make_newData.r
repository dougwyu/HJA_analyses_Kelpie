#### Make prediction grid for sjSDM

## ON LOCAL

getwd()
dir()

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

library(raster)
library(sf)

# wgs84 UTM 10N
utm10N <- 32610

# gis <- file.path(wd, "HJA_scripts/10_eo_data/raw_gis_data") 
gis <- "J:/UEA/Oregon/gis"

dir(gis)
# raster files are here:
dir(file.path(gis, "r_utm"))


## Load rasters

# get elevation raster
# bareearth, projected to UTM 10N wgs84, bilinear, at 30m res, converted to m (/0.30480061) in arcgis
be30 <- raster(file.path(gis, "r_utm/bareEarth_m_30m.tif"))

## TPI (made on ADA)
tpi_fn <- list.files(file.path(gis,"r_utm"), ".*tpi.*\\.tif$", full.names = TRUE)
tpi <- stack(tpi_fn)

terr30_fn <- list.files(file.path(gis,"r_utm"), ".*terr30.*\\.tif$", full.names = TRUE)
terr30 <- stack(terr30_fn)
terr30
names(terr30) <- sub("terr30_", "", names(terr30))

# get original lidar data
lidar_fn <- list.files(file.path(gis,"r_utm"), ".*lidar.*\\.tif$", full.names = TRUE)
lidar_fn
lidarStack <- stack(lidar_fn)
names(lidarStack)
names(lidarStack) <- sub("lidar_metric_mosaic_", "", names(lidarStack))

# get summarised lidar data by variable radius
cover_fn <- list.files(file.path(gis,"r_utm"), ".*l_.*\\.tif$", full.names = TRUE)
cover_fn
coverStack <- stack(cover_fn)
names(coverStack)
names(coverStack) <- sub("l_", "", names(coverStack))

# # load terrain and cover data
# load(file.path(gis, "r_utm/elev_cov.rdata"))
# rm(terr)
# covStack

load(file.path(gis, "r_utm/cut_stack.rdata"))
cutStack
names(cutStack)


## Annual data
std <- brick(file.path(gis, "r_utm/gee/stdDev.tif"))
qnt <- brick(file.path(gis, "r_utm/gee/quantiles.tif"))
cld <- brick(file.path(gis, "r_utm/gee/leastCloud2.tif"))

# set NA value
NAvalue(std) <- -9999
NAvalue(qnt) <- -9999
NAvalue(cld) <- -9999


names(std)
names(qnt)
names(cld)

## Make subset stack  for gee data (used in models)
annualStack <- raster::stack(std[[c("ndmi_stdDev")]],
                             qnt[[c("ndvi_p5", "ndvi_p50", "ndvi_p95", 
                                    "ndmi_p5", "ndmi_p50", "ndmi_p95",
                                    "savi_p50")]],
                             cld[[c("LC08_045029_20180726_B1",
                                    "LC08_045029_20180726_B3","LC08_045029_20180726_B4",
                                    "LC08_045029_20180726_B5",
                                    "LC08_045029_20180726_B7","LC08_045029_20180726_B10")]])


extent(annualStack)
extent(coverStack)
extent(terr30) == extent(be30)

extent(lidarStack)

extent(be30)
extent(cutStack)

extent(tpi)

as.vector(extent(tpi))

origin(annualStack)

origin(covStack)
origin(terr30)
origin(be30)
origin(tpi)

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/comExt.r")
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")

ext <- comExt(annualStack, covStack, terr30, be30, tpi)
ext

ext <- adjExt(ext, d = 30, expand = FALSE)
ext

# interpolate to same extent, origin and resolution
# make template raster

r <- raster(ext, res = c(30,30), crs = proj4string(tpi))
r
origin(r)

beginCluster()
annStack_rs <- raster::resample(annualStack, r)
be30_rs <- resample(be30, r)
tpi_rs <- resample(tpi, r)
terr30_rs <- resample(terr30, r)
coverStack_rs <- resample(coverStack, r)
lidarStack_rs <- resample(lidarStack, r)
cutStack_rs <- resample(cutStack,r)
endCluster()


## create rasters from vector data
hja <- st_read(file.path(gis, "s_nad_utm/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja_bound

hja.utm <- st_transform(hja_bound, crs = utm10N)
hja.utm

HJA <- rasterize(hja.utm, r)
names(HJA) <- "insideHJA"
plot(HJA)

## Years since disturb is in cutStack.. see p0_elevation_cover_predictors.r for this.

# MTOPO

# "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb"


allStack <- stack(annStack_rs, be30_rs, tpi_rs, terr30_rs, coverStack_rs, cutStack_rs, lidarStack_rs, HJA)
allStack
names(allStack)

## CHANGE RESOLUTION HERE/./// or above in resample

# get coords
xyFromCell()


# get xy data (could subset first by vars, but better to save all as matrix for now)

allData <- values(allStack)
dim(allData)
indNA <- rowSums(is.na(allData)) == ncol(allData) # TRUE where all NAs

# remove NAs, convert to data frame and save index to replace after prediction
allData <- data.frame(allData[!indNA, ])

## change categorical to predictor values
allData[, "insideHJA"] <- ifelse(is.na(allData[, "insideHJA"]), "no", "yes")
table(allData[,"insideHJA"])

# cut.r == YrsSinceDisturb

summary(allData)

save(allData, indNA, r, file = file.path("data", "predData.rdata"))


