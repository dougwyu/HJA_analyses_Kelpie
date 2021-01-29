
getwd()
wd <- here::here()
setwd(wd)
dir()


library(sf)
library(raster)


# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N # original coordinates of sample poitns are in this proejction - Marie
nadutm10 <- 26910


gis <- file.path(wd, "HJA_scripts/10_eo_data/raw_gis_data") # originally stored locally
# gis is organised by shape (vector) / raster folders, according to projection
dir(gis)
# raster files are here:
dir(file.path(gis, "r_utm"))

### Load data ###

# load sample points:
# get data
source("Hmsc_CD/local/L1_read_data_v3.r") # reads all data, covariates, etc. from github
rm(Y.train.pa, Y.train.qp, P, otu.pa.csv, otu.qp.csv, otu.ab.csv, X.train, topo.df, otuenv)

head(S.train)
# convert to sf spatial data frame
xy.sf <- st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = nadutm10)
rm(S.train)

# transform to wgs utm to match rasters
xy.utm <- st_transform(xy.sf, crs = utm10N)
rm(xy.sf)
xy.utm


# load p95 lidar raster
covStack <- stack(file.path(gis, "r_utm", c("lidar_metric_mosaic_p95.tif", 
                                       "lidar_metric_mosaic_Cover_2m_4m.tif", 
                                       "lidar_metric_mosaic_Cover_4m_16m.tif")))

names(covStack) <- c("p95", "cov2_4", "cov4_16")
covStack

# extract lidar data at 
cov.pts <- st_as_sf(raster::extract(covStack, xy.utm, sp = T)) # extract, reconvert to sf to get all columns in result
head(cov.pts)
#

head(env.vars)

mData <- env.vars %>%
  dplyr::select(uniqueID, l_p95, l_Cover_2m_4m, l_Cover_4m_16m) %>%
  dplyr::left_join(y = cov.pts) %>%
  as.data.frame()

head(mData)

mData %>%
  dplyr::select(where(is.numeric)) %>%
  pairs()

## OK< all good.

