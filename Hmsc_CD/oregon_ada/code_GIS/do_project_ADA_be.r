options(echo=TRUE) # if you want see commands in output file
getwd() # always run sub from oregon_ada

library(raster)


# wgs84 UTM 10N
utm10N <- 32610
prj4.utm10 <- "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs" #same as above for raster (but now also accepts epsg)
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910

#testing local
gis_in <- gis_out <- "data/gis" # change this to github GIS folder
# gis_in <- "HJA_scripts/10_eo_data/raw_gis_data"
# gis_out <- "HJA_scripts/10_eo_data/processed_gis_data"
dir(gis_in)

## 0. Get prelim data ######
## common extent in utm10N
load(file.path(gis_in, "commonExtent.rdata")) # ext
ext

# CHECK
# r10
# origin(r10)
# r30
# origin(r30)
# extent(r10) == extent(r30)
# r1
# extent(r1) == extent(r10)


### 1. Canopy cover ####

# Create canopy cover (from Lidar) - subtract highest hit from bare earth.
# Create canopy gap/cover binary layer (<4 m is a canopy gap), summarise proportion of canopy cover over 250 and 500 m

# hh <- raster(file.path(gis_in, "latlong_highesthit.tif")) # from Lidar data. Provided by Oregon State University

be <- raster(file.path(gis_in, "latlong_bare_earth.tif"))
# crs: +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs
# at 3.014 ft resolution approx 1m


# test on small data set
# crp.ext <- c(xmin = 850000, xmax = 855000, ymin = 900000, ymax = 905000)
# hh <- crop(hh, crp.ext)
# hh
# be <- crop(be, crp.ext)


## set up raster parallel
# # set number of cores...
# n <- parallel::detectCores() - 1
# n <- 6
# n <- 2

# beginCluster(n = n)


# convert to m
# f1 <- function(x) x * 0.30480

# hh_m <- clusterR(hh, calc, args=list(fun=f1))
# be_m <- clusterR(be, calc, args=list(fun=f1))

# hh <- hh * 0.30480
be_m <- be * 0.30480

# ## canopy height
# s <- stack(hh_m, be_m)
# f2 <- function(h, b) h - b
# ht_m <- clusterR(s, overlay, args = list(fun = f2))
# # ht_m <- hh_m - be_m
# 
# # set values below 0 to 0
# ht_m <- clusterR(ht_m, reclassify, args=list(rcl = c(-Inf, 0, 0)))
# # ht_m[ht_m < 0] <- 0
# ht_m
# # plot(ht)
# 
# ## create gap layer
# ht_gt4m <- clusterR(ht_m, reclassify, args=list(rcl = c(0, 4, 0, 4,Inf,1)))
# # ht_gt4m <- ht_m > 4
# ht_gt4m
# # plot(ht_gt4m)
# 
# ## project to utm10 and resolution of 1m, with origin same as 30 m rasters.

# small test extent
## ext <- extent(projectExtent(ht_m, crs = prj4.utm10))

## CHECK projected extents...
# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/ext2poly.r")
# ext.sf <- ext2poly(r30)
# ext.htgt4 <- ext2poly(ht_gt4m)
# ext.htgt4_utm <- st_transform(ext.htgt4, crs = st_crs(ext.sf))
#
# library(ggplot2)
#
# ggplot()+
#   geom_sf(data = ext.sf, alpha = 0.5, col = "blue", bg = NA)+
#   geom_sf(data = ext.htgt4_utm, alpha = 0.5, col = "red", bg = NA)+
#   coord_sf()

# make template rasters at these resolutions
r1 <- raster(ext, res = c(1,1), crs = prj4.utm10)
r10 <- raster(ext, res = c(10,10), crs = prj4.utm10)
r30 <- raster(ext, res = c(30,30), crs = prj4.utm10)

# beginCluster(n = n)

# # project canopy cover/gap layer
be_utm <- projectRaster(be_m, r1, method = "bilinear", filename = file.path(gis_out, "r_utm/be_utm.tif"),
                        datatype = "FLT4S", overwrite = TRUE)

#endCluster()


## get 30m resolution raster with mean of 1m raster for canopy cover and canopy height
be_r30 <- aggregate(be_utm, fact = 30, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/be_r30.tif"), datatype = "FLT4S", overwrite = TRUE)

be_r10 <- aggregate(be_utm, fact = 10, fun = mean, expand = FALSE,
                    filename = file.path(gis_out, "r_utm/be_r10.tif"), datatype = "FLT4S", overwrite = TRUE)
