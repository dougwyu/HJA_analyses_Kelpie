
## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()


options(echo=TRUE) # if you want see commands in output file
getwd() # always run sub from oregon_ada

library(raster)
library(sf)

# wgs84 UTM 10N
utm10N <- 32610
prj4.utm10 <- "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"

#testing local
# gis <- "../../../../Oregon/gis/marie"
gis <- "data/gis"
dir(gis)


# # canopy cover layer
# # ### Do canopy cover
# hh <- raster(file.path(gis, "latlong_highesthit.tif"))
# hh
# be <- raster(file.path(gis, "latlong_bare_earth.tif"))
# crs: +proj=lcc +lat_0=41.75 +lon_0=-120.5 +lat_1=43 +lat_2=45.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=ft +no_defs
# at 3.014 ft resolution (!!!) approx 1m

# test on small data set
# crp.ext <- c(xmin = 850000, xmax = 855000, ymin = 900000, ymax = 905000)
# hh <- crop(hh, crp.ext)
# hh
# be <- crop(be, crp.ext)

## canopy height
# ht <- hh - be

# set values below 0 to 0
# ht[ht < 0] <- 0
# ht
# plot(ht)

## create gap layer
# ht_gt4m <- ht > 13.1234
# ht_gt4m
# plot(ht_gt4m)

## project to utm10 and resolution of 1m, with origin same as 30 m rasters.
## common extent in utm10N
load(file.path(gis, "commonExtent.rdata")) # ext

# small test extent
# ext <- extent(projectExtent(ht, crs = prj4.utm10))

# make template rasters at these resolutions
r1 <- raster(ext, res = c(1,1), crs = prj4.utm10)
r10 <- raster(ext, res = c(10,10), crs = prj4.utm10)
r30 <- raster(ext, res = c(30,30), crs = prj4.utm10)

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



# CHECK
# r10
# origin(r10)
# r30
# origin(r30)
# extent(r10) == extent(r30)
#
# r1
# extent(r1) == extent(r10)

# # set number of cores...
# n <- 12
# # n <- parallel::detectCores() - 1
# 
# beginCluster(n = 12)
# ht_gt4_utm <- projectRaster(ht_gt4m, r1, method = "ngb", filename = file.path(gis, "ht_gt4_utm.tif"), 
#                             datatype = "INT1U", overwrite = TRUE)
# endCluster()
# 
# ht_gt4_utm

ht_gt4_utm <- raster(file.path(gis, "ht_gt4_utm.tif"))
ht_gt4_utm
dataType(ht_gt4_utm)
#plot(ht_gt4_utm)

# ncell(r1)/(30*30)
# r30

## get 30 resolution raster with mean of 1m raster
# ht_gt4_r30 <- aggregate(ht_gt4_utm, fact = 30, fun = mean, expand = FALSE,
#                         filename = file.path(gis, "ht_gt4_r30.tif"), datatype = "FLT4S", overwrite = TRUE)
# 
# # ht_gt4_r30 <- raster(file.path(gis, "ht_gt4_r30.tif"))
# # plot(ht_gt4_r30)
# rm(ht_gt4_r30); gc()
# 


# do focal windows of 250 and 500 m radius - but based on 30m resolution template.
# get coordinates at 30m resolution grid
xy_r30 <- raster::coordinates(r30)

# try with polygons?
head(data.frame(xy_r30, id = 1:nrow(xy_r30)))
# 
# library(sf)
pts_r30 <- st_as_sf(data.frame(xy_r30, id = 1:nrow(xy_r30))[50000:55000,], 
                    coords = c("x", "y"), crs = utm10N)
# pts_r30 <- st_as_sf(data.frame(xy_r30, id = 1:nrow(xy_r30)), coords = c("x", "y"), crs = utm10N)
pts_r30_b250 <- st_buffer(pts_r30, dist = 250)
pts_r30_b250

# subset  for testing - see time
# pts_r30_b250 <- pts_r30_b250[50000:55000,]

# plot(pts_r30_b250[1,], add = T)

n <- 4
# n <- 2
# 
# extract over radius at 30m res coords
beginCluster(n = n)
vals_gt4_r250 <- raster::extract(ht_gt4_utm, pts_r30_b250, method = "simple", fun = mean, na.rm = T)
endCluster()

save(vals_gt4_r250, file = file.path(gis, "vals_gt4_r250.rdata"))

# get cell numbers and extract?
## tabulate raster


# extract over radius at 30m res coords
# vals_gt4_r250 <- raster::extract(ht_gt4_utm, xy_r30, buffer = 250, method = "simple", fun = mean)

## put back into rasters
gt4_r250 <- r30
values(gt4_r250)[50000:55000] <- vals_gt4_r250
writeRaster(gt4_r250, filename = file.path(gis, "gt4_r30_r250.tif"), datatype = "FLT4S", overwrite = TRUE)
# rm(vals_gt4_r250, gt4_r250)
# gc()

# plot(gt4_r250)

# # same at 500
# vals_gt4_r500 <- raster::extract(ht_gt4_utm, xy_r30, buffer = 500, method = "simple", fun = mean)
# gt4_r500 <- r30
# values(gt4_r500) <- vals_gt4_r500
# writeRaster(gt4_r500, filename = file.path(gis, "gt4_r500.tif"), datatype = "FLT4S", overwrite = TRUE)
# rm(vals_gt4_r500, gt4_r500)
# gc()

#plot(gt4_r250)

