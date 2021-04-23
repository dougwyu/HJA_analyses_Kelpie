
## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()


options(echo=TRUE) # if you want see commands in output file
getwd() # always run sub from oregon_ada

library(raster)

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
#plot(ht_gt4_utm)

# ncell(r1)/(30*30)
# r30

# craete windows
d250 <- (focalWeight(ht_gt4_utm, d = 250, type = "circle")>0)*1
d500 <- (focalWeight(ht_gt4_utm, d = 500, type = "circle")>0)*1

dim(d250)

# do focal over 250 windows
n <- 12
beginCluster(n = n)

ht_gt4_f250 <- focal(ht_gt4_utm, w = d250, fun = function(x) sum(x)/sum(x!=0), pad = TRUE, padValue = 0,
                   filename = file.path(gis, "ht_gt4_f250.tif"), datatype = "FLT4S")

endCluster()
# 


