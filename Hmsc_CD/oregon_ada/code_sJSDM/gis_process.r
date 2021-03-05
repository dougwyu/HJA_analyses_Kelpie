
library(raster)

# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# get elevation raster
# bareearth, projected to UTM 10N wgs84, bilinear, at 10m res, converted to m (/0.30480061) in arcgis
be10 <- raster("data/bareEarth_m_10m.tif")
be10



d9 <- matrix(1,3,3)
d100 <- (focalWeight(be10, d = 100, type = "circle")>0)*1

d250 <- (focalWeight(be10, d = 250, type = "circle")>0)*1
dim(d250)
d500 <- (focalWeight(be10, d = 500, type = "circle")>0)*1
d1k <- (focalWeight(be10, d = 1000, type = "circle")>0)*1
d2k <- (focalWeight(be10, d = 2000, type = "circle")>0)*1

## Do TPI at three different scales
c100 <- ceiling(prod(dim(d100))/2)
c250 <- ceiling(prod(dim(d250))/2)
c500 <- ceiling(prod(dim(d500))/2) # get central cell
c1k <- ceiling(prod(dim(d1k))/2) # get central cell
c2k <- ceiling(prod(dim(d2k))/2) # get central cell

TPI100 <- focal(be10, w=d100, fun=function(x, ...) x[c100] - mean(x[-c100]), pad=TRUE, padValue=NA)
plot(TPI100)

TPI250 <- focal(be10, w=d250, fun=function(x, ...) x[c250] - mean(x[-c250]), pad=TRUE, padValue=NA)
# writeRaster(TPI250, filename = "data/tpi250.tif", datatype = "FLT4S")
# 
TPI500 <- focal(be10, w=d500, fun=function(x, ...) x[c500] - mean(x[-c500]), pad=TRUE, padValue=NA)
# writeRaster(TPI500, filename = "data/tpi500.tif", datatype = "FLT4S")
# 
TPI1k <- focal(be10, w=d1k, fun=function(x, ...) x[c1k] - mean(x[-c1k]), pad=TRUE, padValue=NA)
# writeRaster(TPI1k, filename = "data/tpi1k.tif", datatype = "FLT4S")


# gis <- "J:/UEA/Oregon/gis"
# be10 <- raster(file.path(gis, "r_utm/bareEarth_m_30m.tif"))
# TPI2k <- focal(be10, w=d2k, fun=function(x, ...) x[c2k] - mean(x[-c2k]), pad=TRUE, padValue=NA)
# plot(TPI2k)
# 
# writeRaster(TPI2k, filename = file.path(gis, "r_utm/tpi2k.tif"))

# m250 <- focal(be10, w=d250, fun=mean, pad=TRUE, padValue=NA)
# writeRaster(m250, filename = "data/m250.tif", datatype = "FLT4S")
# 
# m500 <- focal(be10, w=d500, fun=mean, pad=TRUE, padValue=NA)
# writeRaster(m500, filename = "data/m500.tif", datatype = "FLT4S")
# 
# m1k <- focal(be10, w=d1k, fun=mean, pad=TRUE, padValue=NA)
# writeRaster(m1k, filename = "data/m1k.tif", datatype = "FLT4S")



