### Elevation data #####

wd <- here::here()
wd
setwd(wd)
dir()


library(sf)
library(raster)

## resolving mapview issues:
remotes::install_github("r-spatial/mapview")
webshot::install_phantomjs() # 'https://github.com/wch/webshot/releases/download/v0.3.1/phantomjs-2.1.1-windows.zip'
# phantomjs has been installed to C:\Users\55116479\AppData\Roaming\PhantomJS
library(mapview)
citation("mapview")


# EPSG:32610 # WGS 84 / UTM zone 10N
utm10N <- 32610
# EPSG:4269 # NAD 83
nad83 <- 4269

# import site data
# (follows quantiles of prevalence - as in models of 20201209 - see S2_define_models.r)

# get data
source("Hmsc_CD/local/L1_read_data.r")
rm(Y.train.pa, Y.train.qp, P)

head(S.train)
xy.sf <- st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = utm10N)

# transform to NAD (for now, quicker than projecting raster)
xy.nad <- st_transform(xy.sf, crs = nad83)

# import raster elevation data
terr <- brick("J:/UEA/Oregon/gis/r_nad/terrain.tif")
terr

names(terr) <- c("dem", "tri")
plot(terr)

# extract point elevation data
dem.pts <- raster::extract(terr, xy.nad)
colnames(dem.pts) <- c("dem.pt", "tri.pt")
head(dem.pts)

dem.500m <- raster::extract(terr, xy.nad, buffer = 500, fun = mean) # units in m for unprojected data
head(dem.500m)

# add elevation data per point to xy
xy.nad <- cbind(xy.nad, dem.pts)
xy.nad$dem500 <- dem.500m[,"dem"]
xy.nad$diffDem <- abs(xy.nad$dem500 - xy.nad$dem.pt)

plot(xy.nad$tri.pt, xy.nad$diffDem)

xy.nad

hist(xy.nad$dem.pt)
hist(xy.nad$diffDem)

save(xy.nad, file = "Hmsc_CD/local/xy_nad.rdata")
load("Hmsc_CD/local/xy_nad.rdata")

# library(mapview)
#library(leafpop)

mv <- mapview(xy.nad, zcol = c("diffDem", "tri.pt"),
                 map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT"),
                 legend = T)
mv

# save
mapviewOptions(fgb = FALSE)
mapshot(mv, url = "Hmsc_CD/local/plots/dem_map.html", selfcontainted = TRUE) # only works in wd at moment. moved from there.
# mapshot(mv, url = "dem_map.html", selfcontainted = TRUE)


mv2 <- mapview::mapview(xy.nad, zcol = c("diffDem"),
                       map.types = c("Esri.WorldImagery"),
                       legend = T)

x1 <- ~mapshot(mv2, file = "Hmsc_CD/local/plots/dem_map.png")




