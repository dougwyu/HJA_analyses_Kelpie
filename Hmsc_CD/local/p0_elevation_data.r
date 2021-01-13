
#### get elevation data and average over radius around points

library(dplyr)
library(here)
library(glue)
library(sf)
library(raster)

here()
setwd(here())

# get XY data
samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

outputidxstatstabulatefolder <- glue::glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")

otuenv <- read.csv(here("Kelpie_maps", 
                        outputidxstatstabulatefolder, glue("sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_uncorr.csv")))

# M1S1
trap <- "M1"
period <- "S1"
otuenv <- otuenv %>% 
  filter(trap == trap[[1]] & period == period[[1]]) 

XY <- otuenv %>% 
  dplyr::select(UTM_E, UTM_N, SiteName, trap, period) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_"))
  #scale() %>%  # keep orignal for plotting

head(XY)  

# HJ Andrews Experimental Forest. # https://www.davidbuckleyborden.com/hja-experimental-forest
# This is UTM zone: EPSG:3717 # NAD83 / UTM zone 10N

# projection information (EPSG id)
utm10N <- 3717 # NAD83(NSRS2007) UTM 10 N
nad83 <- 4269
# wgs <- 4326

# points as sf and transform to utm and NAD
head(XY)
xy.utm <- st_as_sf(XY, coords = c("UTM_E", "UTM_N"), crs = utm10N)
xy.nad <- st_transform(xy.utm, crs = nad83)
# save
# st_write(xy.nad, "Hmsc_CD/gis/s_nad/m1s 1_points.shp", delete_layer = T)
# for google earth
# xy.wgs <- st_transform(xy.utm, crs = wgs)
# st_write(xy.wgs, "Hmsc_CD/gis/m1s1_points.kml")

## Crop, project, DEM (USGS Elevation 13), calculate roughness and save as multilayer tif.

## elevation raster tile downloaded from here. Stored off github... 350 MB
# https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/n45w123/USGS_13_n45w123.tif
# DEM at ~ 10 m resolution.. 1/3 arc seconds at 45 degrees latitude: 
lat.m <- 111131; lon.m <- 78846
c = sqrt((lat.m/60/60/3)^2 + (lon.m/60/60/3)^2)
sqrt(c^2/2) # approx resolution in m

# 
# # get bounding box of points to crop raster to (to save on processing time for projection) - in raster::extent format
# st_bbox(xy.nad)
# # adjust extent outwards
# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")
# (crop.bb <- adjExt(st_bbox(xy.nad), d = 0.1, outF = "Extent"))
# 
# # get raster
# dem <- raster("J:/UEA/Oregon/gis/r_nad/USGS_13_n45w123.tif")
# dem
# 
# # crop raster to extent of points
# dem.crp <- raster::crop(dem, crop.bb)
# dem.crp
# 
# # # project to utm
# # st_crs(xy.utm)[[2]]
# # 
# # # project (30 m resolution)
# # beginCluster()
# # dem.utm <- projectRaster(dem, crs = CRS(st_crs(xy.utm)[[2]]), res = 10) # at 10 m resolution
# # dem.utm
# # endCluster()
# # 
# # names(dem.utm) <- "dem"
# 
# # get roughness and other terrain metrics
# terrain <- raster::terrain(dem.crp, opt = c("TRI")) #, units = "degrees"
# terrain
# 
# ## save as raster
# writeRaster(stack(dem.crp, terrain), file="Hmsc_CD/gis/r_nad/terrain.tif",
#             datatype = "FLT4S", overwrite = T)


# Import dem, roughness raster
dem.stck <- brick("Hmsc_CD/gis/r_nad/terrain.tif")
names(dem.stck) <- c("dem", "tri")
dem.stck

plot(dem.stck$tri, breaks = seq(0,10,1), col = terrain.colors(10))
plot(st_geometry(xy.nad), add = T, pch = 16, cex = 0.7)


## how wide a radius for mean elvataion?
plot(st_geometry(xy.utm), pch = 19)

# what's min distance between points?
dist <- st_distance(xy.utm)
diag(dist) <- NA
min(dist, na.rm = T)

plot(st_geometry(st_buffer(xy.utm, 500)), add = T) # 500 m radius buffer. 
plot(st_geometry(st_buffer(xy.utm, 1000)), add = T, border= "red")


# get mean values over buffer m radius (projects internally)
tri.pt <- raster::extract(dem.stck$tri, xy.nad)
dem500<- raster::extract(dem.stck$dem, xy.nad, buffer = 500, fun=mean, na.rm = T)

# join to point ID
dem_stats <- data.frame(XY, tri.pt, dem500)
head(dem_stats)

pairs(dem_stats[,c("tri.pt", "dem500")])

save(dem_stats, file = "Hmsc_CD/gis/demStats.rdata")
