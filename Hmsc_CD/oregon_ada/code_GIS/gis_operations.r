
setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(ggplot2)

utm10N <- 32610

## local
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"

## Edited study area
### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))
aoi.pred.sf_edit <- st_make_valid(aoi.pred.sf_edit)


## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

## disturbance
cut.sf <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut.utm <- st_transform(cut.sf, crs = utm10N)
plot(st_geometry(cut.utm))
plot(hja.utm, add = T, col = NA, border = "blue")

median(st_area(cut.utm))
sqrt(median(st_area(cut.utm)))
hist(st_area(cut.utm))

sd(st_area(cut.utm))

hist(cut.sf$YrsSinceDi)
table(cut.sf$YrsSinceDi, useNA = "always")

#
# ## Less than 40 years since disturbance
cut.40 <- subset(cut.utm, YrsSinceDi < 40)
hist(cut.40$YrsSinceDi)


# ## clip to prediction area
cut.intsct <- st_intersects(cut.40, aoi.pred.sf_edit)
cut.40 <- cut.40[lengths(cut.intsct)>0,]

#plot(st_geometry(aoi.pred.sf_edit))
#plot(cut.40, add =T)

# clip to prediction area
tmp <- st_difference(aoi.pred.sf_edit, sf::st_union(cut.40))
cut.40 <- st_difference(aoi.pred.sf_edit, tmp)

st_write(cut.40, file.path(gis_out, "s_utm/cut40.shp"))


