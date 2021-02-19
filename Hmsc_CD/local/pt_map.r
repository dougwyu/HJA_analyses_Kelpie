### MAke map of points... 

getwd()
wd <- here::here()
setwd(wd)
dir()

## Updated to new vars, also changes to elevation_m, canopy_height_m  to _f. 

library(dplyr)

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",minimaprundate,"_",samtoolsfilter,"_", 
                                       samtoolsqual, "_kelpie", kelpierundate,"_", primer,"_vsearch97")

datFile <- paste0("sample_by_species_table_", samtoolsfilter, "_minimap2_", minimaprundate,"_kelpie",
                  kelpierundate,"_FSL_qp.csv")

fn <- file.path(gitHub, outputidxstatstabulatefolder, datFile)

# what file am i using?
basename(fn)

# when was it modified?
file.mtime(fn)

# read complete data set
otuenv <- read.csv(fn, stringsAsFactors = FALSE, na.strings = "NA")

## get predictors and coordiantes, sitenames

# clean up
rm(datFile, gitHub, kelpierundate, minimaprundate, outputidxstatstabulatefolder, primer, samtoolsfilter, samtoolsqual, fn)


# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates
env.vars <- otuenv %>% 
  dplyr::select(!contains("__"), -starts_with("nor")) %>%
  mutate(elevation_m = elevation_f * 0.3048, ## convert to metres
         canopyHeight_m = canopyHeight_f * 0.3048,
         lg_DistStream = log(distToStream_m + 0.001),
         lg_DistRoad = log(distToRoad_m + 0.001),
         lg_YrsDisturb = log(YrsSinceDist + 0.001)) %>%
  dplyr::select(SiteName, UTM_N, UTM_E, clearcut,insideHJA,oldGrowthIndex, elevation_m, canopyHeight_m, 
                precipitation_mm, minT_annual, maxT_annual, mean.NDVI, mean.EVI, mean.green, mean.wet, 
                mean.bright, l_p25, l_p95, l_rumple, lg_DistStream, lg_DistRoad, lg_YrsDisturb, be10, tri, 
                slope, twi, Nss, Ess, ht, ht.r500, cov2_4, cov2_4.r500, cov4_16, cov4_16.r500, mTopo, cut.r1k.pt) %>%
  mutate( #across(where(is.numeric), scale), # scale here # scale when defining models etc.
    clearcut = factor(clearcut),
    insideHJA = factor(insideHJA)) %>%
  distinct(SiteName, .keep_all = TRUE)
  
head(env.vars)

library(sf)

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269

pts.sf <- st_as_sf(env.vars, coords = c("UTM_E", "UTM_N"), crs = nadutm10)

pts.sf

# export all points to shapefile / kml
st_write(pts.sf, "J:/UEA/Oregon/gis/s_nad_utm/all_pts.shp")

# transform to wgs and export to kml
pts.wgs <- st_transform(pts.sf, crs = 4326)
st_write(pts.wgs, "J:/UEA/Oregon/gis/s_wgs/all_pts.kml")

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data

hja <- st_read("J:/UEA/Oregon/gis/s_nad_utm/HJA_Boundary.shp")

hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja_bound

hja.wgs <- st_transform(hja_bound, crs = 4326)
st_write(hja.wgs, "J:/UEA/Oregon/gis/s_wgs/hja.kml")


vars <- c("SiteName", "clearcut","insideHJA","oldGrowthIndex", "elevation_m", "canopyHeight_m",
"precipitation_mm", "minT_annual", "maxT_annual", "mean.NDVI",
"Nss", "Ess", "ht.r500", "cov2_4.r500", "cov4_16.r500", "mTopo", "cut.r1k.pt")

labs <- pts.sf$SiteName

mv1 <- mapview::mapview(pts.sf, zcol = vars, homebutton = FALSE,
              map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT"),
              legend = F)
mv1

## hide the layers by default.. modify leaflet object... 
mv1@map <- mv1@map %>% leaflet::hideGroup(group = vars[-1])

# do boundary map
mv2 <- mapview::mapview(hja_bound, col.regions = "blue", alpha.regions = 0.2, label = FALSE,
                        map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT"), legend= F, homebutton = TRUE)
mv2

mv <- mv1 + mv2
mv

# save
mapview::mapshot(mv, url = "Hmsc_CD/local/plots/pts_map.html", selfcontainted = TRUE)
# only works in wd at moment. moved from there.
mapview::mapshot(mv, url = "pts_map.html", selfcontainted = TRUE)
getwd()

