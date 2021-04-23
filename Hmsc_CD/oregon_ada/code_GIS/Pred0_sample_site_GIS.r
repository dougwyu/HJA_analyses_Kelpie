
### Create spatial object from sample points ####

## ON LOCAL
getwd()
wd <- here::here()
setwd(wd)
dir()


library(sf)

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269

# gis_in <- gis_out <- "data/gis" # change this to github GIS folder
gis_in <- "HJA_scripts/10_eo_data/raw_gis_data"
gis_out <- "HJA_scripts/10_eo_data/processed_gis_data"

dir(gis_in)

# get points
samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",minimaprundate,"_",samtoolsfilter,"_", 
                                       samtoolsqual, "_kelpie", kelpierundate,"_", primer,"_vsearch97")

datFile <- paste0("sample_by_species_table_", samtoolsfilter, "_minimap2_", minimaprundate,"_kelpie",
                  kelpierundate,"_uncorr.csv")

otuenv <- read.csv(file.path(gitHub, outputidxstatstabulatefolder, datFile))

otuenv[1:6,1:10]
coords <- unique(otuenv[,c("SiteName", "UTM_E", "UTM_N")])
xy.sf <- st_as_sf(coords, coords = c("UTM_E", "UTM_N"), crs = nadutm10)

rm(gitHub, otuenv, outputidxstatstabulatefolder, datFile, primer, 
   kelpierundate, minimaprundate, samtoolsfilter, samtoolsqual, coords)

# transform to wgs utm to match rasters
xy.utm <- st_transform(xy.sf, crs = utm10N)
rm(xy.sf)

length(unique(xy.utm$SiteName))

# write
st_write(xy.utm, file.path(gis_out, "s_utm/sample_sites_utm10.shp"), delete_layer = T)
st_write(xy.utm, file.path(gis_out, "s_utm/sample_sites_utm10.kml"), delete_layer = T)
save(xy.utm, file = file.path(gis_out, "sample_sites.rdata"))

