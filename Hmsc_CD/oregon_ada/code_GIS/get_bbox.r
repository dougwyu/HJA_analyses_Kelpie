

library(sf)

## epsg codes
utm10N <- 32610
wgs <- 4326

sites.utm <- st_as_sf(env.vars[,c("UTM_E", "UTM_N", "SiteName", "uniqueID")], coords = c("UTM_E", "UTM_N"), crs = utm10N)
sites.utm

st_bbox(sites.utm)
# xmin    ymin    xmax    ymax 
# 555122 4891311  570775 4908110 

sites.wgs <- st_transform(sites.utm, crs = wgs)
st_bbox(sites.wgs)
# xmin       ymin       xmax       ymax 
# -122.30963   44.17273 -122.11404   44.32382 


# which is here....... 
bbox <- st_as_sfc(st_bbox(sites.wgs))
bbox_file <- tempfile(fileext = ".kml")
st_write(bbox, bbox_file)

shell.exec(bbox_file) # should open in google earth if installed... 


## Areas
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

st_area(hja.utm)/1000000
# 63.65363 km^2

st_area(st_as_sfc(st_bbox(sites.utm)))/1000000
# 262.954747 km^2
