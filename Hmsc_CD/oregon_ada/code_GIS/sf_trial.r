
options(echo=TRUE) # if you want see commands in output file
getwd() # always run sub from oregon_ada

library(raster)
library(sf)

pt1 = st_point(c(0,1))
pt2 = st_point(c(1,2))

st_sfc(pt1, pt2)

d = data.frame(a = 1:2)
d$geom = st_sfc(pt1, pt2)
df = st_as_sf(d)

print(df)

#plot(st_geometry(df))
raster(extent(df))
r_df <- rasterize(df, raster(extent(df)), field = "a")
#plot(r_df)
