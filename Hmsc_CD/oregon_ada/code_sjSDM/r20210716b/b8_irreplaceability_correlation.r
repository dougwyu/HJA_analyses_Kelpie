




### Summarise by logging plots

## Bring in logging and HJA data
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


## Edited study area
### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))
# aoi.pred.sf_edit <- st_make_valid(aoi.pred.sf_edit)

# clip cut to model area and convert to raster
# cut.r <- rasterize(cut.utm, r.msk, field = "YrsSinceDi")
# cut.r <- mask(cut.r, r.msk)
# plot(cut.r)

cut.ints <- cut.utm[lengths(st_intersects(cut.utm, aoi.pred.sf_edit)) > 0, ]

par(mfrow= c(1,1))
plot(beta.r.prob)
plot(cut.ints, col = NA, border = "grey60", add =T)
plot(hja.utm, add = T, border = "blue", col = NA)
# plot(cut.ints[cut.ints$TREATMENT_ == "Clearcut",], col = NA, border = "blue", add =T)

cut_stats <- extract(beta.r, cut.ints, fun = mean, na.rm = T)
cut_stats <- extract(beta.r, cut.ints, fun = max, na.rm = T)

head(cut_stats)

cut.ints$irr_mn <- cut_stats

plot(cut.ints$YrsSinceDi, cut.ints$irr_mn)

plot(cut.ints[, "irr_mn"])

hist(values(beta.r))
plot(beta.r)
plot(hja.utm, add = T, border = "blue", col = NA)
plot(cut.ints, col = NA, border = "grey", add =T)
# plot(cut.ints[cut.ints$TREATMENT_ == "Clearcut",], col = NA, border = "blue", add =T)

