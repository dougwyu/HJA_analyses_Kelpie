### Plot residuals on map



## RUNS ON LOCAL

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
wd <- here::here()
wd
setwd(file.path(wd, "Hmsc_CD"))
dir()

library(ggplot)
library(raster)
library(cowplot)

getwd()

gis <- "J:/UEA/Oregon/gis"
# gis <- file.path(wd, "HJA_scripts/10_eo_data/raw_gis_data") 

load(file.path(resFolder, "residualData.rdata"))
# pred.sp, pred.nsp, resids.sp, resids.nsp, resids.sf,



## load hillshape and elevation

hllshd <- raster(file.path(gis, "r_utm/bE_30m_hlshd.tif"))
be30 <- raster(file.path(gis, "r_utm/bareEarth_m_30m.tif"))

str(getValues(hllshd))
range(getValues(hllshd), na.rm = T)

# prepare for ggplot geom_raster
coords <- xyFromCell(hllshd, seq_len(ncell(hllshd)))
hllshd.df <- utils::stack(as.data.frame(getValues(hllshd)))
hllshd.df <- cbind(coords, hllshd.df)
head(hllshd.df)

range(hllshd.df$values, na.rm=T)


coords <- xyFromCell(be30, seq_len(ncell(be30)))
be30.df <- utils::stack(as.data.frame(getValues(be30)))
be30.df <- cbind(coords, be30.df)
head(be30.df)
range(be30.df$values, na.rm=T)

## prepare points
xy <- cbind(st_coordinates(resids.sf), data.frame(resids.sf))
xy$geometry <- NULL
head(xy)


# can do as stack... but difficult to subset...  as ggplot gets the full range of values to scale 
#  even in subset
# be <- stack(hllshd, be30)
# coords <- xyFromCell(be, seq_len(ncell(be)))
# be.df <- utils::stack(as.data.frame(getValues(be)))
# be.df <- cbind(coords, be.df)
# head(be.df)
# 
# unique(be.df$ind)
# 
# ggplot(be.df) + 
#   geom_tile(aes(x, y, fill = values)) +
#   facet_wrap(~ ind) +
#   scale_fill_gradientn(colours = rev(terrain.colors(225))) +
#   coord_equal()
# 
# 
# ggplot() + 
#   geom_tile(data = subset(be.df, ind = bE_30m_hlshd), aes(x, y, alpha = values)) +
#   scale_alpha(range = c(0.25, 0.55), guide = "none")+
#   coord_equal()
# ggplot() + 
#   geom_tile(data = subset(be.df, ind = bE_30m_hlshd), aes(x, y, fill = values)) +
#   scale_fill_gradientn( colours = grey.colors(50, start = 0.3, end = 0.8))+
#   coord_equal()

ggplot() + 
  geom_tile(data = hllshd.df, aes(x, y, fill = values)) +
  scale_fill_gradientn( colours = grey.colors(50, start = 0.3, end = 0.8))+
  coord_equal()

ggplot() + 
  geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill= "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  coord_equal()



p1 <- ggplot() + 
  geom_tile(data = be30.df, aes(x, y, fill = values)) +
  geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill = "grey20") +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "Elevation") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_point(aes(X,Y, col = mean.nsp), size = 2, alpha = 0.5, data = xy)+
  scale_color_viridis_c(option = "A", name = "Residuals")+
  coord_equal()+
  ggtitle("Mean residuals - non spatial model")


p2 <- ggplot() + 
  geom_tile(data = be30.df, aes(x, y, fill = values)) +
  geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill = "grey20") +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "Elevation") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_point(aes(X,Y, col = mean.sp), size = 2, alpha = 0.5, data = xy)+
  scale_color_viridis_c(option = "A", name = "Residuals")+
  coord_equal()+
  ggtitle("Mean residuals - spatial model")


p3 <- ggplot() + 
  geom_tile(data = be30.df, aes(x, y, fill = values)) +
  geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill = "grey20") +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "Elevation") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_point(aes(X,Y, col = sd.nsp), size = 2, alpha = 0.5, data = xy)+
  scale_color_viridis_c(option = "A", name = "Residuals")+
  coord_equal()+
  ggtitle("Std dev residuals - non spatial model")

p4 <- ggplot() + 
  geom_tile(data = be30.df, aes(x, y, fill = values)) +
  geom_tile(data = hllshd.df, aes(x, y, alpha = values), fill = "grey20") +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "Elevation") +
  scale_alpha(range = c(0.15, 0.65), guide = "none")+
  geom_point(aes(X,Y, col = sd.sp), size = 2, alpha = 0.5, data = xy)+
  scale_color_viridis_c(option = "A", name = "Residuals")+
  coord_equal()+
  ggtitle("Sts dev residuals - spatial model")


plot_grid(p1, p2)
plot_grid(p1, p2, p3, p4, ncol = 2)
ggsave("local/plots/mapped_sp_vs_nsp.png", width = 300, height = 300, units = "mm")





