
## Make raster from predictions


#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

## trial run 
options(echo=TRUE) # if you want see commands in output file
# Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
# library(sjSDM)
# packageVersion("sjSDM")
# getwd() # always run sub from oregon_ada

library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(ggplot2)

utm10N <- 32610

## on ADA
gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210716b"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_test_results.rdata")) # eval.results, sp.res.test, sp.res.train

## Mean AUC per species (and other eval metrics)
str(sp.res.test, max.level = 1)
head(sp.res.test$auc)

sum(is.na(sp.res.test$auc))

## Filter species by auc
auc.filt <- 0.70
sum(sp.res.test$auc > auc.filt, na.rm = T)

# ## extract species over AUC filter
# str(pred.sp, max.level = 1)

# incidence 
incidence <- colSums(otu.pa.csv)/nrow(otu.pa.csv)


spp <- data.frame(species = colnames(get(paste0("otu.", abund, ".csv")))) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID",
                                          "size"),
                  remove = FALSE, sep = "_", convert = TRUE) %>%  ## creates real NAs with convert = T
  mutate(best.name = case_when(is.na(epithet) & is.na(genus) & is.na(family) & is.na(order) ~ class,
                               is.na(epithet) & is.na(genus) & is.na(family) ~ order,
                               is.na(epithet) & is.na(genus) ~ family,
                               is.na(epithet) ~ genus,
                               TRUE ~ paste(genus, epithet, sep = "_")
                               )) %>%
  dplyr::select(-empty)%>%
  mutate(auc = sp.res.test$auc,
         incidence = incidence,
         best.name = paste(best.name, BOLDID, sep = "_"))

head(spp)

sum(is.na(spp$best.name))
sum(grepl("NA_NA", spp$best.name))
head(spp, 30)

sum(is.na(spp$family))

## In ADA

# extrapolated predictions
load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")))
# pred.mn, pred.sd, 

# clamp predictions
load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl



## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 
# paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")

dim(pred.mn)

# load raster templates - reduced areas
load(file.path(gis_out, "templateRaster.rdata")) ## r.msk, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting
# plot(r.msk)
# plot(aoi.pred.sf)

### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "shape/aoi_pred_sf_edit.shp"))
aoi.pred.sf_edit <- st_make_valid(aoi.pred.sf_edit)

pred.in <- pred.mn[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
dim(pred.in)

# clamp predictions filtered by species
pred.in.cl <- pred.mn.cl[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]

## get species names too
spp.in <- spp[sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc), ]
head(spp.in)

## make rasters
# plot(r.aoi.pred)
# x <- data.frame(pred.in)[,1]

rList <- lapply(data.frame(pred.in), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})
# plot(tmp)
rStack.ex <- stack(rList)
names(rStack.ex) <- spp.in$best.name
rStack.ex

## add auc incidence names to stack
names(rStack.ex) <- paste0(spp.in$best.name, " ", "auc=", round(spp.in$auc, 2), " ",  "prev=", round(spp.in$incidence,2))
rm(rList)

## clamp version
rList <- lapply(data.frame(pred.in.cl), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})
# plot(tmp)
rStack.cl <- stack(rList)
names(rStack.cl) <- spp.in$best.name
rStack.cl

## add auc incidence names to stack
names(rStack.cl) <- paste0(spp.in$best.name, " ", "auc=", round(spp.in$auc, 2), " ",  "prev=", round(spp.in$incidence,2))


# threshold?
tr <- 0.5

rStack.bin.ex <- raster::reclassify(rStack.ex, rcl = c(0, tr, 0, tr, 1, 1))
rStack.bin.cl <- raster::reclassify(rStack.cl, rcl = c(0, tr, 0, tr, 1, 1))

rStack.sum.ex <- sum(rStack.ex)
names(rStack.sum.ex) <- "sp sum"

rStack.sum.cl <- sum(rStack.cl)
names(rStack.sum.cl) <- "sp sum"

spRich.ex <- sum(rStack.bin.ex)
names(spRich.ex) <- "sp richness"
spRich.cl <- sum(rStack.bin.cl)
names(spRich.cl) <- "sp richness"
# plot(stack(spRich.ex, rStack.sum.ex))

## Do maps of richness by groups.

## load extra layers for plotting
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)


## Make species richness stack
rStack.bin.cl

spRich_order.cl <- stackApply(rStack.bin.cl, spp.in$order, fun = sum)
names(spRich_order.cl)
names(spRich_order.cl) <- sub("index_", "", names(spRich_order.cl))

spRich_order.ex <- stackApply(rStack.bin.ex, spp.in$order, fun = sum)
names(spRich_order.ex)
names(spRich_order.ex) <- sub("index_", "", names(spRich_order.ex))

# 
# # ## Load ordination results
# # load(file.path(resFolder, "ord_tsne_res.rdata")) #
# # #
# # stck <- raster::stack(rSites1, rSites2)
# # names(stck) <- c("tsne1", "tsne2")
# 
# ## save rasters
# 
# writeRaster(spRich_order.cl, bylayer = T, 
#             filename = file.path(resFolder, "spRich_.tif"), suffix = "names", overwrite = TRUE, datatype = "FLT4S")

writeRaster(rStack.sum.ex, filename = file.path(resFolder, "spSum_ex.tif"), datatype = "FLT4S", overwrite = T)
writeRaster(spRich.ex, filename = file.path(resFolder, "spRich_all_ex.tif"), datatype = "FLT4S", overwrite = T)

writeRaster(rStack.sum.cl, filename = file.path(resFolder, "spSum_cl.tif"), datatype = "FLT4S", overwrite = T)
writeRaster(spRich.cl, filename = file.path(resFolder, "spRich_all_cl.tif"), datatype = "FLT4S", overwrite = T)

# writeRaster(stck, bylayer = T, 
#             filename = file.path(resFolder, "ord.tif"), suffix = "names", overwrite = TRUE, datatype = "FLT4S")

sppFolder <- file.path(resFolder, "spp_tifs_cl")
if(!dir.exists(sppFolder)) dir.create(sppFolder)
writeRaster(rStack.cl, bylayer = T, filename = file.path(sppFolder, "spp_cl.tif"), suffix = "names", datatype = "FLT4S", overwrite = T)




##### Do plots #####

# function to add to each plot
addAll <- function(){

  plot(st_geometry(hja.utm), add = T, col = NA, border = "black")
  #plot(st_geometry(aoi.pred.sf), add = T, col = NA, border = "black")
  plot(st_geometry(xy.utm), add = T, col = "grey40", pch = 3, cex = 0.2)

}

# plot(spRich)


#filter species first for those included by AUC
head(spp)
sapply(spp, function(x) sum(is.na(x)))

sort(table(spp.in$order, useNA = "always"))

top4 <- names(sort(table(spp.in$order), decreasing = T)[1:4])
top4

# plot(rStack[[which(spp.auc$order == "Lepidoptera")]])

source("code_GIS/plotStack.r")

# pdf(file.path(plotsFolder, "coleoptera.pdf"), width = 7, height = 7)
# plotStack(rStack[[which(spp.in$order == "Coleoptera")]], addfun = addAll)
# dev.off()
# 
# pdf(file.path(plotsFolder, "Hymenoptera.pdf"), width = 7, height = 7)
# plotStack(rStack[[which(spp.in$order == "Hymenoptera")]], addfun = addAll)
# dev.off()
# 
# pdf(file.path(plotsFolder, "Lepidoptera.pdf"), width = 7, height = 7)
# plotStack(rStack[[which(spp.in$order == "Lepidoptera")]], addfun = addAll)
# dev.off()
# 
# pdf(file.path(plotsFolder, "Diptera.pdf"), width = 7, height = 7)
# plotStack(rStack[[which(spp.in$order == "Diptera")]])
# dev.off()
# 
# rStack.bin
# 
# spRich_order <- stackApply(rStack.bin, spp.in$order, fun = sum)
# names(spRich_order)
# names(spRich_order) <- sub("index_", "", names(spRich_order))
# 
# 
# pdf(file.path(plotsFolder, "Sp_rich_order.pdf"), width = 7, height = 7)
# plotStack(spRich_order, by = 4, nc = 2, nr= 2)
# dev.off()


# plot(spRich_order[[top4]])

# save(rStack, file = file.path(resFolder, "rasterStacks.rdata"))


# p <- ggplot()+
#   geom_sf(data = hja.utm)+
#   coord_sf(datum = sf::st_crs(32610))
# str(p)
# p$coordinates$limits

st_bbox(aoi.pred.sf)
xlim <- c(553500, 573500)
ylim <- c(4889000, 4910000)


st_bbox(aoi.pred.sf_edit)
xlim <- c(554650, 572250)
ylim <- c(4890780, 4908610)

cut.sf <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut.utm <- st_transform(cut.sf, crs = utm10N)
# # plot(cut.utm)
#
# ## Less than 40 years since disturbance
cut.40 <- subset(cut.utm, YrsSinceDi < 40)


# ## clip to prediction area
cut.intsct <- st_intersects(cut.40, aoi.pred.sf_edit)
cut.40 <- cut.40[lengths(cut.intsct)>0,]

#plot(st_geometry(aoi.pred.sf_edit))
#plot(cut.40, add =T)

# clip to prediction area
tmp <- st_difference(aoi.pred.sf_edit, sf::st_union(cut.40))
cut.40 <- st_difference(aoi.pred.sf_edit, tmp)

#plot(st_geometry(aoi.pred.sf_edit))
#plot(cut.40, add = T)


# plot(r.aoi.pred)
# plot(hja.utm, add = T, col = NA)
# plot(xy.utm, add = T, pch = 16, col = "black")

# ggplot()+
#   geom_sf(data = cut.40, bg = NA)+
#   geom_sf(data = hja.utm)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)
#
hllshd <- raster(file.path(gis_out, "r_utm/bE_30m_hlshd.tif"))


# mask to species predictrions area
aoi.pred.sf_edit

#plot(hllshd)
#plot(spRich)

tmp.msk <- rasterize(aoi.pred.sf_edit, hllshd)
#plot(tmp.msk)

hllshd <- mask(hllshd, tmp.msk)
hllshd <- crop(hllshd, spRich.ex)

# ## prepare raster data
allR <- stack(spRich_order.cl, spRich.cl, rStack.sum.cl)
names(allR)
coords <- xyFromCell(allR, seq_len(ncell(allR)))
df1 <- as.data.frame(values(allR))
df1 <- cbind(coords, df1)
# head(df1)

## rearrange
df2 <- df1 %>%
  tidyr::pivot_longer(cols = -c(x, y), names_to = "Order", values_to = "Sp_rich")

head(df2)

## extrapolated
allR.ex <- stack(spRich_order.ex, spRich.ex, rStack.sum.ex)
names(allR.ex)
coords <- xyFromCell(allR.ex, seq_len(ncell(allR.ex)))
df3 <- as.data.frame(values(allR.ex))
df3 <- cbind(coords, df3)
# head(df1)

## rearrange
df4 <- df3 %>%
  tidyr::pivot_longer(cols = -c(x, y), names_to = "Order", values_to = "Sp_rich")

head(df4)

### hillshade raster - prepare for ggplot
coords <- xyFromCell(hllshd, seq_len(ncell(hllshd)))
hsd.df <- data.frame(hsd = values(hllshd))
hsd.df <- cbind(coords, hsd.df)
head(hsd.df)

## save plotting data

save(hsd.df, df1, df2, df3, df4, cut.40, hja.utm, xy.utm, spRich.ex, spRich.cl,
     rStack.sum.ex, rStack.sum.cl,rStack.bin.ex, rStack.bin.cl, spp.in, spp, xlim, ylim,
     file = file.path(resFolder, "plotData.rdata"))



# ## Sp richness in landscape
#
# # p<- ggplot() +
# #   geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
# #   scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
# #   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
# #   scale_alpha(range = c(0.25, 0.65), guide = "none")+
# #   geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
# #   geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
# #   coord_sf(datum = sf::st_crs(32610), ylim = c(4890000,4910000), xlim = c(550000,572000))
# #
# # ggsave(file.path(plotsFolder, "spRich_map.pdf"), p)




## Sp richness in HJA
p1 <- ggplot() +
  geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
  #geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spRich_map_HJA.pdf"), p1)



head(df1)
## Sp summed in HJA
p2 <- ggplot() +
  geom_tile(data = df1, aes(x, y, fill = sp.sum)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "summed probability") +
  #geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spSum_map_HJA.pdf"), p2)


### clamp version

p3 <- ggplot() +
  geom_tile(data = df3, aes(x, y, fill = sp.richness)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness - cl") +
  #geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spRich_map_HJA_cl.pdf"), p1)



head(df1)
## Sp summed in HJA
p4 <- ggplot() +
  geom_tile(data = df3, aes(x, y, fill = sp.sum)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "summed probability - cl") +
  #geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  #scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spSum_map_HJA_cl.pdf"), p2)


# library(cowplot)
# 
# p5 <- plot_grid(p1, p3, p2, p4, nrow = 2, ncol= 2, labels = c("richness - no clamping", "clamping", "summed - no clamping", "clamping") )
# ggsave(file.path(plotsFolder, "clamp_comparision.pdf"), p5)



## Sp richness by order
# 
# p3 <- ggplot(subset(df2, Order %in% top4), aes(x, y)) +
#   geom_tile(aes(fill = Sp_rich)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20", inherit.aes = FALSE) +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
#   geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
#   geom_sf(data = cut.40, bg = NA, inherit.aes = FALSE)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   facet_wrap(~Order)
# 
# ggsave(file.path(plotsFolder, "spRich_map_x_Order.pdf"), p3, width = 7, height = 7)
#
#
# ### Figure plots
#
#
# st_bbox(aoi.pred.sf_edit)
xlim <- c(554650, 572250)
ylim <- c(4890780, 4908610)

head(df2)
# ## Prepare ordination data
#
# # ## Load ordination results
# load(file.path(resFolder, "ord_tsne_res.rdata")) #
# #
# # ## function to load site scores into raster
# # # makeR <- function(r, siteScores, NAs) {
# # #
# # #   rSites <- raster(r)
# # #   rSites[] <- NA
# # #   rSites[NAs] <- siteScores
# # #   rSites
# # #
# # # }
# # #
# # # rSites1 <- makeR(r, tsne$Y[,1], NAs)
# # # rSites2 <- makeR(r, tsne$Y[,2], NAs)
# #
# stck <- raster::stack(rSites1, rSites2)
# names(stck) <- c("tsne1", "tsne2")


# 
# # prepare for ggplot
# coords <- xyFromCell(stck, seq_len(ncell(stck)))
# tsne.df <- data.frame(values(stck))
# tsne.df <- cbind(coords, tsne.df)
# head(tsne.df)
# 
# unique(df1$sp.richness)
# 
# #
# myTheme <- theme(panel.grid.major = element_line(color = '#cccccc'
#                                                  ,linetype = 'dashed'
#                                                  ,size = .3),
#                  axis.title = element_blank(),
#                  legend.key = element_rect(fill = "grey80"))
# 
# myTheme2 <- theme(panel.grid.major = element_line(color = '#cccccc'
#                                                   ,linetype = 'dashed'
#                                                   ,size = .3),
#                   axis.title = element_blank(),
#                   legend.key = element_rect(fill = "grey80"),
#                   axis.text = element_blank())
# 
# 
# 
# ## species richness
# p1 <- ggplot() +
#   geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
#   scale_fill_fermenter(palette = "YlGnBu", type = "seq", name = "Species richness")+ ## binned scale
#   # scale_fill_gradientn(colours = heat.colors(225), name = "Species richness") +
#   # scale_fill_viridis_c(name = "Species richness", option = "C", alpha = 0.7)+
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, aes(color = "A"), bg = NA, inherit.aes = FALSE, show.legend = "polygon")+
#   geom_sf(data= xy.utm, aes(color = "B"), size = 2.5, inherit.aes = FALSE, show.legend = "point")+
#   geom_sf(data = cut.40, alpha = 0.3, bg = NA, aes(color = "C"), show.legend = "polygon")+
#   scale_color_manual(values = c("A" = "blue", "B"= "white", "C" = "grey50"),
#                      labels = c("HJA border", "Sample sites", "Logged within 40 years"),
#                      name = NULL,
#                      guide = guide_legend(override.aes = list(linetype = c("solid", "blank", "solid"),
#                                                               shape = c(NA, 16, NA))))+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme
# 
# #unlink(file.path(plotsFolder, "results_p1.png"))
# ggsave(file.path(plotsFolder, "results_p1.pdf"), p1, width = 300, height = 300, units= "mm")
# 
# ## ordination
# p2 <- ggplot() +
#   geom_tile(data = tsne.df, aes(x, y, fill = tsne1)) +
#   scale_fill_viridis_c(name = "T-SNE axis 1", option = "B", alpha = 0.7)+
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
#   geom_sf(data= xy.utm, col = "white", size = 2.5, inherit.aes = FALSE)+
#   geom_sf(data = cut.40, col = "grey50", bg = NA)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme
# 
# #unlink(file.path(plotsFolder, "results_tsne1.png"))
# ggsave(file.path(plotsFolder, "results_tsne1.pdf"), p2, width = 300, height = 300, units= "mm")
# 
# ## ordination
# p2 <- ggplot() +
#   geom_tile(data = tsne.df, aes(x, y, fill = tsne2)) +
#   scale_fill_viridis_c(name = "T-SNE axis 2", option = "B", alpha = 0.7)+
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
#   geom_sf(data= xy.utm, col = "white", size = 2.5, inherit.aes = FALSE)+
#   geom_sf(data = cut.40, col = "grey50", bg = NA)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme
# 
# #unlink(file.path(plotsFolder, "results_tsne2.png"))
# ggsave(file.path(plotsFolder, "results_tsne2.pdf"), p2, width = 300, height = 300, units= "mm")
# 
# 

#
