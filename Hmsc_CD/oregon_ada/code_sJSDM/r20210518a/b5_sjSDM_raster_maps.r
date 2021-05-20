
## Make raster from predictions


#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

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
library(raster)
library(sf)
library(ggplot2)

utm10N <- 32610

## on ADA
gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"

baseFolder <- "code_sjSDM/r20210518a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_results.rdata")) # sp.mn.test
rm(eval.results, sp.mn.train, sp.res.test, sp.res.train)

## Mean AUC per species (and other eval metrics)
str(sp.mn.test, max.level = 1)
head(sp.mn.test$auc)

## Filter species by auc
auc.filt <- 0.65
sum(sp.mn.test$auc > auc.filt)

## extract species over AUC filter
str(pred.sp, max.level = 1)

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
  mutate(auc = sp.mn.test$auc,
         incidence = incidence)

head(spp)

sum(is.na(spp$best.name))
sum(grepl("NA_NA", spp$best.name))
head(spp, 30)

sum(is.na(spp$family))






load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) # pred.mn, pred.sd, 

## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 

dim(pred.mn)

# load raster templates - reduced areas
load(file.path(gis_out, "templateRaster.rdata")) ## r, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting

pred.in <- pred.mn[,sp.mn.test$auc > auc.filt]
dim(pred.in)

## get species names too
spp.in <- spp[sp.mn.test$auc > auc.filt, ]
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

rStack <- stack(rList)
names(rStack) <- spp.in$best.name

## add auc incidence names to stack
names(rStack) <- paste0(spp.in$best.name, " ", "auc=", round(spp.in$auc, 2), " ",  "prev=", round(spp.in$incidence,2))


# threshold?
tr <- 0.5

rStack.bin <- raster::reclassify(rStack, rcl = c(0, tr, 0, tr, 1, 1))

rStack.sum <- sum(rStack)
names(rStack.sum) <- "sp sum"

spRich <- sum(rStack.bin)
names(spRich) <- "sp richness"
# plot(stack(spRich, rStack.sum))

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

##### Do plots #####

# function to add to each plot
addAll <- function(){
  
  plot(st_geometry(hja.utm), add = T, col = NA, border = "black")
  #plot(st_geometry(aoi.pred.sf), add = T, col = NA, border = "black")
  plot(st_geometry(xy.utm), add = T, col = "grey40", pch = 3, cex = 0.2)
  
}

plot(spRich)


#filter species first for those included by AUC
head(spp)
sapply(spp, function(x) sum(is.na(x)))

sort(table(spp.in$order, useNA = "always"))

top4 <- names(sort(table(spp.in$order), decreasing = T)[1:4])
top4

# plot(rStack[[which(spp.auc$order == "Lepidoptera")]])


pdf(file.path(plotsFolder, "coleoptera.pdf"), width = 7, height = 7)
plotStack(rStack[[which(spp.auc$order == "Coleoptera")]], addfun = addAll)
dev.off()


pdf(file.path(plotsFolder, "Hymenoptera.pdf"), width = 7, height = 7)
plotStack(rStack[[which(spp.auc$order == "Hymenoptera")]], addfun = addAll)
dev.off()

pdf(file.path(plotsFolder, "Lepidoptera.pdf"), width = 7, height = 7)
plotStack(rStack[[which(spp.auc$order == "Lepidoptera")]], addfun = addAll)
dev.off()

pdf(file.path(plotsFolder, "Diptera.pdf"), width = 7, height = 7)
plotStack(rStack[[which(spp.auc$order == "Diptera")]])
dev.off()


spRich_order <- stackApply(rStack.bin, spp.auc$order, fun = sum)
names(spRich_order)
names(spRich_order) <- sub("index_", "", names(spRich_order))

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

load(file.path(gis_out, "sample_sites.rdata"))

cut.sf <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut.utm <- st_transform(cut.sf, crs = utm10N)
# plot(cut.utm)

cut.40 <- subset(cut.utm, YrsSinceDi > 40)

# plot(r.aoi.pred)
# plot(hja.utm, add = T, col = NA)
# plot(xy.utm, add = T, pch = 16, col = "black")

# ggplot()+
#   geom_sf(data = cut.40, bg = NA)+
#   geom_sf(data = hja.utm)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

hllshd <- raster(file.path(gis_out, "r_utm/bE_30m_hlshd.tif"))


## prepare raster data
allR <- stack(spRich_order, spRich, rStack.sum)
names(allR)
coords <- xyFromCell(allR, seq_len(ncell(allR)))
df1 <- as.data.frame(values(allR))
df1 <- cbind(coords, df1)
# head(df1)

## rearrange
df2 <- df1 %>%
  tidyr::pivot_longer(cols = -c(x, y), names_to = "Order", values_to = "Sp_rich")

head(df2)

coords <- xyFromCell(hllshd, seq_len(ncell(hllshd)))
hsd.df <- data.frame(hsd = values(hllshd))
hsd.df <- cbind(coords, hsd.df)
head(hsd.df)

## Sp richness in landscape

# p<- ggplot() + 
#   geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
#   geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
#   coord_sf(datum = sf::st_crs(32610), ylim = c(4890000,4910000), xlim = c(550000,572000))
# 
# ggsave(file.path(plotsFolder, "spRich_map.pdf"), p)

## Sp richness in HJA
p<- ggplot() + 
  geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spRich_map_HJA.png"), p)


head(df1)
## Sp summed in HJA
p<- ggplot() + 
  geom_tile(data = df1, aes(x, y, fill = sp.sum)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "summed probability") +
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

ggsave(file.path(plotsFolder, "spSum_map_HJA.png"), p)


## Sp richness by order

p<- ggplot(subset(df2, Order %in% top4), aes(x, y)) + 
  geom_tile(aes(fill = Sp_rich)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20", inherit.aes = FALSE) +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  geom_sf(data = cut.40, bg = NA, inherit.aes = FALSE)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
  facet_wrap(~Order)

ggsave(file.path(plotsFolder, "spRich_map_x_Order.png"), p, width = 7, height = 7)


