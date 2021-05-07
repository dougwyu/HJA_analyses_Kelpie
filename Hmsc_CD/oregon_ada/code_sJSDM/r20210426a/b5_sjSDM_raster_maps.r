
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
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada

library(dplyr)
library(raster)
library(sf)

utm10N <- 32610

gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"

resFolder <-"code_sjSDM/r20210426a/results"
abund <- "pa"


# load model data - for speices classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, minocc, noSteps, vars, varsName, abund, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv

spp <- data.frame(species = colnames(get(paste0("otu.", abund, ".csv")))) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID",
                                          "size"),
                  remove = FALSE, sep = "_") %>%
  dplyr::select(-empty)
head(spp)


#load(file.path(resFolder, "predictions.rdata")) # pred.sp, pred.mn, pred.sd, 
load(file.path(gis_out, "r_oversize/predictions.rdata")) # pred.mn, pred.sd, 

# load raster templates
load("data/gis/templateRaster.rdata") ## r, indNA

## load species AUC resutls for filtering
load(file.path(resFolder, "sp_results.rdata")) # sp.mn.test
rm(eval.results, sp.mn.train, sp.res.test, sp.res.train)

## Mean AUC per species (and other eval metrics)
str(sp.mn.test, max.level = 1)

## Filter species by auc
auc.filt <- 0.65
sum(sp.mn.test$auc > auc.filt)

## extract species over AUC filter
str(pred.sp, max.level = 1)

pred.in <- pred.sp[,sp.mn.test$auc > auc.filt]
dim(pred.in)

## make rasters

rList <- lapply(data.frame(pred.in), function(x) {
  
  tmp <- r
  tmp[indNA] <- x
  tmp
  
})

rStack <- stack(rList)

# threshold?
tr <- 0.5

rStack.bin <- raster::reclassify(rStack, rcl = c(0, tr, 0, tr, 1, 1))

rStack.sum <- sum(rStack)
names(rStack.sum) <- "sp sum"

spRich <- sum(rStack.bin)
names(spRich) <- "sp richness"
plot(stack(spRich, rStack.sum))

## Do maps of richness by groups.

#filter species first for those included by AUC
head(spp)
sapply(spp, function(x) sum(is.na(x)))

spp.auc <- spp[sp.mn.test$auc > auc.filt,]
sort(table(spp.auc$order, useNA = "always"))

top4 <- names(sort(table(spp.auc$order), decreasing = T)[1:4])
top4

plot(rStack[[spp.auc$order == "Lepidoptera"]])



spRich_order <- stackApply(rStack.bin, spp.auc$order, fun = sum)
names(spRich_order)
names(spRich_order) <- sub("index_", "", names(spRich_order))

plot(spRich_order[[top4]])

save(rStack, file = file.path(resFolder, "rasterStacks.rdata"))


library(ggplot2)

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

p <- ggplot()+
  geom_sf(data = hja.utm)+
  coord_sf(datum = sf::st_crs(32610))

str(p)
p$coordinates$limits

xlim <- c(559000, 572000)
ylim <- c(4894000, 4904000)

load(file.path(gis_out, "sample_sites.rdata"))

cut.sf <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut.utm <- st_transform(cut.sf, crs = utm10N)
plot(cut.utm)

cut.40 <- subset(cut.utm, YrsSinceDi > 40)

ggplot()+
  geom_sf(data = cut.40, bg = NA)+
  geom_sf(data = hja.utm)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)

hllshd <- raster(file.path(gis_out, "r_utm/bE_30m_hlshd.tif"))

## prepare raster data
allR <- stack(spRich_order, spRich, rStack.sum)
names(allR)
coords <- xyFromCell(allR, seq_len(ncell(allR)))
df1 <- as.data.frame(values(allR))
df1 <- cbind(coords, df1)
head(df1)

## rearrange
df2 <- df1 %>%
  tidyr::pivot_longer(cols = -c(x, y), names_to = "Order", values_to = "Sp_rich")

head(df2)

coords <- xyFromCell(hllshd, seq_len(ncell(hllshd)))
hsd.df <- data.frame(hsd = values(hllshd))
hsd.df <- cbind(coords, hsd.df)
head(hsd.df)

## Sp richness in landscape

p<- ggplot() + 
  geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225)), name = "sp richness") +
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
  coord_sf(datum = sf::st_crs(32610), ylim = c(4890000,4910000), xlim = c(550000,572000))

ggsave("../local/plots/spRich_map.png", p)  


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

ggsave("../local/plots/spRich_map_HJA.png", p)  



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

ggsave("../local/plots/spSum_map_HJA.png", p) 




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

ggsave("../local/plots/spRich_map_x_Order.png", p, width = 300, height = 300, units = "mm")



getwd()



