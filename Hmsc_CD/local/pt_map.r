### MAke map of points... 

getwd()
wd <- here::here()
setwd(wd)
wd
dir()

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")
source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/GIS/adjExt.r")

## Updated to new vars, also changes to elevation_m, canopy_height_m  to _f. 

library(dplyr)
library(sf)

# wgs84 UTM 10N
utm10N <- 32610
# EPSG:26910  NAD83 / UTM zone 10N
nadutm10 <- 26910
# EPSG:4269 # NAD 83
# nad83 <- 4269


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

fn

# what file am i using?
basename(fn)

# when was it modified?
file.mtime(fn)

# read complete data set - has coordinates, names, all spp and all predictors
otuenv <- read.csv(fn, stringsAsFactors = FALSE, na.strings = "NA")

# clean up
rm(datFile, gitHub, kelpierundate, minimaprundate, outputidxstatstabulatefolder, primer, samtoolsfilter, samtoolsqual, fn)

#### Explore species at traps M1 and M2 #####

## check S1, S2, M1 and M2 distribution of sample locations
pts.chk <- otuenv %>%
  select(c(UTM_N, UTM_E, SiteName, period, trap))

pts.chk

# Check all sitenames have unique coordinates:  
uniqueSites <- unique(pts.chk[, c("UTM_N", "UTM_E", "SiteName")])
length(uniqueSites$SiteName)

## check what periods, traps at each site
pts.chk2 <- otuenv %>%
  filter(period == "S1") %>%
  select(c(UTM_N, UTM_E, SiteName, trap)) %>%
  tidyr::pivot_wider(names_from = trap, values_from = trap, values_fn = length)
pts.chk2
# pts.chk2 <- reshape2::dcast(subset(pts.chk, period == "S1"), SiteName ~ trap, fun.aggregate = length)

data.frame(pts.chk2)
colSums(pts.chk2[,c("M1", "M2")], na.rm = T)

# paired sites
data.frame(pts.chk2[which(pts.chk2$M1 == 1 & pts.chk2$M2 == 1),])
nrow(pts.chk2[which(pts.chk2$M1 == 1 & pts.chk2$M2 == 1),]) # 32 sites
data.frame(pts.chk2[which(pts.chk2$M1 == 1 & pts.chk2$M2 == 1),])$SiteName

## m2 not at m1?
data.frame(pts.chk2[which(is.na(pts.chk2$M1) & pts.chk2$M2 == 1),])
# 1 4897312 566184 HOBO-036 NA  1


# species numbers - prevalence and overlap between M1 and M2
# keep OTUs with >=5 incidences
minocc <- 6 # set to high number (e.g. 20) for testing
otu.qp.csv <- otuenv %>% 
  #dplyr::filter(period == "S1" & trap == "M1")%>%
  dplyr::filter(period == "S1")%>%
  dplyr::select(contains("__")) ## file above is already qp

dim(otu.qp.csv)
#rowSums(otu.qp.csv>0)

otu.qp.csv <- otu.qp.csv[ , colSums(otu.qp.csv > 0) >= minocc] ## if using both M1 and M2, then this gives minocc across all samples
otus <- colnames(otu.qp.csv)
length(unique(otus))


head(colSums(otu.qp.csv > 0))
max(colSums(otu.qp.csv > 0))

# presence absence  per site per trap
sp.chk <- otuenv %>% 
  dplyr::filter(period == "S1")%>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0)

# otus * m1 sites + otus * m2 sites
1210*88 + 1210*33
sp.chk

# All OTUs, and presence at M1 or M2
sp.chk <- otuenv %>% 
  dplyr::filter(period == "S1")%>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T))
# all OTUs
1210*2
sp.chk

## OTUs > minocc at M1 and M2, but filtering for minocc by M1 and M2 separately.
sp.chk <- otuenv %>% 
  dplyr::filter(period == "S1")%>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc)
## species >= minocc
sp.chk # nrow is 304 (222 + 82)

# Number of Otus per site
sum(sp.chk$trap == "M1") # 222
sum(sp.chk$trap == "M2") # 82

otuenv %>% 
  dplyr::filter(period == "S1")%>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>%
  group_by(trap) %>%
  summarise(nOTU = n())


## Species by M1 and M2, with minocc calculated per trap
sp.chk <- otuenv %>% 
  dplyr::filter(period == "S1")%>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0)

sp.chk

nrow(sp.chk)
# total OTUs is 225

# Shared species
otuenv %>% 
  dplyr::filter(period == "S1") %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0) %>%
  summarise(Total = n(),
            M1_total = sum(M1, na.rm =T),
            M2_total = sum(M2, na.rm = T),
            shared = sum(M1 & M2, na.rm = T),
            M1_only = sum(M1[is.na(M2)]),
            M2_only = sum(M2[is.na(M1)]))


# Total M1_total M2_total shared M1_only M2_only
# <int>    <int>    <int>  <int>   <int>   <int>
#  225      222       82     79     143       3
## get predictors and coordiantes, sitenames


## Make species richness per sites per trap - and add to spatial 
sp.rch <- otuenv %>% 
  dplyr::filter(period == "S1") %>%
  dplyr::select(SiteName, period, trap, UTM_N, UTM_E, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(UTM_N, UTM_E, SiteName, period, trap) %>%
  summarise(spRich = sum(value, na.rm = T)) %>%
  tidyr::pivot_wider(names_from = c(period, trap), values_from = spRich, names_sort = TRUE) %>%
  mutate(Trap = case_when(!is.na(S1_M1) & is.na(S1_M2) ~ "M1",
                            is.na(S1_M1) & !is.na(S1_M2) ~ "M2",
                            !is.na(S1_M1) & !is.na(S1_M2) ~ "M1M2",
                              TRUE ~ "CHECK"))
  
  # mutate(Period = case_when(all(!is.na(S1_M1) | !is.na(S1_M2), !is.na(S2_M1) | !is.na(S2_M2)) ~ "S1S2",
  #                           all(is.na(S1_M1), is.na(S1_M2)) ~ "S2",
  #                           all(is.na(S2_M1), is.na(S2_M2)) ~ "S1",
  #                           TRUE ~ "CHECK"))

sp.rch
data.frame(sp.rch)

# max richness
sp.rch %>%
  ungroup() %>%
  summarise(across(contains(c("S1")), max, na.rm = T))


# Make spatial
rch.sf <- st_as_sf(sp.rch, coords = c("UTM_E", "UTM_N"), crs = nadutm10)

# export all points to shapefile / kml
st_write(rch.sf, "J:/UEA/Oregon/gis/s_nad_utm/sp_rich_pts.shp", delete_layer = T)

# transform to wgs and export to kml
rch.wgs <- st_transform(rch.sf, crs = 4326)
st_write(rch.wgs, "J:/UEA/Oregon/gis/s_wgs/sp_rich_all_pts.kml", delete_layer =T)

#transform to UTM
rch.utm <- st_transform(rch.sf, crs = utm10N)
st_write(rch.utm, "J:/UEA/Oregon/gis/s_utm/sp_rich_S1_m1m2.shp", delete_layer =T)

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

pts.sf <- st_as_sf(env.vars, coords = c("UTM_E", "UTM_N"), crs = nadutm10)

pts.sf

# export all points to shapefile / kml
# st_write(pts.sf, "J:/UEA/Oregon/gis/s_nad_utm/all_pts.shp")

# transform to wgs and export to kml
pts.wgs <- st_transform(pts.sf, crs = 4326)
# st_write(pts.wgs, "J:/UEA/Oregon/gis/s_wgs/all_pts.kml")

## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data

hja <- st_read("J:/UEA/Oregon/gis/s_nad_utm/HJA_Boundary.shp")

hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja_bound

hja.wgs <- st_transform(hja_bound, crs = 4326)
# st_write(hja.wgs, "J:/UEA/Oregon/gis/s_wgs/hja.kml")
# st_write(hja.wgs, "J:/UEA/Oregon/gis/s_wgs/hja.shp")
# 
# hja_utm <- st_transform(hja.wgs, crs = utm10N)
# st_write(hja_utm, "J:/UEA/Oregon/gis/s_utm/hja.shp")

## Make a bounding box for study area on GEE
plot(st_geometry(pts.sf))

ext.std<- st_bbox(pts.sf)
str(bb)

bb.sf <- adjExt(st_buffer(pts.sf, 5000), d = 1000)
plot(st_geometry(bb.sf), col = NA)
plot(st_geometry(pts.sf), add = T)

st_write(bb.sf,  "J:/UEA/Oregon/gis/s_wgs/bb_pts.shp")
st_write(bb.sf,  "J:/UEA/Oregon/gis/s_wgs/bb_pts.kml")

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

rch.vars <- c("S1_M1","S1_M2","Trap")
mv3 <- mapview::mapview(rch.sf, zcol = rch.vars, homebutton = FALSE,
                        map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT"),
                        legend = T, alpha.regions = 1)
mv3@map <- mv3@map %>% leaflet::hideGroup(group = rch.vars[-1])

mv3 + mv2

mv <- mv1 + mv2
mv

# save
mapview::mapshot(mv, url = "Hmsc_CD/local/plots/pts_map.html", selfcontainted = TRUE)
# only works in wd at moment. moved from there.
mapview::mapshot(mv, url = "pts_map.html", selfcontainted = TRUE)
getwd()

## Species richness interpolated raster
## inverse distance weighted (IDW)

#### gstat examples
library(gstat)
library(raster)

rch.utm <- st_transform(rch.sf, crs = utm10N)
rch.ch <- st_buffer(st_convex_hull(st_union(rch.utm)),2500)

rch <- data.frame(x = st_coordinates(rch.utm)[,1],
                  y = st_coordinates(rch.utm)[,2],
                            rch.utm)
rch$geometry <- NULL
head(rch)



# get template raster
r <- raster("HJA_scripts/10_eo_data/raw_gis_data/r_utm/lidar_metric_mosaic_Cover_2m_4m.tif")
r <- raster(r)
r
res(r) <- c(500,500)

# set up model
idwS1 <- gstat(id = "S1_M1", formula = S1_M1~1, locations = ~x+y, 
              data = rch[complete.cases(rch$S1_M1),], nmax=7, set=list(idp = 0.5))
idwS1

s1 <- interpolate(r, idwS1)

s1 <- mask(s1, as(rch.ch, "Spatial"))
plot(s1, main = "S1_M1 sp richness")
#plot(rch.utm[,"S1_M1"], pch = 16, add = T, pal = rev(terrain.colors(10)))
plot(rch.utm[,"S1_M1"], pch = 16, add = T, pal = rev(heat.colors(10)))
plot(hja_bound, add = T, col = NA, border = "blue")

png("Hmsc_CD/local/plots/sp_rich_interpolation.png", height = 200, width = 200, units = "mm", res = 100)
plot(s1, main = "S1_M1 sp richness")
#plot(rch.utm[,"S1_M1"], pch = 16, add = T, pal = rev(heat.colors(10)))
plot(hja_bound, add = T, col = NA, border = "blue")
dev.off()


load("Hmsc_CD/oregon_ada/results/ecocopula/ecocopula_modlist_pilot.rdata")
# modList
# mod 2 ia otu.pa ~ be500 * oldGrowthIndex + tri + insideHJA + lg_YrsDisturb)
site_res <- modList[[2]]$site[,c("Factor1", "Factor2", "UTM_E", "UTM_N")]
head(site_res)
head(rch)

eco_rch <- merge(site_res, data.frame(rch.sf, st_coordinates(rch.sf)), by.x = c("UTM_E", "UTM_N"), by.y = c("X", "Y"))
head(eco_rch)
plot(eco_rch$Factor1, eco_rch$S1_M1)
cor.test(eco_rch$Factor1, eco_rch$S1_M1, method = "spearman")

tmp <- merge(rch.utm, eco_rch[,c("SiteName", "Factor1", "Factor2")], by = "SiteName")
plot(tmp$S1_M1 ~ tmp$Factor1)


plot(s1)
# plot(rch.utm[,"S1_M1"], pch = 16, add = T, pal = rev(terrain.colors(10)))
plot(hja_bound, add = T, col = NA, border = "blue")
plot(tmp[,"Factor1"], add =T, pch = 16, pal = heat.colors(10))



# alpha = 1 - no transparency... 
symbols(1,1, circle = 1, bg = rgb(0,0,0,1))
symbols(1.1,1.1, circle = 1, bg = rgb(0,0,0,1), add = T)

symbols(1,1, circle = 1, bg = rgb(0,0,0,0.5))
symbols(1.1,1.1, circle = 1, bg = rgb(0,0,0,0.5), add = T)

## get US states for ref map
usa1.wgs <- raster::getData("GADM", country= "USA", path= "J:/UEA/Oregon/gis/s_wgs", level = 1)

usa1.sf <- st_as_sf(usa1.wgs)
st_write(usa1.sf, dsn = "J:/UEA/Oregon/gis/s_wgs", layer = "gadm_usa.shp", driver = "ESRI Shapefile")
