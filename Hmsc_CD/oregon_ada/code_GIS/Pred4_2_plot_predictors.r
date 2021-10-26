

## Local
wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"
# wd <- "J:/UEA/gitHRepos/HJA_analyses_Kelpie"
setwd(wd)
# dir()
getwd()


library(raster)
library(sf)

utm10N <- 32610
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"

### 1. Load data #### 

# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))
# get names and name groups
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck
names(allBrck)

## Load predictor names and codes and rename
tableS1 <- read.csv("Hmsc_CD/oregon_ada/code_GIS/table_S_predictors_20210820.csv")
head(tableS1)

sum(tableS1$PredictorCode %in% allNames)
tableS1$PredictorCode[!tableS1$PredictorCode %in% allNames]
#  "lg_cover2m_4m"  "lg_cover2m_max" "lg_cover4m_16m" "lg_DistRoad"    "lg_DistStream"
# "l_Cover_2m_4m"   "l_Cover_2m_max"    "l_Cover_4m_16m","DistRoad","DistStream"

## get only 58 candidate predictors
allBrck <- subset(allBrck, tableS1$PredictorCode[tableS1$PredictorCode %in% allNames])
names(allBrck)

# subset names file and reorder
tableS1 <- subset(tableS1, PredictorCode %in% names(allBrck))
unique(tableS1$Group)

tableS1$Group <- factor(tableS1$Group, 
                        levels = c("Lidar - Canopy","Topography","Landsat - annual","Landsat - bands","Anthropogenic"))
tableS1 <- tableS1[order(tableS1$Group, tableS1$PredictorName), ]
tableS1[,c(1:2,4)]

## reorder raster stack
allBrck <- subset(allBrck, tableS1$PredictorCode)
names(allBrck)

## Use these names for plot (check these symbols work as names on raster... ... NO!)
plotNames <- tableS1$PredictorName[match(names(allBrck), tableS1$PredictorCode)]



## templates 
load("Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata") # aoi.pred.sf
rm(r.msk, indNA, r.aoi.pred)

# get edited prediction area
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))

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
  plot(st_geometry(aoi.pred.sf_edit), add = T, col = NA, border = "black")
  plot(st_geometry(xy.utm), add = T, col = "black", pch = 3, cex = 0.5)
  
}

pdf("Hmsc_CD/local/plots/predictors_trial.pdf", width = 7, height = 10)
plot(allBrck[[1:6]], nc = 2, nr = 3, addfun = addAll, main = plotNames[1:6])
dev.off()

66%/%12

by <- 6
end <- 1:(nlayers(allBrck)%/%by)*by
start <- (end - by) + 1

start
end

#i = 1

pdf("Hmsc_CD/local/plots/predictors_all_plots.pdf", width = 7, height = 10)

for(i in seq_along(start)){
  
  plot(allBrck[[ start[i]:end[i] ]], nc = 2, nr = 3, addfun = addAll, main= plotNames[ start[i]:end[i] ])
}

dev.off()

### some individual plots

n <- 36

extent(allBrck[[i]])
hcl.pals()

png("Hmsc_CD/local/plots/predictors_small%02d.png", width = 500, height = 500)
par(bg = NA, mar = c(0,0,0,0), oma = c(0,0,0,0))
for(i in 1:n) image(allBrck[[i]], axes = F,  
                    col = hcl.colors(12, c("YlOrRd",
                                           "Viridis",
                                           "Plasma",
                                           "Inferno",
                                           "Terrain",
                                           "Terrain 2")[i%%6+1], rev = TRUE))
dev.off()

