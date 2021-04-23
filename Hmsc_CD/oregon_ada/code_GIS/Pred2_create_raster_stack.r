
## Local
# wd <- here::here()
# wd
# setwd(wd)
# dir()
getwd()

library(raster)
library(sf)

gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"

## Eventual data files: 
# gis_in <- "HJA_scripts/10_eo_data/raw_gis_data"
# gis_out <- "HJA_scripts/10_eo_data/processed_gis_data"

### 1. Load processed data ####

## Add all rasters to a stack, export as single multilayered tif and import as rasterbrick.

# Could also get all tifs (list.files() and then stack(filenames)

## 1 Canopy
load(file.path(gis_out, "be_ht.rdata"))
## 2 Cut
load(file.path(gis_out, "cut_stack.rdata"))
## 3. Topography
load(file.path(gis_out, "terr30.rdata"))
## 4. Lidar
load(file.path(gis_out, "lidarStack.rdata"))
## 5. Streams/Roads
load(file.path(gis_out, "admStck.rdata"))
## 6. Annual indices
load(file.path(gis_out, "annualStack.rdata"))
## 7. Temperature
load(file.path(gis_out, "tp.rdata"))

## Stack up 
allStck <- stack(be_ht, cutStack, terr30, lidarStck, admStack, annualStack, tp_r30)
allNames <- names(allStck)
allNames

# get name groups
nameList <- list(be_ht = be_ht.names, cutStack = cut.names, terr30 = terr30.names, 
                  lidarStck = lid.names, admStack = adm.names, annualStack = annual.names, tp_r30 = tp.names)

writeRaster(allStck, bylayer = F, filename = file.path(gis_out, "r_utm/allStack.tif"), overwrite = TRUE)
save(allNames, nameList, file = file.path(gis_out, "allNames.rdata"))

rm(be_ht, cutStack, terr30, lidarStck, admStack, annualStack, tp_r30,
   adm.names, annual.names, be_ht.names, cut.names, lid.names, terr30.names, tp.names, all.names)

## 2. Extract point values ####

# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck

rm(allNames, nameList)

## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

allVars <- data.frame(SiteName = xy.utm$SiteName, extract(allBrck, xy.utm))
head(allVars)

allVars$insideHJA <- ifelse(allVars$insideHJA==0, "no", "yes")
allVars$insideHJA <- factor(allVars$insideHJA, levels = c("no", "yes"))
table(allVars$insideHJA)

str(allVars)

summary(allVars)

save(allVars, file = file.path(gis_out, "envVars.rdata"))
save(allVars, file = file.path("Hmsc_CD/oregon_ada/data", "envVars.rdata"))

# load(file.path(gis_out, "envVars.rdata"))

