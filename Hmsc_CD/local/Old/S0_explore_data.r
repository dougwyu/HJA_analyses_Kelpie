## S1 Data read from github ####

## Christian D
## 20/11/2020  - start date

# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"  - 
# if you send the .sub from here. whereever you send the sub from is wd
# with folders Oregon, etc... 

# tilde expansion relative to HOME system env variable ## on ADA, HOME is "/gpfs/home/hsp20azu"
setwd("~/Oregon_winscp")  # ./ relative to working directory 
dir()

# scaled data ready to use.. for both pa and qp (eg for JSDM) here:
# https://github.com/dougwyu/HJA_analyses_Kelpie/tree/master/sjSDM/R-git/data/kelpie_data/for_adagpu/data_20201119_5minocc_gismslidar

# ALL Predictors
# unscaled data, corrected for FSL, divided into Sessions and Traps, 
# with species removed not at subsetted sites:
# https://github.com/dougwyu/HJA_analyses_Kelpie/tree/master/Kelpie_maps/outputs_minimap2_20200929_F2308_q48_kelpie20200927_BF3BR2_vsearch97

# Files beginning otuenv_M1S2, etc.. 

## Data from github..

# gitFolder <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps/outputs_minimap2_20200929_F2308_q48_kelpie20200927_BF3BR2_vsearch97/"
#   
files <- c("otuenv_M1S1_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M2S1_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M1S2_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M2S2_minimap2_20200929_kelpie20200927.csv")
# 
# dfiles <- file.path(gitFolder, files)
# 
# dataF <- lapply(dfiles, read.csv)
# 
# ## Also copied the 4 csv files to data folder.  - quicker to read in 
# mapply(function(x,y) write.csv(x, file = file.path("data",y)), dataF, files)
# rm(dfiles, gitFolder, dataF)

# train data only S1 M1
# test data S1 M2

SXY.train <- read.csv(file.path("data", files[1]))
# str(SXY.train)
cols <- colnames(SXY.train)
cols[1:120]

sppCols <- cols[grepl("__", cols)]
# which are the species columns?
range(which(cols %in% sppCols))

preds <- cols[7:101]

Y.train <- SXY.train[, sppCols]
S.train <- SXY.train[, 1:6]
X.train <- SXY.train[,preds]

head(S.train)
str(S.train)

## Where are the sites?
# Willamette National Forest 
# HJ Andrews Experimental Forest. # https://www.davidbuckleyborden.com/hja-experimental-forest
# EPSG:32610 # WGS 84 / UTM zone 10N

# library(sf)
# library(mapview)

S.sf <- sf::st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = 32610)

mv <- mapview::mapview(S.sf, 
                       map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT", "Thunderforest.Outdoors"),
                       legend = F)
mv

## Distance between sites

S.train.dist <- sf::st_distance(S.sf)
#S.train.dist
diag(S.train.dist) <- NA

# pairwise distances
min(S.train.dist, na.rm = T)
plot(sort(S.train.dist[upper.tri(S.train.dist)]))
hist(S.train.dist[upper.tri(S.train.dist)])


## Filter and explore predictors
preds
# remove scaled predictors for pairs plots, univariate analysis etc
preds[grepl("scale", preds)]

preds <- preds[!grepl("scale", preds)]

scale(X.train[,"NDVI_20180717"])[1:10]
X.train[,"nor.NDVI_20180717"][1:10] # normalised between 0,1

X.train[,"NDVI_20180717"][1:10]
X.train[,"nor.NDVI_20180717"][1:10]
X.train[,"mean.NDVI"][1:10]

X.train[,"B1_20180726"][1:10]

preds[grepl("NDVI", preds)]

# remove non normalised Veg INdices
# preds[grepl("(?<!nor)\\.NDVI|(?<!nor)\\.EVI|^NDVI|^EVI", preds, perl = T)]
preds <- preds[!grepl("^NDVI|^EVI", preds, perl = T)]
preds

# landsat bands
sort(preds[grepl("B\\d*_\\d+", preds)])

pairs(X.train[,preds[grepl("B1_", preds)]])
pairs(X.train[,preds[grepl("B3_", preds)]])
pairs(X.train[,preds[grepl("B10_", preds)]])
pairs(X.train[,preds[grepl("W_", preds)]])

pairs(X.train[, c("B1_20180717", "mean.NDVI", "mean.EVI")])


preds[grepl("B1_", preds)]

X.train[1:15, preds[grepl("B1_", preds)]]
summary(X.train[, preds[grepl("B1_", preds)]])

# values are reflectance?? DN??  corrected reflectance? Why such differnt values? But consistent anywway..

# Can bands be combined bandwise? like ndiv
X.train$B1_mean <- apply(X.train[, preds[grepl("B1_", preds)]], 1, mean)
pairs(X.train[, colnames(X.train)[grepl("B1_", colnames(X.train))]])
# doesn't matter what the value is... 

X.train$B_mean <- apply(X.train[, preds[grepl("B_", preds)]], 1, mean) ## 8???
X.train$B10_mean <- apply(X.train[, preds[grepl("B10_", preds)]], 1, mean)
X.train$B11_mean <- apply(X.train[, preds[grepl("B11_", preds)]], 1, mean)
X.train$B2_mean <- apply(X.train[, preds[grepl("B2_", preds)]], 1, mean)
X.train$B3_mean <- apply(X.train[, preds[grepl("B3_", preds)]], 1, mean)
X.train$B4_mean <- apply(X.train[, preds[grepl("B4_", preds)]], 1, mean)
X.train$B5_mean <- apply(X.train[, preds[grepl("B5_", preds)]], 1, mean)
X.train$B6_mean <- apply(X.train[, preds[grepl("B6_", preds)]], 1, mean)
X.train$B7_mean <- apply(X.train[, preds[grepl("B7_", preds)]], 1, mean)


## point way out by itself... 
pairs(X.train[, preds[grepl("NDVI", preds)]])

summary(X.train[,"mean.NDVI"])

hist(X.train[,"mean.NDVI"])

which.min(X.train[,"mean.NDVI"])
plot(sort(X.train[,"mean.NDVI"]))

S.train$ndvi_low <- X.train[,"mean.NDVI"] < 0.8

table(S.train$ndvi_low)

library(sf)
head(S.train)
sites <- st_as_sf(S.train, coords = c("UTM_E", "UTM_N"), crs = 32610)


mapview::mapview(sites, zcol = "ndvi_low", map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT", "Thunderforest.Outdoors"))


plot(X.train[,"mean.NDVI"],  X.train[,"mean.EVI"])
cor(X.train[,"mean.NDVI"],  X.train[,"mean.EVI"])


preds <- colnames(X.train)

# remvoe scaled
preds <- preds[!grepl("scale", preds)]

# remove individual sat bands
preds[grepl("^.\\d{0,2}_\\d+$", preds)]
preds <- preds[!grepl("^.\\d{0,2}_\\d+$", preds)]

# remove individual veg indices
preds[grepl("(NDVI|EVI)_\\d*", preds)]
preds <- preds[!grepl("(NDVI|EVI)_\\d*", preds)]

preds

## Look at covers
pairs(X.train[, preds[grepl("^l_", preds)]])
# all returns variables are all correlated with first returns - remove

preds <- preds[!grepl("_all", preds)]

pairs(X.train[, preds[grepl("^l_", preds)]])
cor(X.train[, preds[grepl("^l_", preds)]]) > 0.7

# remove l_95 # correlated to l_25 amnd rumple
preds <- preds[!preds %in% c("l_p95")]

preds
summary(X.train[, preds])

table(X.train$clearcut)
table(X.train$insideHJA)

getwd()
source("../code_local/pairs_helpers.r")
## all landsat data
pairs(X.train[, preds[grepl("mean", preds)]], diag.panel = panel.hist, lower.panel = panel.cor)

# out B, B2, B3, B5, B7, B10, mean.bright, 
preds <- preds[!preds %in% c("B_mean", "B2_mean", "B3_mean", 
                             "B5_mean","B6_mean", "B7_mean", "B10_mean", "B11_mean", "mean.bright")]

pairs(X.train[-45, preds[grepl("mean", preds)]], diag.panel = panel.hist, lower.panel = panel.cor)

# these extreme ndvi values (bare ground... maybe) are increasing the correlations... 
pairs(X.train[X.train[,"mean.NDVI"] > 0.8, 
              preds[grepl("mean", preds)]], diag.panel = panel.hist, lower.panel = panel.cor)



preds
pairs(X.train[, preds[3:11]], diag.panel = panel.hist, lower.panel = panel.cor)

# remove max, min temperatures and leave elevation

# log distances stream, road, years since distr

summary(log10(X.train$distToRoad_m))
summary(log(X.train$YrsSinceDist))

hist(log10(X.train$YrsSinceDist))
hist(sqrt(X.train$YrsSinceDist))

summary(X.train$YrsSinceDist)

X.train$lg_DistRoad <- log10(X.train$distToRoad_m)
X.train$lg_YrsDisturb <- log10(X.train$YrsSinceDist)

## LOG COVERS
summary(log10(X.train$l_Cover_2m_max))
summary(log10(X.train$l_Cover_2m_4m))

hist(log10(X.train$l_Cover_2m_max))
hist(log10(X.train$l_Cover_2m_4m))

summary(log10(X.train$l_Cover_4m_16m))
summary(X.train$l_Cover_4m_16m)
summary(log10(X.train$l_Cover_4m_16m+0.01))

hist(log10(X.train$l_Cover_4m_16m))
hist(X.train$l_Cover_4m_16m)
hist(log10(X.train$l_Cover_4m_16m+0.01))


X.train$lg_cover2m_max <- log10(X.train$l_Cover_2m_max)
X.train$lg_cover2m_4m <- log10(X.train$l_Cover_2m_4m)
X.train$lg_cover4m_16m <- log10(X.train$l_Cover_4m_16m+0.01)

preds <- preds[!preds %in% c("distToRoad_m", "YrsSinceDist", "maxT_annual", "minT_annual", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m")]
preds <- c(preds, "lg_DistRoad", "lg_YrsDisturb", "lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m")
preds


pairs(X.train[X.train[,"mean.NDVI"] > 0.8, preds[-1:-2]], diag.panel = panel.hist, lower.panel = panel.cor)


# cover2m_max... some outliers, and Yrs-disturb very skwewed by the 200 years... 
## years could be categorical...  
table(cut(X.train$YrsSinceDist, breaks = 3))
hist(X.train$YrsSinceDist)
## binary, gt200
X.train$YrsDist_gt00 <- X.train$YrsSinceDist > 100
table(X.train$YrsDist_gt00)

preds <- c(preds, "YrsDist_gt00")
preds
# 
# [1] "clearcut"         "insideHJA"        "oldGrowthIndex"   "elevation_m"      "canopyHeight_m"  
# [6] "precipitation_mm" "distToStream_m"   "mean.NDVI"        "mean.EVI"         "mean.green"      
# [11] "mean.wet"         "l_p25"            "l_rumple"         "B1_mean"          "B4_mean"         
# [16] "lg_DistRoad"      "lg_YrsDisturb"    "lg_cover2m_max"   "lg_cover2m_4m"    "lg_cover4m_16m"  
# "YrsDist_gt00"


## separate by remote vs field measured...  can we predict bioD from just remote sensed




