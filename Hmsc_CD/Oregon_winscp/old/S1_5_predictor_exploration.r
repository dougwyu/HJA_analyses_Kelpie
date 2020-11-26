### Exploring predictors

# setwd("J:/UEA/Oregon")

sP

## Only local: 
setwd("J:/UEA/Oregon")
dir()
setwd("~/Oregon_winscp")


## Check coordinates
load("data/allData.Rdata")

head(S)
# Oregon UTM is 10N
# Willamette National Forest 
# HJ Andrews Experimental Forest. # https://www.davidbuckleyborden.com/hja-experimental-forest


# EPSG:32610
# WGS 84 / UTM zone 10N

addmargins(table(SXY[,c("trap", "session")]))


library(sf)
library(mapview)

S.sf <- st_as_sf(S, coords = c("Route_x", "Route_y"), crs = 32610)

mv <- mapview(S.sf, zcol = c("session"), 
              col.regions = c("blue", "green"), 
              map.types = c("Esri.WorldImagery", "OpenStreetMap.HOT", "Thunderforest.Outdoors"))


getwd()
mapshot(mv, "../Oregon_notes/mapview_oregon.html")

## Distance between sites

S.train.dist <- st_distance(subset(S.sf, session == "S1" & trap == "M1"))
S.train.dist
diag(S.train.dist) <- NA

min(S.train.dist, na.rm = T)
plot(sort(S.train.dist[upper.tri(S.train.dist)]))
hist(S.train.dist[upper.tri(S.train.dist)])





## predictor choices - see exel doc for description
# cbind(X_cols_use)
# [1,] "elevation"         
# [2,] "canopy.ht"         
# [3,] "min.T"             
# [4,] "max.T"             
# [5,] "precipitation"     
# [6,] "metre.road"        -- log
# [7,] "metre.stream"      ---log
# [8,] "yrs.disturb.min"   ---log
# [9,] "hja"    # hja inside (no logging really, nearer to primary forest) or outside experimental forest (logging)
# make a domindant land  cover variable around each point. % of structure... 

# [10,] "lysis.ratio"   #  ignore    lysis buffer ratio. different samples are different sizes - control for this?
# [11,] "spike.mean"        # 
# [12,] "l_Cover_2m_4m"     
# [13,] "l_Cover_2m_4m_all" 
# [14,] "l_Cover_2m_max"    
# [15,] "l_Cover_2m_max_all"
# [16,] "l_Cover_4m_16m"    # understorey
# [17,] "l_p25"             
# [18,] "l_p25_all"         
# [19,] "l_p95"             
# [20,] "l_p95_all"         
# [21,] "l_rumple"          
# [22,] "mean.NDVI"         
# [23,] "mean.EVI"          
# [24,] "mean.bright"       
# [25,] "mean.green"        
# [26,] "mean.wet" 