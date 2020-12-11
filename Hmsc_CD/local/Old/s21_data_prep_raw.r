## unscaled data prep ####


library(dplyr)
library(vegan)

## FROM https://github.com/dougwyu/HJA_analyses_Kelpie/blob/master/sjSDM/R-git/R/ada/1_dataprep.Rmd

baseURL <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"



# folder and csv file , as below:
# "outputs_minimap2_20200929_F2308_q48_kelpie20200927_BF3BR2_vsearch97"
# "otuenv_M1S1_minimap2_20200929_kelpie20200927.csv"

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"
trap <- "M1"
period <- "S1"


outputidxstatstabulatefolder <- paste0("outputs_minimap2_", minimaprundate, "_", samtoolsfilter, "_", samtoolsqual, "_", "kelpie", kelpierundate, "_", primer, "_", "vsearch97")

otuenvfilename <- paste0("otuenv_", trap,period,"_minimap2_", minimaprundate, "_kelpie", kelpierundate, ".csv")

otuenv <- read.csv(file.path(baseURL, outputidxstatstabulatefolder, otuenvfilename))

# non species names
non.otu.cols <- colnames(otuenv)[!grepl("__", colnames(otuenv))]

## already scaled.... 
head(otuenv[, non.otu.cols])

XY.csv <- subset(otuenv, select = c("UTM_E", "UTM_N"))
head(XY.csv) # no scaling now

# # scale XY data
# XY.csv <- otuenv %>% 
#   select(UTM_E, UTM_N) %>% 
#   scale() %>% 
#   as_tibble()

library(dplyr)

# otu.data
# keep OTUs with >=5 incidences
minocc <- 5 # set to high number (e.g. 20) for testing
otu.qp.csv <- otuenv %>% select(contains("__"))
otu.qp.csv <- otu.qp.csv[ , vegan::specnumber(otu.qp.csv, MARGIN = 2) >= minocc] 
# convert to 0/1 data
otu.pa.csv <- otu.qp.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
# otu.pa.csv <- vegan::decostand(otu.qp.csv, method = "pa")
# summcomparedf <- summary(comparedf(otu.pa.csv, otu.pa.csv_test))

otu.pa.csv[1:10, 1:5]
otu.qp.csv[1:10, 1:5]

# default is all environmental covariates: GIS + MS + Lidar
env.scale.csv <- otuenv %>% 
  select(!contains("__")) %>% 
  select(-SiteName, -trap, -period, -UTM_E, -UTM_N, 
         -clearcut, -oldGrowthIndex, -starts_with("nor")) %>% 
  mutate(insideHJA = ifelse(insideHJA == "yes", 1, 0))

names(env.scale.csv)

corrplot::corrplot(cor(env.scale.csv), method = "ellipse", type = "lower", tl.cex = 0.5)
# GIS + MS + LiDAR:  gismslidar
msdate <- c("20180717") # alternative "20180726" # date of MS data 

env.scale.csv <- env.scale.csv %>%
  select(insideHJA, elevation_m, canopyHeight_m, minT_annual, precipitation_mm, distToRoad_m, distToStream_m, YrsSinceDist, contains(all_of(msdate)), l_Cover_2m_max, l_Cover_2m_4m, l_Cover_4m_16m, l_p25, l_rumple)

head(env.scale.csv)
corrplot(cor(env.scale.csv), method = "ellipse", type = "lower", tl.cex = 0.5)

# 
# # set variables
# prepdate <- 20201030 # run date
# minocc <- 5 # minimum occupancy (incidence) per OTU
# envvar <- "gismslidar" # gismslidarmin, gismslidar, gis, ms, lidar, mslidar
# (datafolder <- glue("data_{prepdate}_{minocc}minocc_{envvar}"))
# dir_create(here("data", "kelpie_data", "for_adagpu", datafolder))
# write_csv(env.scale.csv, here("data", "kelpie_data", "for_adagpu", datafolder,  "scale.env.csv"))
# write_csv(XY.csv, here("data", "kelpie_data", "for_adagpu", datafolder, "XY.csv"))
# write_csv(otu.qp.csv, here("data", "kelpie_data", "for_adagpu", datafolder,  "otu.qp.csv"))
# write_csv(otu.pa.csv, here("data", "kelpie_data", "for_adagpu", datafolder, "otu.pa.csv"))
