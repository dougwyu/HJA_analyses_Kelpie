#### Read data on Ada  #####

## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()
getwd()

# wd set on pipeline on Ada

library(dplyr)
# library(here) # not on ADA
library(glue)

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- glue::glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")

datFile <- glue("sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_uncorr.csv")

otuenv <- read.csv(file.path(gitHub, outputidxstatstabulatefolder, datFile))

# M1S1
trap <- "M1"
period <- "S1"
otuenv <- otuenv %>% 
  dplyr::filter(trap == trap[[1]] & period == period[[1]]) 

# bring in DEM stats
load("data/demStats.rdata") # temporary location for moment... 

# keep OTUs with >=5 incidences
# original read number abundance
minocc <- 5 # set to high number (e.g. 20) for testing
otu.ab.csv <- otuenv %>% dplyr::select(contains("__"))
otu.ab.csv <- otu.ab.csv[ , colSums(otu.ab.csv > 0) >= minocc]

# log(FSL) correction and scale to quasiprobability
otu.qp.csv <- otu.ab.csv %>% 
  mutate(across(contains("__"), 
                ~ .x /(otuenv$COISpike_sum*otuenv$lysis_ratio))) %>% 
  mutate(across(contains("__"), ~ log(.x + 0.001))) %>% 
  mutate(across(contains("__"), ~ scales::rescale(.x))) # {scales}
max(otu.qp.csv) == 1 # should be TRUE

# otu.qp.csv[1:10, 1:5]

# convert to presence/absence data
otu.pa.csv <- otu.ab.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
min(colSums(otu.pa.csv)) == minocc # should be TRUE

Y.train.pa <- otu.pa.csv
Y.train.qp <- otu.qp.csv


# env covariates
# otuenv %>% 
#   dplyr::select(!contains("__"), -UTM_E, -UTM_N, -starts_with("nor")) %>% 
#   names(.)
#  [1] "SiteName"           "trap"               "period"            
#  [4] "lysis_ratio"        "COISpike_sum"       "clearcut"          
#  [7] "insideHJA"          "oldGrowthIndex"     "elevation_m"       
# [10] "canopyHeight_m"     "minT_annual"        "maxT_annual"       
# [13] "precipitation_mm"   "distToRoad_m"       "distToStream_m"    
# [16] "YrsSinceDist"       "B1_20180717"        "B2_20180717"       
# [19] "B3_20180717"        "B4_20180717"        "B5_20180717"       
# [22] "B6_20180717"        "B7_20180717"        "B10_20180717"      
# [25] "B11_20180717"       "NDVI_20180717"      "EVI_20180717"      
# [28] "B_20180717"         "G_20180717"         "W_20180717"        
# [31] "B1_20180726"        "B2_20180726"        "B3_20180726"       
# [34] "B4_20180726"        "B5_20180726"        "B6_20180726"       
# [37] "B7_20180726"        "B10_20180726"       "B11_20180726"      
# [40] "NDVI_20180726"      "EVI_20180726"       "B_20180726"        
# [43] "G_20180726"         "W_20180726"         "B1_20180802"       
# [46] "B2_20180802"        "B3_20180802"        "B4_20180802"       
# [49] "B5_20180802"        "B6_20180802"        "B7_20180802"       
# [52] "B10_20180802"       "B11_20180802"       "NDVI_20180802"     
# [55] "EVI_20180802"       "B_20180802"         "G_20180802"        
# [58] "W_20180802"         "B1_20180818"        "B2_20180818"       
# [61] "B3_20180818"        "B4_20180818"        "B5_20180818"       
# [64] "B6_20180818"        "B7_20180818"        "B10_20180818"      
# [67] "B11_20180818"       "NDVI_20180818"      "EVI_20180818"      
# [70] "B_20180818"         "G_20180818"         "W_20180818"        
# [73] "mean.NDVI"          "mean.EVI"           "mean.bright"       
# [76] "mean.green"         "mean.wet"           "mean.NDVI.scale"   
# [79] "mean.EVI.scale"     "mean.green.scale"   "mean.bright.scale" 
# [82] "mean.wet.scale"     "l_Cover_2m_max"     "l_Cover_2m_max_all"
# [85] "l_Cover_2m_4m"      "l_Cover_2m_4m_all"  "l_Cover_4m_16m"    
# [88] "l_p25"              "l_p25_all"          "l_p95"             
# [91] "l_p95_all"          "l_rumple"          

# how many otu in total?
# sum(grepl("__", colnames(otuenv)))

# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates
env.vars <- otuenv %>% 
  dplyr::select(!contains("__"), -UTM_E, -UTM_N, -starts_with("nor")) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_"),
         elevation_m = elevation_m * 0.3048, ## convert to metres???
         B1_median = apply(across(starts_with("B1_")), 1, median),
         B2_median = apply(across(starts_with("B2_")), 1, median),
         B3_median = apply(across(starts_with("B3_")), 1, median),
         B4_median = apply(across(starts_with("B4_")), 1, median),
         B5_median = apply(across(starts_with("B5_")), 1, median),
         B6_median = apply(across(starts_with("B6_")), 1, median),
         B7_median = apply(across(starts_with("B7_")), 1, median),
         B10_median = apply(across(starts_with("B10_")), 1, median),
         B11_median = apply(across(starts_with("B11_")), 1, median),
         lg_DistStream = log(distToStream_m + 0.001),
         lg_DistRoad = log(distToRoad_m + 0.001),
         lg_YrsDisturb = log(YrsSinceDist + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001)) %>%
  dplyr::select(uniqueID, clearcut,insideHJA,oldGrowthIndex, elevation_m, canopyHeight_m, precipitation_mm, minT_annual, maxT_annual, mean.NDVI, mean.EVI, mean.green, mean.wet, mean.bright, l_p25, l_rumple, B1_median, B2_median,B3_median,B4_median,B5_median,B6_median,B7_median,B10_median,B11_median,lg_DistStream, lg_DistRoad, lg_YrsDisturb, lg_cover2m_max, lg_cover2m_4m, lg_cover4m_16m) %>% 
  dplyr::left_join(y = dem_stats[,c("uniqueID", "dem500", "tri.pt")], by = "uniqueID") %>%
  mutate(across(where(is.numeric), scale), # scale here
         clearcut = factor(clearcut),
         insideHJA = factor(insideHJA)) %>%  
  dplyr::select(-uniqueID)

# str(env.vars)

# reduce variables with Variance INflation Factor, 
# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/Eco/viffer.r")
#  
# vif <- env.vars %>%
#    dplyr::select(where(is.numeric), -elevation_m) %>%
#    viffer(z = 2)
#  
# vif

# GVIF
# canopyHeight_m   1.796164
# precipitation_mm 1.887097
# minT_annual      1.892609
# mean.green       1.659580
# mean.wet         1.331805
# l_rumple         1.647866
# lg_DistStream    1.138292
# lg_DistRoad      1.256308
# lg_YrsDisturb    1.641602
# lg_cover2m_max   1.625411
# lg_cover2m_4m    1.590426
# lg_cover4m_16m   1.307569
# tri.pt           1.236871

# X.train <- env.vars %>%
#   dplyr::select(all_of(rownames(vif)),clearcut, insideHJA)

# vif.vars <- c("canopyHeight_m", "precipitation_mm","minT_annual","mean.green","mean.wet","l_rumple","lg_DistStream",  "lg_DistRoad","lg_YrsDisturb","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m", "tri.pt")


# variable selection made a S2_define_models
X.train <- env.vars

## Study design data
S.train <- otuenv %>% 
  dplyr::select(SiteName,trap,period, UTM_E, UTM_N) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_"))

head(S.train)

spp <- data.frame(species = colnames(Y.train.pa)) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                               "genus", "epithet", "BOLD", "BOLDID",
                               "size"),
                     remove = FALSE, sep = "_") %>%
  select(-empty)

head(spp)

# convert to NAs
for(c in seq_along(spp)[-1]) spp[,c] <- sub("NA", NA, spp[,c])

# Add dummy family and genus
spp$family[is.na(spp$family)] <- sprintf("fam%03d", 1:sum((is.na(spp$family))))
spp$genus[is.na(spp$genus)] <- sprintf("gen%03d", 1:sum((is.na(spp$genus))))

head(spp)

# convert to factors for ape
spp <- spp[order(spp$class, spp$order, spp$family, spp$genus),]
tax.cols <- c("class", "order", "family", "genus", "epithet", "species")
for(i in tax.cols) spp[,i] <- factor(spp[,i])

head(spp)

P <- ape::as.phylo(~class/order/family/genus/species, data = spp, collapse = F)

P$edge.length = rep(1, length(P$edge)) # make all lengths eqaul between tax levels
ape::is.rooted(P)


all(P$tip.label %in% colnames(Y.train.pa))
all(P$tip.label %in% colnames(Y.train.qp))

# save(Y.train.pa, Y.train.qp, X.train, S.train, P, file = "data/allData_vif.rdata")

rm(c, datFile, gitHub, i, kelpierundate, minimaprundate, minocc, outputidxstatstabulatefolder, period, primer, samtoolsfilter, samtoolsqual, tax.cols, trap)

rm(dem_stats, spp)
rm(otu.ab.csv, otuenv)


