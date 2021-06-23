#### Read data on Ada  #####

## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()
getwd()

# wd set on pipeline on Ada

library(dplyr)

 # model settings:
abund <- "pa"
period <- "S1"
trap <- "M1"

# name <- "Dip"
name <- NULL
spChoose <- paste0(period, trap, name)

### 1. Get data from github #####

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",
                                       minimaprundate,"_",
                                       samtoolsfilter,"_", 
                                       samtoolsqual, 
                                       "_kelpie", 
                                       kelpierundate,
                                       "_",
                                       primer,
                                       "_vsearch97")

datFile <- paste0("sample_by_species_table_", 
                  samtoolsfilter, 
                  "_minimap2_",
                  minimaprundate,
                  "_kelpie",
                  kelpierundate,
                  "_FSL_qp.csv")

# file path:
fn <- file.path(gitHub, outputidxstatstabulatefolder, datFile)

# what file am i using?
basename(fn)

# when was it modified? - only if stored locally. 
file.mtime(fn)

# read complete data set
otuenv <- read.csv(fn, stringsAsFactors = FALSE, na.strings = "NA")

# keep OTUs with >= minocc incidences AND with presnece at both M1 or M2
minocc <- 6 # set to high number (e.g. 20) for testing

## get Species columns by M1 and M2, with minocc calculated per trap
## can choose below whether to include just species shared between M1 and M2
spM <- otuenv %>% 
  dplyr::filter(period == !!period) %>%
  dplyr::select(SiteName, trap, contains("__")) %>%
  tidyr::pivot_longer(cols = contains("__"), names_to = "OTU", values_drop_na = FALSE) %>%
  mutate(value = value>0) %>% # change to PA
  group_by(OTU, trap) %>%
  summarise(nSites = sum(value, na.rm = T)) %>% # Number of sites at which present
  filter(nSites >= minocc) %>% # filter by minocc
  ungroup() %>%
  tidyr::pivot_wider(names_from = trap, values_from = nSites, values_fn = function(x) sum(x)>0) %>%
  filter(M1) %>% # CHOOOSE HERE FOR SINGLE. OR SHARED TRAP SPECIES GFROUP: filter(M1 & M2)
  # tidyr::separate(col = OTU, into = c("ID", "empty", "class", "order", "family",
  #                                     "genus", "epithet", "BOLD", "BOLDID", "size"),
  #                 remove = FALSE, sep = "_", convert = TRUE) %>%
  # #dplyr::filter(order == "Diptera")%>%
  dplyr::select(OTU)

nrow(spM)


# filter species here to those in sp.M1m2$OTU - already filtered for minocc
otu.qp.csv <- otuenv %>% 
  dplyr::filter(period == !!period & trap == !!trap) %>%
  dplyr::select(spM$OTU) ## 

# convert to presence/absence data
otu.pa.csv <- otu.qp.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
min(colSums(otu.pa.csv)) >= minocc # should be TRUE

Y.train.pa <- otu.pa.csv
Y.train.qp <- otu.qp.csv

## clean up
rm(datFile, gitHub, kelpierundate, minimaprundate, outputidxstatstabulatefolder, primer, samtoolsfilter, samtoolsqual, fn, spM)

# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates

## Load new version of vars
load("data/envVars.rdata")
# load(file.path(resFolder, "data/envVars.rdata"))

head(allVars)
# colnames(allVars)

env.vars <- otuenv %>% 
  dplyr::filter(period == !!period & trap == !!trap) %>%
  dplyr::select(trap, period, UTM_E, UTM_N, SiteName) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_")) %>%
  left_join(y = allVars, by = "SiteName") %>%
  mutate(lg_DistStream = log(DistStream + 0.001),
         lg_DistRoad = log(DistRoad + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001))

## All vars
# [1] "SiteName"         "aspect30"             "cut_msk"   "cut_40msk"               

## total no of vars: 90
# "ht30"                     "gt4_r30"                  "gt4_250"                 
# [5] "gt4_500"     "cut_r"                               
# [9] "cut_r1k"                  "cut_r500"                 "cut_r250"                 "cut40_r1k"               
# [13] "cut40_r500"               "cut40_r250"               "be30"                     "tri30"                   
# [17] "slope30"                  "Nss30"                    "Ess30"                   
# [21] "twi30"                    "tpi250"                   "tpi500"                   "tpi1k"                   
# [25] "l_Cover_2m_4m"            "l_Cover_2m_max"           "l_Cover_4m_16m"           "l_p25"                   
# [29] "l_p95"                    "l_rumple"                 "DistStream"               "DistRoad"                
# [33] "insideHJA"                "ndmi_stdDev_r100"         "ndmi_stdDev_r250"         "ndmi_stdDev_r500"        
# [37] "nbr_stdDev_r100"          "nbr_stdDev_r250"          "nbr_stdDev_r500"          "ndvi_p5_r100"            
# [41] "ndvi_p5_r250"             "ndvi_p5_r500"             "ndvi_p50_r100"            "ndvi_p50_r250"           
# [45] "ndvi_p50_r500"            "ndvi_p95_r100"            "ndvi_p95_r250"            "ndvi_p95_r500"           
# [49] "ndmi_p5_r100"             "ndmi_p5_r250"             "ndmi_p5_r500"             "ndmi_p50_r100"           
# [53] "ndmi_p50_r250"            "ndmi_p50_r500"            "ndmi_p95_r100"            "ndmi_p95_r250"           
# [57] "ndmi_p95_r500"            "LC08_045029_20180726_B1"  "LC08_045029_20180726_B3"  "LC08_045029_20180726_B4" 
# [61] "LC08_045029_20180726_B5"  "LC08_045029_20180726_B7"  "LC08_045029_20180726_B10" "maxT_annual"             
# [65] "meanT_annual"             "minT_annual"              "precipitation_mm" 


# varsName <- "vars6"
# vars <- c("ht30","gt4_r30","gt4_250",
#           "cut_r250","cut40_r1k","cut40_r250", 
#           "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k",
#           "l_Cover_2m_max", "l_Cover_4m_16m", "l_rumple",
#           "DistStream","DistRoad",
#           "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500", "ndvi_p5_r100","ndvi_p5_r500",
#           "ndvi_p50_r100","ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500",
#           "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10",
#           "minT_annual","precipitation_mm")


# varsName <- "vars8" # be30 replaces minT
# vars <- c("Ess30", "DistRoad", "Nss30", "DistStream", "tri30", "LC08_045029_20180726_B5", "l_rumple", "cut_40msk", "insideHJA", "cut_msk", "cut40_r1k", "l_Cover_2m_4m", "l_Cover_4m_16m", "cut_r1k", "be30", "ndmi_p95_r100", "ndvi_p50_r500", "cut40_r250", "nbr_stdDev_r250", "precipitation_mm", "gt4_250", "tpi250", "cut_r250", "LC08_045029_20180726_B1", "ndvi_p50_r100", "gt4_r30", "twi30", "ndvi_p5_r100", "tpi1k", "ndvi_p5_r500", "l_p25", "LC08_045029_20180726_B10")

# str(env.vars)
# head(env.vars)
# cat(paste(colnames(env.vars), collapse = '", "'))

# varsName <- "vars6t10"
# vars <- c("gt4_r30",
#           "cut_r250","cut40_r1k",
#           "be30", "Nss30","tpi250",
#           "l_Cover_4m_16m", 
#           "DistStream","DistRoad",
#           "ndmi_stdDev_r100","nbr_stdDev_r100",  "ndvi_p5_r100")

varsName <- "vars8t10"
vars <- c("be30","Ess30","LC08_045029_20180726_B5", "tri30","cut_r250","cut_msk","Nss30","precipitation_mm","DistRoad","nbr_stdDev_r250")

 
# check names 
all(vars %in% colnames(env.vars))

## Save model paraameters- data saved in pipeline... 
save(otu.pa.csv, otu.qp.csv, otuenv, env.vars,spChoose,minocc, vars, varsName, abund, trap, period,
    file = file.path(resFolder, "data", paste0("modelParams_",abund,".rdata")))


##  variable selection made a S2_define_models
X.train <- env.vars[,vars]

## Study design data
S.train <- otuenv %>% 
  dplyr::select(SiteName,trap,period, UTM_E, UTM_N) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_")) %>%
  filter(period == !!period, trap == !!trap)

head(S.train)

spp <- data.frame(species = colnames(Y.train.pa)) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                               "genus", "epithet", "BOLD", "BOLDID",
                               "size"),
                     remove = FALSE, sep = "_", convert = TRUE) %>%
  dplyr::select(-empty)

head(spp)

# convert to NAs
#for(c in seq_along(spp)[-1]) spp[,c] <- sub("NA", NA, spp[,c])

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

rm(i, tax.cols, spp)
rm(otuenv)


