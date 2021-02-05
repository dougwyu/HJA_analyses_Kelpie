#### Read data LOCAL  #####

# wd set beforehand - this is sourced only

library(dplyr)

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

# what file am i using?
basename(fn)

# when was it modified?
file.mtime(fn)

# read complete data set
otuenv <- read.csv(fn, stringsAsFactors = FALSE, na.strings = "NA")

# Filter M1S1
trap <- "M1"
period <- "S1"
otuenv <- otuenv %>% 
  dplyr::filter(trap == trap[[1]] & period == period[[1]]) 

# clean up
rm(datFile, gitHub, kelpierundate, minimaprundate, outputidxstatstabulatefolder, period, primer, samtoolsfilter, samtoolsqual, trap, fn)


# keep OTUs with >=5 incidences
minocc <- 5 # set to high number (e.g. 20) for testing
otu.qp.csv <- otuenv %>% dplyr::select(contains("__")) ## file above is already qp
otu.qp.csv <- otu.qp.csv[ , colSums(otu.qp.csv > 0) >= minocc]

# convert to presence/absence data
otu.pa.csv <- otu.qp.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
min(colSums(otu.pa.csv)) == minocc # should be TRUE

rm(minocc)

# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates
env.vars <- otuenv %>% 
  dplyr::select(!contains("__"), -UTM_E, -UTM_N, -starts_with("nor")) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_"),
         elevation_m = elevation_f * 0.3048, ## convert to metres
         canopyHeight_m = canopyHeight_f * 0.3048,
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
  dplyr::select(uniqueID, clearcut,insideHJA,oldGrowthIndex, elevation_m, canopyHeight_m, precipitation_mm, minT_annual,
                maxT_annual, mean.NDVI, mean.EVI, mean.green, mean.wet, mean.bright, l_p25, l_p95, l_rumple, B1_median,
                B2_median,B3_median,B4_median,B5_median,B6_median,B7_median,B10_median,B11_median,lg_DistStream,
                lg_DistRoad, lg_YrsDisturb, lg_cover2m_max, lg_cover2m_4m, lg_cover4m_16m, l_Cover_2m_4m,l_Cover_4m_16m,
                be10, tri, slope, twi, Nss, Ess, ht, ht.r250, ht.r500, ht.r1k, cov2_4, cov2_4.r250, cov2_4.r500, cov2_4.r1k,
                cov4_16, cov4_16.r250, cov4_16.r500, cov4_16.r1k, be500, mTopo, cut.r1k.pt,B1_20180717, B2_20180717,
                B3_20180717, B4_20180717, B5_20180717, B6_20180717, B7_20180717, B10_20180717, B11_20180717, NDVI_20180717,
                EVI_20180717, B_20180717, G_20180717, W_20180717) %>%
  mutate(
    #across(where(is.numeric), scale), # scale here # scale when defining models etc.
    clearcut = factor(clearcut),
    insideHJA = factor(insideHJA))
#dplyr::select(-uniqueID)

# str(env.vars)
# head(env.vars)
summary(env.vars)

# reduce variables with Variance INflation Factor, 
# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/Eco/viffer.r")
#  
# vif <- env.vars %>%
#    dplyr::select(where(is.numeric), -elevation_m) %>%
#    viffer(z = 2)
#  
# vif

# X.train <- env.vars %>%
#   dplyr::select(all_of(rownames(vif)),clearcut, insideHJA)

# vif.vars <- c("canopyHeight_m", "precipitation_mm","minT_annual","mean.green","mean.wet","l_rumple","lg_DistStream",  "lg_DistRoad","lg_YrsDisturb","lg_cover2m_max","lg_cover2m_4m", "lg_cover4m_16m", "tri.pt")


## Study design data
Sp.data <- otuenv %>% 
  dplyr::select(SiteName,trap,period, UTM_E, UTM_N) %>%
  mutate(uniqueID = paste(SiteName, trap, period, sep = "_"))

head(Sp.data)

spp <- data.frame(species = colnames(otu.pa.csv)) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family","genus", 
                                          "epithet", "BOLD", "BOLDID","size"),
                  remove = FALSE, sep = "_") %>%
  dplyr::select(-empty)

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


all(P$tip.label %in% colnames(otu.pa.csv))
all(P$tip.label %in% colnames(otu.qp.csv))

# save(Y.train.pa, Y.train.qp, X.train, S.train, P, file = "data/allData_vif.rdata")

rm(c, i, tax.cols, spp)
rm(otuenv)

