### Species richness plots comparing high and low prevalence species

wd <- here::here()
wd
setwd(wd)
dir()

# (follows quantiles of prevalence - as in models of 20201209 - see S2_define_models.r)

# get data
samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

gitHub <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/Kelpie_maps"

outputidxstatstabulatefolder <- paste0("outputs_minimap2_",minimaprundate,"_",samtoolsfilter,"_",samtoolsqual,"_kelpie",kelpierundate,"_",primer,"_vsearch97")

datFile <- paste0("sample_by_species_table_", samtoolsfilter, "_minimap2_", minimaprundate, "_kelpie", kelpierundate, "_uncorr.csv")

otuenv <- read.csv(file.path(gitHub, outputidxstatstabulatefolder, datFile))

# M1S1
trap <- "M1"
period <- "S1"
otuenv <- otuenv %>% 
  dplyr::filter(trap == trap[[1]] & period == period[[1]]) 

# bring in DEM stats - 
# load("Hmsc_CD/oregon_ada/data/demStats.rdata") # temporary location for moment... 
# head(dem_stats) # dem500: averaged elevation over 500 m; tri.pt: topographic roughness index at each sample location

# keep OTUs with >=5 incidences
# original read number abundance
minocc <- 0 # set to high number (e.g. 20) for testing
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

## Make two data sets, with < and >= 5 occurrences
colLT5 <- colSums(otu.pa.csv) < 5
colGT5 <- colSums(otu.pa.csv) >=5

sum(colLT5) # 942 spp have less than 5
sum(colGT5) # 268 

sum(colLT5, colGT5) == ncol(otu.pa.csv)

# make groups
otu.lt5 <- otu.pa.csv[, colLT5]
otu.gt5 <- otu.pa.csv[, colGT5]

# correlation
corT <- cor.test(rowSums(otu.lt5),rowSums(otu.gt5), method = "spearman")
str(corT)

# scatter plot of lowest prevalence species against all other quantiles
png("Hmsc_CD/local/plots/spRich_cor_minocc5.png", width = 250, height = 200, units= "mm", res = 200)
plot(rowSums(otu.lt5),rowSums(otu.gt5), 
     xlab = "sp richness (< 5 occurrences)", 
     ylab = "sp richness (>= 5 occurrences)" , pch = 16)
text(32, 70, paste("Spearman correlation\np =", round(corT$p.value, 5),
                   "rho =", round(corT$estimate,3),
                   "n =", nrow(otu.pa.csv)), adj = 0)
dev.off()

# correlation test
corT


# Spearman's rank correlation rho
# 
# data:  rowSums(otu.lt5) and rowSums(otu.gt5)
# S = 72749, p-value = 0.000584
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3594023 

