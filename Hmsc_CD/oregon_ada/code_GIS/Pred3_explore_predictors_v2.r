
## Local
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")
setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie")
# dir()
getwd()

library(raster)
library(sf)

library(dplyr)


# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"

gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"
gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"

## Load new data for prediction and new scaled data
# load("data/newData_scaled.rdata")
# newData.sc, xy.sites.sc, allVars.sc

## too big for github - data is on ADA and locally:
load("D:/CD/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_scaled.rdata")

head(allVars.sc)
head(newData.sc)

all.vars <- c("ht30", "gt4_r30", "gt4_250", "gt4_500", "cut_r1k", "cut_r500", "cut_r250", "cut40_r1k", "cut40_r500", 
              "cut40_r250", "be30", "tri30","slope30", "Nss30", "Ess30", "twi30", "tpi250", "tpi500", "tpi1k", "l_p25",
              "l_p95", "l_rumple", "ndmi_stdDev_r100","ndmi_stdDev_r250", "ndmi_stdDev_r500", "nbr_stdDev_r100",
              "nbr_stdDev_r250", "nbr_stdDev_r500", "ndvi_p5_r100", "ndvi_p5_r250","ndvi_p5_r500", "ndvi_p50_r100",
              "ndvi_p50_r250", "ndvi_p50_r500", "ndvi_p95_r100", "ndvi_p95_r250", "ndvi_p95_r500", "ndmi_p5_r100",
              "ndmi_p5_r250", "ndmi_p5_r500", "ndmi_p50_r100", "ndmi_p50_r250", "ndmi_p50_r500", "ndmi_p95_r100",
              "ndmi_p95_r250", "ndmi_p95_r500","LC08_045029_20180726_B1", "LC08_045029_20180726_B3", 
              "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7",
              "LC08_045029_20180726_B10", "lg_DistStream", "lg_DistRoad", "lg_cover2m_max", 
              "lg_cover2m_4m", "lg_cover4m_16m")

sum(all.vars %in% colnames(newData.sc))

## filter new data for prediction to same vars and take sample (for plotting)
set.seed(99)
rndVars <- newData.sc %>%
  slice_sample(n = 500) %>%
  dplyr::select(all_of(all.vars))
dim(rndVars)

library(corrplot)

cols <- colorRampPalette(c("red", "grey70", "blue"))

corrplot(cor(rndVars[,1:10]), method = "ellipse", type= "lower", is.corr = TRUE, diag = FALSE, col = cols(5))


## plot by groups
nameList <- list(
  be_ht = c("ht30", "gt4_r30","gt4_250","gt4_500"),
  cut = c("cut_r1k","cut_r500","cut_r250","cut40_r1k","cut40_r500","cut40_r250"),
  terr = c("be30","tri30","slope30","Nss30","Ess30","twi30","tpi250", "tpi500","tpi1k"),
  lidar = c("l_p25","l_p95","l_rumple","lg_cover2m_4m","lg_cover2m_max","lg_cover4m_16m"),
  adm = c("lg_DistStream","lg_DistRoad"),# ,"insideHJA"),
  stdDevLS = c("ndmi_stdDev_r100","ndmi_stdDev_r250", "ndmi_stdDev_r500", "nbr_stdDev_r100",
               "nbr_stdDev_r250", "nbr_stdDev_r500"),
  p5_50_95_LS = c("ndvi_p5_r100", "ndvi_p5_r250","ndvi_p5_r500", "ndvi_p50_r100",
               "ndvi_p50_r250", "ndvi_p50_r500", "ndvi_p95_r100", "ndvi_p95_r250", "ndvi_p95_r500", "ndmi_p5_r100",
               "ndmi_p5_r250", "ndmi_p5_r500", "ndmi_p50_r100", "ndmi_p50_r250", "ndmi_p50_r500", "ndmi_p95_r100",
               "ndmi_p95_r250", "ndmi_p95_r500"),
  rawLS = c("LC08_045029_20180726_B1", "LC08_045029_20180726_B3", 
               "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7",
               "LC08_045029_20180726_B10")
  )

# check names
all.vars %in% unlist(nameList)

getwd()
pdf("Hmsc_CD/local/plots/var10_corplots.pdf")

op <- par(xpd = NA)

# do corrpolot by predictor groups
lapply(nameList, function(x){
  
  corrplot(cor(rndVars[,x]), method = "ellipse", type= "lower", is.corr = TRUE, 
           diag = FALSE, col = cols(5), mar = c(1,3,6,2))
  
})

dev.off()

shell.exec(file.path(getwd(), "Hmsc_CD/local/plots/var10_corplots.pdf"))

## KEEP these based on correlations

# # "ht30","gt4_r30","gt4_250"
# 
# #  "cut_r250", "cut_r1k", "cut40_r1k","cut40_r250"
# 
# #  "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k" 
# 
# # "lg_cover2m_max", "lg_cover4m_16m", "l_rumple"
# 
# #  "lg_DistStream","lg_DistRoad"
# 
# # "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500"
# 
# #  "ndvi_p5_r100","ndvi_p5_r500"         
# # "ndvi_p50_r100"   "ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500", "ndmi_p95_r100", "ndmi_p95_r500"
# 
# # "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10"



vars <- c(
  "ht30","gt4_r30","gt4_250",
  "cut_r250", "cut_r1k", "cut40_r1k","cut40_r250",
  "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k",
  "lg_cover2m_max", "lg_cover4m_16m", "l_rumple",
  "lg_DistStream","lg_DistRoad",
  "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500",
  "ndvi_p5_r100","ndvi_p5_r500",
  "ndvi_p50_r100","ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500", "ndmi_p95_r100", "ndmi_p95_r500",
  "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10")

seq(-1,1,0.1)
length(seq(-1,1,0.1))
cols.in <- c(rep("red",3), rep("grey70", 14), rep("blue",3))

seq(-1,1,length.out = 9)
cols.in <- c(rep("red",1), rep("grey70", 6), rep("blue",1))

cols.in <- cols(5)


pdf("Hmsc_CD/local/plots/var10_corplot_selection.pdf")

# corrplot(cor(rndVars[,vars]), method = "number", type= "lower", is.corr = TRUE, number.cex = 0.3, addCoefasPercent = T, cl.cex = 0.3,
#          diag = FALSE, col = cols, mar = c(1,3,6,2))
corrplot(cor(rndVars[,vars]), method = "ellipse", type= "lower", is.corr = TRUE, 
         addCoef.col = "black", number.cex = 0.3, addCoefasPercent = T,cl.cex = 0.3, 
         diag = FALSE, col = cols.in, mar = c(1,3,6,2))

dev.off()

shell.exec(file.path(getwd(), "Hmsc_CD/local/plots/var10_corplot_selection.pdf"))


