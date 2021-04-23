
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


# load as brick
allBrck <- brick(file.path(gis_out, "r_utm/allStack.tif"))

# get names and name groups
load(file.path(gis_out, "allNames.rdata"))
names(allBrck) <- allNames
allBrck

rm(allNames)

## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm

## OR create random points 

# check extent of NAs across all rasters
freq(allBrck$cut_r)  # only layer with lots of NAs

# get index to remove NAs from data (without year since cut, as no year is NA) . 
## OJOO More efficient to make brick without this.
indNA <- complete.cases(values(dropLayer(allBrck, "cut_r")))
sum(indNA)
# rNA <- sum(dropLayer(allBrck, "cut_r")) > 0

# Plot to check extent of all predictors with values:
rNA <- raster(allBrck)
rNA[] <- indNA

plot(rNA)
plot(xy.utm, add = T, pch = 16, col = "black")
rm(rNA)

## make random selection of points to extract values
coords <- coordinates(allBrck)
coords <- coords[indNA,]

set.seed(1757)
coords <- coords[sample(seq_len(nrow(coords)), size = 500), ]

rndVars <- extract(allBrck, coords)
head(rndVars)

save(rndVars, file = file.path("Hmsc_CD/oregon_ada/data", "rndVars.rdata"))

# remove insideHJA, and binary, and year since distrubt
rndVars <- subset(rndVars, select = -c(insideHJA, cut_r, cut_msk, cut_40msk))

library(corrplot)

cols <- colorRampPalette(c("red", "grey90", "blue"))

corrplot(cor(rndVars[,1:10]), method = "ellipse", type= "lower", is.corr = TRUE, diag = FALSE, col = cols(5))

# remove names from nameList
x <- nameList[[2]]

nameList2 <- lapply(nameList, function(x){
    x[! x %in% c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")]
  }) 

getwd()
pdf("Hmsc_CD/local/plots/var_corplots.pdf")

op <- par(xpd = NA)

# do corrpolot by predictor groups
lapply(nameList2, function(x){
  
  corrplot(cor(rndVars[,x]), method = "ellipse", type= "lower", is.corr = TRUE, 
           diag = FALSE, col = cols(5), mar = c(1,3,6,2))
  
})

dev.off()

shell.exec(file.path(getwd(), "Hmsc_CD/local/plots/var_corplots.pdf"))

## KEEP 

# "ht30","gt4_r30","gt4_250"
#  "cut_r250","cut40_r1k","cut40_r250"
#  "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k" 
# "l_Cover_2m_max", "l_Cover_4m_16m", "l_rumple"
#  "DistStream","DistRoad"

# "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500", "ndvi_p5_r100","ndvi_p5_r500"         
# "ndvi_p50_r100"   "ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500"
# "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10"

# "minT_annual"      "precipitation_mm"


vars <- c("ht30","gt4_r30","gt4_250",
          "cut_r250","cut40_r1k","cut40_r250", 
          "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k",
          "l_Cover_2m_max", "l_Cover_4m_16m", "l_rumple",
          "DistStream","DistRoad",
          "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500", "ndvi_p5_r100","ndvi_p5_r500",
          "ndvi_p50_r100","ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500",
          "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10",
          "minT_annual","precipitation_mm")


pdf("Hmsc_CD/local/plots/var_corplot_selection.pdf")

corrplot(cor(rndVars[,vars]), method = "ellipse", type= "lower", is.corr = TRUE, 
         diag = FALSE, col = cols(5), mar = c(1,3,6,2))


dev.off()

shell.exec(file.path(getwd(), "Hmsc_CD/local/plots/var_corplot_selection.pdf"))




