
library(dplyr)
library(ggplot2)
library(raster)

setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"

## create a reduced prediction area - convex hull around (all points + HJA) + buffer
## Load sample site points
load(file.path(gis_out, "sample_sites.rdata"))
xy.utm; xy.all.utm

## bring in raster template to check new data
load("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada/data/gis/templateRaster.rdata") # r.msk, indNA, 
rm(r.aoi.pred, aoi.pred.sf)


## check scale of background to range of variables at sample sites
load("J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_unscaled.rdata") # allVars, newData, indNA
head(allVars) # extracted predictor values at all sample point sites/
colnames(newData) # all data (without NAs) from all predictor layers across study area

# remove these 
# varsOut <- c("SiteName", "trap", "period", "geometry", "cut_r", "cut_msk", "cut_40msk", "insideHJA", 
#               "aspect30","maxT_annual","meanT_annual","minT_annual","precipitation_mm")


# Just use these -- - just those in model
## Final set of VIF chosen predictors
vars11 <- c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
            "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
            "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
            "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
            "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA
# 
# tmp1 <- allVars %>%
#   filter(period == "S1") %>%
#   dplyr::select(all_of(vars11)) %>%
#   tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
#   mutate(set = "samples")%>%
#   arrange(predictors)
#   
# head(tmp1)
# 
# compVars <- newData %>%
#   dplyr::select(all_of(vars11)) %>%
#   tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
#   mutate(set = "background")%>%
#   arrange(predictors)%>%
#   rbind(tmp1)
# 
# rm(tmp1)

# head(compVars)
# save(compVars, file = "code_GIS/compVars.rdata")
load("code_GIS/compVars.rdata")

sumVars <- compVars %>%
  group_by(predictors, set)%>%
  summarise(min = min(value),
            max = max(value),
            pcL = quantile(value, prob = 0.05),
            pcH = quantile(value, prob = 0.95),
            valL = seq(min,max, length.out = 100)[5],
            valH = seq(min,max, length.out = 100)[95])%>%
  #filter(set == "samples") %>%
  tidyr::pivot_wider(names_from = set, values_from = c(min, max, pcL, pcH, valL, valH))%>%
  mutate(colL = if_else(max_samples >= valH_background, "green", "red"),   # pcH_background
         colH = if_else(min_samples <= valL_background, "green", "red"))# pcL_background
  # select(c(1,3,5,6,8,10,12,14,15))

sumVars
#sumVars[,c(1,3,5,6,8,10,12,14,15)]


ggplot(compVars, aes(x = value))+
  geom_histogram()+
  geom_vline(data = sumVars, aes(xintercept = max_samples, colour = colH))+
  geom_vline(data = sumVars, aes(xintercept = min_samples, col = colL))+
  scale_colour_identity(guide = "legend", labels = c("outside", "within"), 
                        name = "5/95 % data value")+ # name = "5/95 % data value" name ="5/95 percentile"
  #geom_vline(data = sumVars, aes(xintercept = pcL_background), col = "blue")+
  #geom_vline(data = sumVars, aes(xintercept = pcH_background), col = "blue")+
facet_wrap(~predictors, scales = "free", ncol = 5, nrow = 6) # 

0.025* nrow(newData)

getwd()
ggsave("../local/plots/sample_background_quantile.png")
ggsave("../local/plots/sample_background_data_value.png")


### Make new data masking values out of range of sample sites.
## IS specific to variables used in model. So needs to be made for each new set of vars...

sumVars
# probably faster to copy max/min values as extra column, but used more memory in RAM>

## uses sumVars - specific to vars and can change cutoff value there.
newData_ir <- newData %>%
  dplyr::select(all_of(vars11)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value") %>%
  rowwise() %>%
  mutate(value2 = 
           case_when(
             value < sumVars$min_samples[sumVars$predictors == predictors] |
               value > sumVars$max_samples[sumVars$predictors == predictors] ~ NA_real_, # must be same type... 
             TRUE ~ value
           ))%>%
  arrange(predictors)
  

newData_clamp <- newData %>%
  dplyr::select(all_of(vars11)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value") %>%
  rowwise() %>%
  mutate(value2 = 
           case_when(
             value < sumVars$min_samples[sumVars$predictors == predictors] 
             ~ sumVars$min_samples[sumVars$predictors == predictors],
             value > sumVars$max_samples[sumVars$predictors == predictors] 
             ~ sumVars$max_samples[sumVars$predictors == predictors],
             TRUE ~ value
           ))%>%
  arrange(predictors)%>%
  as.data.frame()

newData_clamp

save(newData_ir,newData_clamp, file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_ir.rdata")
load("J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_ir.rdata")

newData[1:10, vars11]

summary(newData_ir[newData_ir$predictors == "be30",])

## Change to wide format
newData_cl_w <- do.call(data.frame, 
                        lapply(vars11, function(x) newData_clamp[newData_clamp$predictors == x, "value2"]))
colnames(newData_cl_w) <- vars11
head(newData_cl_w)

newData_ir_w <- do.call(data.frame, 
                        lapply(vars11, function(x) newData_ir[newData_ir$predictors == x, "value2"]))
colnames(newData_ir_w) <- vars11
head(newData_ir_w)

save(newData_cl_w, newData_ir_w, 
     file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_clamp_w.rdata")


## Make raster stack to check area of NAs

x <- vars11[1]

## make rasters
rList <- lapply(vars11, function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- newData_ir$value2[newData_ir$predictors == x]
  tmp
  
})

# plot(tmp)

rStack_ir <- stack(rList)
names(rStack_ir) <- vars11
rStack_ir

plot(rStack_ir[[1:10]], colNA = "black")

## sum to get all NAs
smStack_ir <- is.na(sum(rStack_ir))
plot(smStack_ir, colNA  = "black")

# how many extra NAs?
sum(indNA) - sum(is.na(values(smStack_ir)))

## Clamp stack
rList <- lapply(vars11, function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- newData_clamp$value2[newData_clamp$predictors == x]
  tmp
  
})

# plot(tmp)

rStack_clamp <- stack(rList)
names(rStack_clamp) <- vars11
rStack_clamp

plot(rStack_clamp[[1:10]], colNA = "black")

rStack_ir[[1]]
rStack_clamp[[1]]

save(rStack_ir, rStack_clamp, file = "J:/UEA/Oregon/gis/processed_gis_data/r_oversize/rStack_ir_clamp.rdata")

## Make newData (in wide format)
# load("J:/UEA/Oregon/gis/processed_gis_data/r_oversize/rStack_ir_clamp.rdata") # rStack_ir, rStack_clamp

# newData_cl_w <- raster::values(rStack_clamp)
# head(newData_cl_w)
# newData_cl_w <- data.frame(newData_cl_w[indNA, ])
# 

## Scale whole data set - apart from categorical predictors
allBrck.sc <- scale(dropLayer(allBrck, c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")))
# stores scale parameters in the @data slot
allBrck.sc # in memory
# str(allBrck.sc)
inMemory(allBrck.sc[[1]])

sapply(1:nlayers(allBrck.sc), function(x) inMemory(allBrck.sc[[x]]))

## add back categorical - but bring into memory first
catRasters <- readAll(allBrck[[c("insideHJA", "cut_r" , "cut_msk", "cut_40msk")]])
catRasters[[1]]

allBrck.sc <- addLayer(allBrck.sc, catRasters)
names(allBrck.sc)
