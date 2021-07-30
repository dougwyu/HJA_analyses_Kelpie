
library(dplyr)
library(ggplot2)

setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

## check scale of background to range of variables at sample sites

load("J:/UEA/Oregon/gis/processed_gis_data/r_oversize/newData_unscaled.rdata") # allVars, newData, indNA
head(allVars)
colnames(newData)

# remove these - just those in model
varsOut <- c("SiteName", "trap", "period", "geometry", "cut_r", "cut_msk", "cut_40msk", "insideHJA", 
              "aspect30","maxT_annual","meanT_annual","minT_annual","precipitation_mm")

vars11 <- c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
            "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
            "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
            "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
            "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA

tmp1 <- allVars %>%
  filter(period == "S1") %>%
  dplyr::select(all_of(vars11)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
  mutate(set = "samples")%>%
  arrange(predictors)
  
head(tmp1)

compVars <- newData %>%
  dplyr::select(all_of(vars11)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "predictors", values_to = "value")%>%
  mutate(set = "background")%>%
  arrange(predictors)%>%
  rbind(tmp1)

head(compVars)
save(compVars, file = "code_GIS/compVars.rdata")
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
         colH = if_else(min_samples <= valL_background, "green", "red"))%>%  # pcL_background
  select(c(1,3,5,6,8,10,12,14,15))

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
