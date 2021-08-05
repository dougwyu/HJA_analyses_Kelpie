## TSNE Correlation #####


library(corrplot)


# only on local
# getwd()
# wd <- here::here()
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

## on ADA
gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"

# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210716b"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)

#### Load TSNE data 
load(file.path(resFolder, "ord_tsne_res_cl_p50.rdata"))
# tsne, r, rSites1, rSites2, NAs, 

str(tsne, max.level = 1)

### Load irreplaceability data ##



## Load predictor data
## Load new data for prediction and new scaled data
load(file.path(gis_out, "r_oversize/newData_unscaled.rdata")) # allVars, newData, indNA, 

## Final set of VIF chosen predictors
vars11 <- c("gt4_500", "cut_r1k","cut_r250","cut40_r1k","cut40_r250","be30","tri30","Nss30",
            "Ess30","twi30","tpi250","tpi1k","l_rumple","nbr_stdDev_r100","ndvi_p5_r100",
            "ndvi_p5_r500","ndvi_p50_r100","ndvi_p50_r500","ndmi_p95_r100",
            "LC08_045029_20180726_B1","LC08_045029_20180726_B5","lg_DistStream",
            "lg_DistRoad","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m") # insideHJA


dim(tsne$Y)
dim(newData[,vars11])

mod.cor <- cor(tsne$Y, newData[,vars11])
mod.cor

getwd()
pdf(file.path(plotsFolder, "tsne_predictors_corrplot.pdf"), width = 8, height = 8)
corrplot::corrplot(t(mod.cor), 
           is.corr = F,
           method = "ellipse",
           cl.pos = "n",
           col.lim = c(-1,1),
           title = "tsne",
           oma = c(0,0,0,0),
           #oma = c(2,2,5,1),
           mar = c(0,0,1,0),
           addCoef.col = "black",
           addCoefasPercent = T,
           number.cex = 0.5)

dev.off()

shell.exec(file.path(getwd(), plotsFolder, "tsne_predictors_corrplot.pdf"))

## scatter plots

library(ggplot2)
library(dplyr)

set.seed(45)
ind <- sample(1:nrow(newData), size = 5000)

# p1 <- newData %>%
#   select(all_of(vars11))%>%
#   slice(ind)%>%
#   #rename_with(~gsub("(.*)", "\\1__", .x))%>%
#   mutate(tsne1 = tsne$Y[ind,1],
#          tsne2 = tsne$Y[ind,2])%>%
#   tidyr::pivot_longer(cols = !contains("tsne"), names_to = "predictor", values_to = "pred_val")%>%
#   tidyr::pivot_longer(cols = contains("tsne"), names_to = "tsne", values_to = "tsne_val")%>%
#   ggplot(aes(x = pred_val, y = tsne_val))+
#   geom_point(shape = ".")+
#   facet_grid(cols = vars(predictor), rows = vars(tsne), scales = "free", space="free")
# 
# p1
# 
# ggsave(file.path(plotsFolder, "tsne_predictors_scatterplot.pdf"), p1)
# shell.exec(file.path(getwd(), plotsFolder, "tsne_predictors_scatterplot.pdf"))

pdf(file.path(plotsFolder, "tsne_predictors_scatterplot.pdf"), width = 2.5, height = 30)
par(mfcol=c(length(vars11), ncol(tsne$Y)), mar = c(2,2,2,2))
sapply(vars11, function(x) smoothScatter(
  newData[ind,x], tsne$Y[ind,1], 
  xlab = "", ylab = "", main = x, cex.main = 0.5, axes = T, cex.axis = 0.4, cex = 0.1, pch = 16))
sapply(vars11, function(x) smoothScatter(
  newData[ind,x], tsne$Y[ind,2], 
  xlab = "", ylab = "", main = x, cex.main = 0.5, axes = T, cex.axis = 0.4, cex = 0.1, pch = 16))

dev.off()

shell.exec(file.path(getwd(), plotsFolder, "tsne_predictors_scatterplot.pdf"))
