### Ordination plots 


library(raster)
# library(dplyr)

## only on local
getwd()
wd <- here::here()
setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
getwd()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

## on ADA
# gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210627a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
abund <- "pa"

dir(resFolder)
dir(file.path(gis_out, "r_oversize"))

## Load ordination results
load(file.path(resFolder, "ord_tsne_res.rdata")) # 
load(file.path(resFolder, "ord_tsne_res_30.rdata")) # 
## load(file.path(resFolder, "ord_pca_res.rdata")) #  too large for github
load(file.path(gis_out, "r_oversize", "ord_pca_res.rdata")) #  too large for github

## plot ordination axes

plot(tsne$Y, pch = ".", asp = 1)
vegan:::plot.cca(pca, scaling = "sites")

## function to load site scores into raster
makeR <- function(r, siteScores, NAs) {
  
  rSites <- raster(r)
  rSites[] <- NA
  rSites[NAs] <- siteScores
  rSites
  
}

rtsne1 <- makeR(r, tsne$Y[,1], NAs)
rtsne2 <- makeR(r, tsne$Y[,2], NAs)

rtsne1_5c <- makeR(r, tsne$Y[,1], NAs)
rtsne2_5c <- makeR(r, tsne$Y[,2], NAs)

pcaR1 <- makeR(r, vegan::scores(pca, 1, "sites"), NAs)
pcaR2 <- makeR(r, vegan::scores(pca, 2, "sites"), NAs)

stck <- raster::stack(rtsne1, rtsne2, pcaR1, pcaR2)
names(stck) <- c("tsne1", "tsne2", "pca1", "pca2")

pdf(file.path(plotsFolder, "ordination_plots.pdf"), width = 8, height = 8)
plot(stck)
dev.off()


stck <- raster::stack(rtsne1, rtsne2, rtsne1_5c, rtsne2_5c)
pdf(file.path(plotsFolder, "ordination_plots_tsne.pdf"), width = 8, height = 8)
plot(stck)
dev.off()
