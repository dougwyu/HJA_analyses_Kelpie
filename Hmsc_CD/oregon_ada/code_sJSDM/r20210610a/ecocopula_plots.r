#### plots

library(dplyr)
library(raster)
library(sf)
library(ggplot2)

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

gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"

baseFolder <- "code_sjSDM/r20210610a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
abund <- "pa"


# load raster ecocopula 
load(file.path(resFolder, "ecocop_res_f5.rdata")) # mod.pa, mod.pa.ord, rSites.pa
load(file.path(resFolder, "rast_template_data.rdata"))

# ## make site scores into raster
rSites.pa_fac1 <- rSites.pa
rSites.pa_fac2 <- r


rSites.pa_fac2[] <- NA
rSites.pa_fac2[indNA2] <- mod.pa.ord$scores[,"Factor2"]
rSites.pa

cols <- colorRampPalette(c("darkblue", "white", "darkblue"))


f12<- stack(rSites.pa_fac1, rSites.pa_fac2)
names(f12) <- c("factor 1", "factor 2")

cellStats(f12, "range")

pdf(file.path(plotsFolder, "site_scores_vars8.pdf"), width = 8, height = 5)
plot(f12, col = cols(255), breaks = seq(-3, 3, length.out = 256), legend = FALSE, addfun = addAll)
dev.off()

pdf(file.path(plotsFolder, "site_scores_vars8_terrCol.pdf"), width = 8, height = 5)
plot(f12, addfun = addAll)
dev.off()


# plot(rSites.pa)

# # make into df for ggplot
# coords <- xyFromCell(rSites, seq_len(ncell(rSites)))
# df1 <- as.data.frame(values(rSites))
# df1 <- cbind(coords, df1)

# Species scores
# sp_res <- data.frame(mod.ord$loadings, species = colnames(dataN$otu.pa))

# save(mod.prob, mod.prob.ord, rSites.pa, file = file.path(resFolder, "ecocop_res_prob.rdata"))
