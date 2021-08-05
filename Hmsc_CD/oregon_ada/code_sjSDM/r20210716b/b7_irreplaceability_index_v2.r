## Make raster from predictions


#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD/oregon_ada"))
# dir()

## trial run 
options(echo=TRUE) # if you want see commands in output file

library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(ggplot2)

utm10N <- 32610

## on ADA
# gis_out <- gis_in <- "data/gis"


######### Load data  ##############


## local
gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210716b"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

# load raster templates - reduced areas
load(file.path(gis_out, "templateRaster.rdata")) ## r.msk, indNA aoi.pred.sf, r.aoi.pred - reduced area for plotting
# plot(r.msk)
# plot(aoi.pred.sf)


# load model data - for species classification
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_test_results.rdata")) # eval.results, sp.res.test, sp.res.train


########### Filter species for analysis ##########


## Mean AUC per species (and other eval metrics)
str(sp.res.test, max.level = 1)
head(sp.res.test$auc)

sum(is.na(sp.res.test$auc))

## Filter species by auc
auc.filt <- 0.70
sum(sp.res.test$auc > auc.filt, na.rm = T)

# ## extract species over AUC filter
# str(pred.sp, max.level = 1)

# incidence 
incidence <- colSums(otu.pa.csv)/nrow(otu.pa.csv)

# get species names, BOLD ID etc
spp <- data.frame(species = colnames(get(paste0("otu.", abund, ".csv")))) %>%
  tidyr::separate(col = species, into = c("OTU", "empty", "class", "order", "family",
                                          "genus", "epithet", "BOLD", "BOLDID",
                                          "size"),
                  remove = FALSE, sep = "_", convert = TRUE) %>%  ## creates real NAs with convert = T
  mutate(best.name = case_when(is.na(epithet) & is.na(genus) & is.na(family) & is.na(order) ~ class,
                               is.na(epithet) & is.na(genus) & is.na(family) ~ order,
                               is.na(epithet) & is.na(genus) ~ family,
                               is.na(epithet) ~ genus,
                               TRUE ~ paste(genus, epithet, sep = "_")
  )) %>%
  dplyr::select(-empty)%>%
  mutate(auc = sp.res.test$auc,
         incidence = incidence,
         best.name = paste(best.name, BOLDID, sep = "_"))

head(spp)

sum(is.na(spp$best.name))
sum(grepl("NA_NA", spp$best.name))
head(spp, 30)

sum(is.na(spp$family))


### Load prediction results
## In ADA
# load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) # pred.mn, pred.sd, 

## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
load(file.path(gis_out, "r_oversize", paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata"))) 
paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")

dim(pred.mn)

pred.in <- pred.mn[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
dim(pred.in)

## get species names too
spp.in <- spp[sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc), ]
head(spp.in)


##### Make rasters for 'site' based analysis ('summarised by blocks') ########

## make rasters
# plot(r.aoi.pred)
# x <- data.frame(pred.in)[,1]

rList <- lapply(data.frame(pred.in), function(x) {
  
  tmp <- r.msk
  tmp[indNA] <- x
  tmp
  
})

# plot(tmp)

rStack <- stack(rList)
# names(rStack) <- spp.in$best.name
rStack
## add auc incidence names to stack
names(rStack) <- paste0(spp.in$best.name, " ", "auc=", round(spp.in$auc, 2), " ",  "prev=", round(spp.in$incidence,2))


# Use a small threshold to remove all cells with probabilities less than 0.1, but remaining cells 
## retain their raw probablity value
tr2 <- 0.1
rStack.ct <- raster::reclassify(rStack, rcl = c(0, tr2, 0))
plot(rStack.ct, 1:6)

save(rStack, rStack.ct, file = file.path(resFolder, "spStacks.rdata"))
load(file.path(resFolder, "spStacks.rdata"))

sapply(1:nlayers(rStack),function(x) raster::inMemory(rStack[[x]]))
rm(rStack); gc()

## Create block grid
plot(r.msk)

rBlock <- r.msk # from template raster

## aggregate to 120m res (4 x 30 m), 4 x 4 cells = 1 site
rBlock <- aggregate(r.msk, fact = 4)
## assign block id numbers
rBlock[] <- 1:ncell(rBlock)
names(rBlock) <- "ID"

# assign block values to original raster template
rBlckID <- disaggregate(rBlock, 4, method = "")

# plot small region to check
# xlim <- c(xmin(rBlckID), xmin(rBlckID) + res(rBlckID)[1]*20)
# ylim <- c(ymax(rBlckID) - res(rBlckID)[1]*20, ymax(rBlckID))
# plot(rBlckID, col = sample(rainbow(2550)), xlim = xlim, ylim = ylim)
# text(rBlckID)

# ## Check raster aggregate
# n <- 20
# r1 <- raster(matrix(1:n^2, ncol = n),xmn=0, xmx=n, ymn=0, ymx=n)
# r1
# 
# r2 <- aggregate(r1, fact = 5)
# r2
# r2[] <- 1:ncell(r2)
# 
# plot(r1, col = sample(rainbow(255)))
# plot(r2)
# 
# r3 <- disaggregate(r2, 5)
# r3
# plot(r3)
# plot(stack(r1, r3))

# plot(r.msk, xlim = c(555000,560000), ylim = c(4905000,4907500), colNA = "black")
# plot(rBlckID, add = T, col = sample(gray(seq(0,1, 0.00001), 0.2)))


# add to stack
rProb <- stack(rStack, rBlckID)
rCt <- stack(rStack.ct, rBlckID)

# Get matrix of cells(sites) x species + last column
rValsProb <- values(rProb)
rValsCT <- values(rCt)

save(rValsProb, rValsCT, file= file.path(resFolder, "rVals.rdata"))
load(file.path(resFolder, "rVals.rdata"))

dim(rVals)
# head(rVals)
rVals[2000:2010, 80:ncol(rVals)]
rVals[1:36, 80:ncol(rVals)]

sum(rVals[rVals[,"ID"] == 1, 1], na.rm = T)

## Summarise by block
blckValsProb <- aggregate(rValsProb[,1:(ncol(rValsProb)-1)], list(rValsProb[,"ID"]), sum, na.rm = F)
blckValsCT <- aggregate(rValsCT[,1:(ncol(rValsCT)-1)], list(rValsCT[,"ID"]), sum, na.rm = F)

save(blckValsProb, blckValsCT, file = file.path(resFolder, "blckVals.rdata"))
load(file.path(resFolder, "blckVals.rdata"))

# blckVals[2000:2010, 80:ncol(blckVals)]
# str(blckVals)
# colnames(blckVals)
# range(blckVals[,2:88], na.rm= T) #] max of 16 (from 4x4 aggregation of original raster)
# hist(unlist((blckVals[,2:88])))
# hist(unlist((blckValsProb[,2:88])))
# range(blckVals$Group.1)

# remove ID col (group.1 - first column). Each row is a cell/site - convert to matrix.
x_in_prob <- as.matrix(blckValsProb[, 2:ncol(blckVals)])
x_in_ct <- as.matrix(blckValsCT[, 2:ncol(blckValsCT)])

dim(x_in)
dim(x_in_prob)

x_in_prob[1:100, 1:10]
blckValsProb[1:100, 1:10]

# Rs distribution
hist(colSums(x_in_ct, na.rm = T))
hist(colSums(x_in_prob, na.rm = T))

## Median Rs, range
range(colSums(x_in, na.rm = T))
median(colSums(x_in, na.rm = T))

range(colSums(x_in_prob, na.rm = T))
median(colSums(x_in_prob, na.rm = T))

median(colSums(x_in_ct, na.rm = T))

# ## total available area:
sum(complete.cases(rVals)) == sum(!is.na(values(r.msk)))
# 247,743 pixels (30x30m)

# Max prob sum is ~175,000  max would be 1 in every available pixel

beta.r.prob <- irrAB(x = x_in_prob, pc = 0.90, type = "total", r = rBlock)
writeRaster(beta.r.prob, file = file.path(resFolder, "beta_r_prob.tif"))

tmp_stck <- stack(beta.r.bin, beta.r.prob)
names(tmp_stck) <- c("beta_binary_90pc", "beta_prob_90pc")
plot(tmp_stck)

beta.r.prob_ct <- irrAB(x = x_in_ct, pc = 0.90, type = "total", r = rBlock)
plot(beta.r.prob_ct)
writeRaster(beta.r.prob_ct, file = file.path(resFolder, "beta_r_prob_ct.tif"))

