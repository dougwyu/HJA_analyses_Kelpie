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


# threshold?
tr <- 0.5

rStack.bin <- raster::reclassify(rStack, rcl = c(0, tr, 0, tr, 1, 1))

tr2 <- 0.2
rStack.ct <- raster::reclassify(rStack, rcl = c(0, tr2, 0))
plot(rStack.ct, 1:6)

save(rStack.bin, rStack, rStack.ct, file = file.path(resFolder, "spStacks.rdata"))
load(file.path(resFolder, "spStacks.rdata"))

sapply(1:nlayers(rStack),function(x) raster::inMemory(rStack[[x]]))
sapply(1:nlayers(rStack.bin),function(x) raster::inMemory(rStack[[x]]))
rm(rStack); gc()

## Create block grid
plot(r.msk)

rBlock <- r.msk # from template raster

## aggregate to 120m res
rBlock <- aggregate(r.msk, fact = 4)
## assign block id numbers
rBlock[] <- 1:ncell(rBlock)
names(rBlock) <- "ID"

# assign block values to original raster template
rBlckID <- disaggregate(rBlock, 4, method = "")

# plot small region to check
xlim <- c(xmin(rBlckID), xmin(rBlckID) + res(rBlckID)[1]*20)
ylim <- c(ymax(rBlckID) - res(rBlckID)[1]*20, ymax(rBlckID))

plot(rBlckID, col = sample(rainbow(2550)), xlim = xlim, ylim = ylim)
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

plot(r.msk, xlim = c(555000,560000), ylim = c(4905000,4907500), colNA = "black")
plot(rBlckID, add = T, col = sample(gray(seq(0,1, 0.00001), 0.2)))


# add to stack
rBin <- stack(rStack.bin, rBlckID)
rProb <- stack(rStack, rBlckID)
rCt <- stack(rStack.ct, rBlckID)

# Get matrix of cells(sites) x species + last column
rVals <- values(rBin)
rValsProb <- values(rProb)
rValsCT <- values(rCt)

save(rVals, rValsProb, rValsCT, file= file.path(resFolder, "rVals.rdata"))
load(file.path(resFolder, "rVals.rdata"))

dim(rVals)
# head(rVals)
rVals[2000:2010, 80:ncol(rVals)]
rVals[1:36, 80:ncol(rVals)]

sum(rVals[rVals[,"ID"] == 1, 1], na.rm = T)

## Summarise by block
blckVals <- aggregate(rVals[,1:(ncol(rVals)-1)], list(rVals[,"ID"]), sum, na.rm = F)
blckValsProb <- aggregate(rValsProb[,1:(ncol(rValsProb)-1)], list(rValsProb[,"ID"]), sum, na.rm = F)
blckValsCT <- aggregate(rValsCT[,1:(ncol(rValsCT)-1)], list(rValsCT[,"ID"]), sum, na.rm = F)

save(blckVals, blckValsProb, blckValsCT, file = file.path(resFolder, "blckVals.rdata"))

blckVals[2000:2010, 80:ncol(blckVals)]
str(blckVals)
colnames(blckVals)
range(blckVals[,2:88], na.rm= T) #] max of 16 (from 4x4 aggregation of original raster)
hist(unlist((blckVals[,2:88])))
hist(unlist((blckValsProb[,2:88])))
range(blckVals$Group.1)

# remove ID col (group.1 - first column). Each row is a cell/site - convert to matrix.
x_in <- as.matrix(blckVals[, 2:ncol(blckVals)])
x_in_prob <- as.matrix(blckValsProb[, 2:ncol(blckVals)])
x_in_ct <- as.matrix(blckValsCT[, 2:ncol(blckValsCT)])

dim(x_in)
dim(x_in_prob)

x_in_prob[1:100, 1:10]
blckValsProb[1:100, 1:10]

# Rs distribution
hist(colSums(x_in, na.rm = T))
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


# alpha <- x_in

# dim(x_in)
# range(x_in)
# 
# ## Total available area by species
# # Rs <- colSums(rVals, na.rm = T)[-ncol(rVals)] # remove last ID column
# Rs <- colSums(x_in, na.rm = T)
# 
# hist(Rs)
# 
# median(Rs)
# range(Rs) # 0 184574
# # 1852
# pcRs <- 0.25 # set target as percentage of median area, for example
# t <- pcRs* median(Rs)
# 
# # or each species with % of range
# t <- floor(pcRs * Rs)
# 
# # or set areas
# t <- 1500
# 
# Rs_t <- Rs - t
# 
# ## Alpha irreplaceability
# ## Cases -
# ## If site does not have species, already is 0
# 
# ## If target is 0 (only if specific species targets)
# alpha[alpha, t ==0] <- 0
# 
# ## assign 1 where all cells are needed to complete target
# # alpha[, t >= Rs] # this filter columns for which target is larger than total area
# # [alpha[, t >= Rs] > 0] # this filters whole matrix for cell > 0
# alpha[, t >= Rs][alpha[, t >= Rs] > 0] <- 1
# range(alpha)
# 
# # for remainder...
# tmp <- t(t(alpha[, t < Rs]) / Rs_t[t < Rs])
# dim(tmp)
# range(tmp)
# tmp[tmp > 1] <- 1 # adjust back to 1 if some sites have more than in total - target
# 
# alpha[, t < Rs] <- tmp
# range(alpha)
# 
# # check matrix operations
# # m <- matrix(c(5,3,6,10,20,30,1,2,3), ncol = 3)
# # v <- c(10, 20, 30)
# # m
# # # m*v
# # t(t(m) * v) # multiply each row by v
# 
# 
# 
# ## Beta irreplaceability
# # rowwise
# beta <- apply(alpha, 1, function(x) 1 - prod(1- x))
# range(beta)
# 
# ## put back into raster
# beta.r <- raster(rBlock)
# 
# beta.r[] <- beta
# 
# plot(beta.r)

#### Function version

irrAB <- function(x, pc = NULL, type= c("total", "median"), tr = NULL, r = NULL, res = c("both", "alpha", "beta")){
  
  # x is a matrix of sites by species, with values of species pop/area, etc within sites
  # pc is % threshold of total area
  # tr is a vector of thresholds or single treshold.
  # r is an optional r template to write results to and return only beta raster
  # res value to return
  
  res <- match.arg(res)
  type <- match.arg(type)
  
  ## Total available area by species
  Rs <- colSums(x, na.rm = T) # range(Rs)
  
  # each species with % of range
  if(!is.null(pc)) {
    
    tr <- switch(type, 
                 total = pc * Rs,
                 median = pc * median(Rs))
  }
  Rs_t <- Rs - tr
  
  ## Alpha irreplaceability
  ## If site does not have species, already is 0
  ## If target is 0 (only if specific species targets)
  if(length(tr) == ncol(x)) x[, tr ==0] <- 0
  
  ## assign 1 where all cells are needed to complete target (where target is greater than total available )
  # x[, tr >= Rs] # this filter columns for which target is larger than total area
  # [x[, tr >= Rs] > 0] # this filters whole matrix for cell > 0
  x[, tr >= Rs][x[, tr >= Rs] > 0] <- 1
  
  # for remainder... 
  tmp <- t(t(x[, tr < Rs]) / Rs_t[tr < Rs]) # divide each row by Rs_t for those species with tr < Rs
  tmp[tmp > 1] <- 1 # adjust back to 1 if some sites have more than in total - target # range(tmp, na.rm =T)
  x[, tr < Rs] <- tmp
  
  ## Beta irreplaceability
  # rowwise
  beta <- apply(x, 1, function(z) 1 - prod(1- z)) # range(beta, na.rm = T)
  
  ## put back into raster
  if(!is.null(r)) {
    beta.r <- raster(r)
    beta.r[] <- beta
    return(beta.r)
  } else switch(res, 
                alpha = return(x),
                beta = return(beta),
                both =  return(list(alpha=x, beta=beta))
  )
  
}

## Binary maps
abs <- lapply(seq(.25, .95, .05), function(i) irrAB(x = x_in, pc = i, type = "median", res = "beta"))
abs <- lapply(seq(500,10000, 1000), function(i) irrAB(x = x_in, tr = i, type = "median", res = "beta"))
# median 1852, max 184574
par(mfrow = c(3,5), mar = c(2,2,2,2))
sapply(abs, hist)

abs <- lapply(seq(500,10000, 1000), function(i) irrAB(x = x_in, tr = i, type = "median", r = rBlock))
stck <- stack(abs)
names(stck) <- paste("treshold", seq(500,10000, 1000))
plot(stck)

abs <- lapply(seq(0.5,0.95, 0.05), function(i) irrAB(x = x_in, pc = i, type = "total", r = rBlock))
stck <- stack(abs)
names(stck) <- paste("pc", seq(0.5,0.95, 0.05))
plot(stck)


##Prob maps
abs <- lapply(seq(.25, 1, .05), function(i) irrAB(x = x_in_prob, pc = i, type = "total", res = "beta"))
abs <- lapply(seq(1000,10000, 1000), function(i) irrAB(x = x_in_prob, tr = i, res = "beta"))
# median  33609.72

par(mfrow = c(3,5), mar = c(2,2,2,2))
sapply(abs, hist)

pc.seq <- seq(.25, 1, .05)
abs.r <- lapply(pc.seq, function(i) irrAB(x = x_in_prob, pc = i, type = "total", r = rBlock))
tr.seq <- seq(500,10000, 1000)
abs.r <- lapply(tr.seq, function(i) irrAB(x = x_in_prob, tr = i, r = rBlock))

stck <- stack(abs.r)

names(stck) <- paste("tr", tr.seq)
names(stck) <- paste("pc", pc.seq)
plot(stck)

### CT Maps
abs <- lapply(seq(.25, 1, .05), function(i) irrAB(x = x_in_ct, pc = i, type = "total", res = "beta"))
abs <- lapply(seq(1000,10000, 1000), function(i) irrAB(x = x_in_ct, tr = i, res = "beta"))
# median  33609.72

par(mfrow = c(3,5), mar = c(2,2,2,2))
sapply(abs, hist)

pc.seq <- seq(.25, 1, .05)
abs.r <- lapply(pc.seq, function(i) irrAB(x = x_in_ct, pc = i, type = "total", r = rBlock))
tr.seq <- seq(500,10000, 1000)
abs.r <- lapply(tr.seq, function(i) irrAB(x = x_in_prob, tr = i, r = rBlock))

stck <- stack(abs.r)

names(stck) <- paste("tr", tr.seq)
names(stck) <- paste("pc", pc.seq)
plot(stck)





beta.r <- irrAB(x = x_in, tr = 25, type = "median", r = rBlock)

beta.r.bin <- irrAB(x = x_in, pc = 0.90, type = "total", r = rBlock)
beta.r.prob <- irrAB(x = x_in_prob, pc = 0.90, type = "total", r = rBlock)

tmp_stck <- stack(beta.r.bin, beta.r.prob)
names(tmp_stck) <- c("beta_binary_90pc", "beta_prob_90pc")
plot(tmp_stck)

### Summarise by logging plots

## Bring in logging and HJA data
## bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
hja <- st_read(file.path(gis_in, "shape/HJA_Boundary.shp"))
hja_bound <- subset(hja, FP_NAME == "H.J. Andrew Experimental Forest")
hja.utm <- st_transform(hja_bound, crs = utm10N)

## disturbance
cut.sf <- st_read(file.path(gis_in, "shape/disturbance.shp"))
cut.utm <- st_transform(cut.sf, crs = utm10N)
plot(st_geometry(cut.utm))
plot(hja.utm, add = T, col = NA, border = "blue")


## Edited study area
### bring in manually edited prediction area outline to replace above
aoi.pred.sf_edit <- st_read(file.path(gis_out, "s_utm/aoi_pred_sf_edit.shp"))
# aoi.pred.sf_edit <- st_make_valid(aoi.pred.sf_edit)

# clip cut to model area and convert to raster
# cut.r <- rasterize(cut.utm, r.msk, field = "YrsSinceDi")
# cut.r <- mask(cut.r, r.msk)
# plot(cut.r)

cut.ints <- cut.utm[lengths(st_intersects(cut.utm, aoi.pred.sf_edit)) > 0, ]

par(mfrow= c(1,1))
plot(beta.r.prob)
plot(cut.ints, col = NA, border = "grey60", add =T)
plot(hja.utm, add = T, border = "blue", col = NA)
# plot(cut.ints[cut.ints$TREATMENT_ == "Clearcut",], col = NA, border = "blue", add =T)

cut_stats <- extract(beta.r, cut.ints, fun = mean, na.rm = T)
cut_stats <- extract(beta.r, cut.ints, fun = max, na.rm = T)

head(cut_stats)

cut.ints$irr_mn <- cut_stats

plot(cut.ints$YrsSinceDi, cut.ints$irr_mn)

plot(cut.ints[, "irr_mn"])

hist(values(beta.r))
plot(beta.r)
plot(hja.utm, add = T, border = "blue", col = NA)
plot(cut.ints, col = NA, border = "grey", add =T)
# plot(cut.ints[cut.ints$TREATMENT_ == "Clearcut",], col = NA, border = "blue", add =T)

