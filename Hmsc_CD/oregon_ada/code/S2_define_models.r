#### Model set up #### 

## Christian D
## 12/11/2020

## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("~/Oregon_winscp")
# dir()
# rm(list = ls())

library(Hmsc)

# load(file = file.path("data", "allData.Rdata")) # S,X,Y,P & Tax (where P is based on Tax)

# Check for absent (0) or ubiquitous species (1).
range(colMeans(Y.train > 0))
min(colSums(Y.train > 0))
# =5.

# Data reduction for the simple pilot model
raretaxa <- which(colSums(Y.train > 0) < 10)
length(raretaxa)
Y.train <- Y.train[, -raretaxa]
rm(raretaxa)

# any sites with no species?
sum(rowSums(Y.train) == 0)

# hist(colMeans(Y.train > 0), main = "prevalence")

# Species richness across sites.
# hist(rowSums(Y.train))

names(X.train)

# Mirkka: At Otso's suggestion, I have moved the S variable session into the XData, in order to include this as a fixed effect representing possible changes in community structure resulting from the effect of the regional fire.

# ONly single session now in this data.

# GIS XFormula set up for the pilot model, plus sampling session:.
# elevation, canopy.ht, yrs.disturb.min

# non.rs <- colnames(X.train)[!grepl("^l_|^mean", colnames(X.train))]
# pairs(X.train[,non.rs])

# non rs,  
preds <- c("elevation_m", "insideHJA", "mean.wet", "lg_YrsDisturb", "mean.EVI", "lg_cover2m_max","lg_cover2m_4m", "l_rumple")
# top 8 univariate/jacknife predictors, without correlated (B1, B4, precipitation, NDVI, green, l_p25)
# check names
# all(preds %in% colnames(X.train))
# pairs(X.train[,preds])
# hist(X.train[, "mean.EVI"])
# plot(X.train[, "mean.EVI"], X.train[, "mean.NDVI"])
ind <- which(X.train[,"mean.NDVI"] < 0.8)
# 45 less than 0.7
X.train_adj_evi <- X.train[-ind,]
S.train_adj_evi <- S.train[-ind,]
Y.train_adj_evi <- Y.train[-ind,]

XFormula <- as.formula(paste0("~ ", paste0(preds, collapse = " + ")))
XFormula

rm(preds, ind)

# Note: X covariates are scaled in hmsc by default.
# 94 sites with unique spatial coordinates, sampled during two time periods
# and with Malaise traps M1

# S$site_trap_period is a unique sample code in the StudyDesign matrix, but is currently not included in the models.
head(S.train)
length(unique(S.train$SiteName))
length(unique(S.train$site_trap_period))

studyDesign <- data.frame(
  site = as.factor(S.train$SiteName),
  unique_sample = as.factor(S.train$site_trap_period)
)

studyDesign_adj <- data.frame(
  site = as.factor(S.train_adj_evi$SiteName),
  unique_sample = as.factor(S.train_adj_evi$site_trap_period)
)


xy <- data.frame(S.train[match(unique(S.train$SiteName), S.train$SiteName), c("UTM_E", "UTM_N")])
xy_adj <- data.frame(S.train_adj_evi[match(unique(S.train_adj_evi$SiteName), S.train_adj_evi$SiteName), c("UTM_E", "UTM_N")]) 

# check for duplicated coordiantes>
sum(duplicated(xy))
sum(duplicated(xy_adj))

rownames(xy) <- unique(S.train$SiteName)
rownames(xy_adj) <- unique(S.train_adj_evi$SiteName)

# head(xy); tail(xy)

# d <- dist(xy)
# ((87*87) - 87)/2
# min(d)
# hist(d)


rL.site <- HmscRandomLevel(sData = xy)
rL.site_adj <- HmscRandomLevel(sData = xy_adj)

# str(rL.site)

sum(apply(Y.train, 2, max) > 1) # Y is presence absence only data

## Make models

models <- list(m1 = Hmsc(Y = Y.train,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit"),
               
               m2 = Hmsc(Y = Y.train,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign,
                         ranLevels = {
                           list("site" = rL.site)
                         }),
               
               m3 = Hmsc(Y = Y.train_adj_evi,
                         XData = X.train_adj_evi, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit"),
               
               m4 = Hmsc(Y = Y.train_adj_evi,
                         XData = X.train_adj_evi, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign_adj,
                         ranLevels = {
                           list("site" = rL.site_adj)
                         })
               )

names(models) <- c("pa_nsp_tp8", "pa_sp_tp8", "pa_nsp_tp8_adj_evi", "pa_sp_tp8_adj_evi")
modelnames <- c("pa_nsp_tp8", "pa_sp_tp8", "pa_nsp_tp8_adj_evi", "pa_sp_tp8_adj_evi")

models

rm(S.train, X.train, SXY.train, xy, Y.train, P, XFormula, rL.site, studyDesign, studyDesign_adj, rL.site_adj,
   S.train_adj_evi, X.train_adj_evi, Y.train_adj_evi, xy_adj)
## now in pipeline
# save(models, modelnames, resFolder, file = file.path(modFolder, "unfitted_models.rdata"))

