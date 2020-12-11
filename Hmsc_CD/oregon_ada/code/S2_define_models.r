#### Model set up #### 

## Christian D
## 12/11/2020

## Only local: 
# setwd("J:/UEA/Oregon/Oproject/oregon_ada")
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

## Make data sets with  highest species prevalence
sp.prev <- colSums(Y.train)/nrow(Y.train)
hist(sp.prev)

table(cut(sp.prev, quantile(sp.prev), include.lowest = T, right = F), useNA= "always")

# make groups
qGrps <- cut(sp.prev, quantile(sp.prev), include.lowest = T, right = F, labels = 1:4)
table(qGrps, useNA = "always")

Y.q4 <- Y.train[,qGrps == 4]
Y.q3 <- Y.train[,qGrps == 3]

set.seed(2666)
Y.qr1 <- Y.train[,sample(seq_along(Y.train), size = sum(qGrps == 4))]
Y.qr2 <- Y.train[,sample(seq_along(Y.train), size = sum(qGrps == 3))]

rm(Y.train) # make sure it's not anywhere below

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

XFormula <- as.formula(paste0("~ ", paste0(preds, collapse = " + ")))
XFormula

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


xy <- data.frame(S.train[match(unique(S.train$SiteName), S.train$SiteName), c("UTM_E", "UTM_N")])

# check for duplicated coordiantes>
sum(duplicated(xy))
rownames(xy) <- unique(S.train$SiteName)

# head(xy); tail(xy)

# d <- dist(xy)
# ((87*87) - 87)/2
# min(d)
# hist(d)

rL.site <- HmscRandomLevel(sData = xy)

# str(rL.site)

## Make models

models <- list(m1 = Hmsc(Y = Y.q4,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign,
                         ranLevels = {
                           list("site" = rL.site)
                         }),
               
               m2 = Hmsc(Y = Y.q3,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign,
                         ranLevels = {
                           list("site" = rL.site)
                         }),
               
               m2 = Hmsc(Y = Y.qr1,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign,
                         ranLevels = {
                           list("site" = rL.site)
                         }),
               
               m2 = Hmsc(Y = Y.qr2,
                         XData = X.train, XFormula = XFormula,
                         phyloTree = P,
                         distr = "probit",
                         studyDesign = studyDesign,
                         ranLevels = {
                           list("site" = rL.site)
                         })
               )

names(models) <- c("pa_sp_tp8_q4", "pa_sp_tp8_q3", "pa_sp_tp8_qr1", "pa_sp_tp8_qr2")
modelnames <- c("pa_sp_tp8_q4", "pa_sp_tp8_q3", "pa_sp_tp8_qr1", "pa_sp_tp8_qr2")

models

rm(S.train, X.train, SXY.train, xy, P, XFormula, rL.site, studyDesign, preds, Y.q4, Y.q3, Y.qr1, Y.qr2, qGrps, sp.prev)
## now in pipeline
# save(models, modelnames, resFolder, file = file.path(modFolder, "unfitted_models.rdata"))

