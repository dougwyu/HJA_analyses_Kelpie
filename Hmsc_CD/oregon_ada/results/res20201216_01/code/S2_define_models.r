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
# raretaxa <- which(colSums(Y.train > 0) < 5)
# length(raretaxa)
# Y.train <- Y.train[, -raretaxa]
# rm(raretaxa)

# any sites with no species?
sum(rowSums(Y.train) == 0)

# Species richness across sites.
# hist(rowSums(Y.train))

# names(X.train)

# GIS XFormula set up for the pilot model, plus sampling session:.
# elevation, canopy.ht, yrs.disturb.min

# non.rs <- colnames(X.train)[!grepl("^l_|^mean", colnames(X.train))]
# pairs(X.train[,non.rs])

# non rs,  
preds <- colnames(X.train)
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
length(unique(S.train$uniqueID))

studyDesign <- data.frame(
  site = as.factor(S.train$SiteName),
  unique_sample = as.factor(S.train$uniqueID)
)


xy <- data.frame(S.train[match(unique(S.train$SiteName), S.train$SiteName), c("UTM_E", "UTM_N")])

# check for duplicated coordiantes>
sum(duplicated(xy))
rownames(xy) <- unique(S.train$SiteName)

# head(xy); tail(xy)

rL.site <- HmscRandomLevel(sData = xy)

## Make models
models <- list(
  
  # m1 = Hmsc(Y = Y.train,
  #                        XData = X.train, XFormula = XFormula,
  #                        phyloTree = P,
  #                        distr = "probit",
  #                        studyDesign = studyDesign,
  #                        ranLevels = {
  #                          list("site" = rL.site)
  #                        })
  m2 = Hmsc(Y = Y.train,
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit")
  )

# names(models) <- c("pa_sp_vif", "pa_nsp_vif")
# modelnames <- c("pa_sp_vif", "pa_nsp_vif")

names(models) <- c("pa_nsp_vif")
modelnames <- c("pa_nsp_vif")


models

rm(S.train, X.train, P, XFormula, rL.site, studyDesign, preds)
## now in pipeline
# save(models, modelnames, resFolder, file = file.path(modFolder, "unfitted_models.rdata"))

