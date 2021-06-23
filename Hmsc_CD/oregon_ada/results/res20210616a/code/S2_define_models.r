#### Model set up #### 

## Christian D
## 12/11/2020

## NEW TOPO VARS 14 Jan

## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
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
range(colMeans(Y.train.pa > 0))
min(colSums(Y.train.pa > 0))
# =5.

# Data reduction - reduce species - done in read data
# raretaxa <- which(colSums(Y.train.pa > 0) < 10)
# length(raretaxa)
# Y.train.pa2 <- Y.train.pa[, -raretaxa] # reduced species
# rm(raretaxa)

ncol(Y.train.pa) # full species

# any sites with no species?
sum(rowSums(Y.train.pa) == 0)

# Species richness across sites.
# hist(rowSums(Y.train.pa))

names(X.train)

XFormula <- as.formula(paste0("~ ", paste0(colnames(X.train), collapse = " + ")))
XFormula
# 

# Note: X covariates are scaled in hmsc by default.
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
  
  m1 = Hmsc(Y = Y.train.pa, # 
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit"),
  
  m2 = Hmsc(Y = Y.train.qp,
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit"),
  
  m3 = Hmsc(Y = Y.train.pa,
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            }),
  
  m4 = Hmsc(Y = Y.train.qp,
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            })
  
  )

names(models) <- c("pa_nsp_vars8t10", "qp_nsp_vars8t10", "pa_sp_vars8t10", "qp_sp_vars8t10")
modelnames <-  c("pa_nsp_vars8t10", "qp_nsp_vars8t10", "pa_sp_vars8t10", "qp_sp_vars8t10")

models

rm(S.train, X.train, P, rL.site, studyDesign, XFormula)

## now in pipeline
# save(models, modelnames, resFolder, file = file.path(modFolder, "unfitted_models.rdata"))

