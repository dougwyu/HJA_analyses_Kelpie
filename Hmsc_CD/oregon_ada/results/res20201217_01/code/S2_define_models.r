#### Model set up #### 

## Christian D
## 12/11/2020

## Only local: 
# setwd("J:/UEA/gitRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("~/Oregon_winscp")
# dir()
# rm(list = ls())

library(dplyr)
library(Hmsc)
# load(file = file.path("data", "allData_vif.Rdata")) # S,X,Y,P & Tax (where P is based on Tax)


# Check for absent (0) or ubiquitous species (1).
range(colMeans(Y.train.pa > 0))
min(colSums(Y.train.pa > 0))

range(colMeans(Y.train.qp))
min(colSums(Y.train.qp > 0))
# =5.

# Data reduction for the simple pilot model
# raretaxa <- which(colSums(Y.train > 0) < 5)
# length(raretaxa)
# Y.train <- Y.train[, -raretaxa]
# rm(raretaxa)

# any sites with no species?
sum(rowSums(Y.train.pa) == 0)
sum(rowSums(Y.train.qp) == 0)

## reduce number of species to those with 
raretaxa <- which(colSums(Y.train.pa > 0) < 8) # 106 spp removed with <8
length(raretaxa)

Y.train.pa <- Y.train.pa[, -raretaxa]
Y.train.qp <- Y.train.qp[, -raretaxa]
rm(raretaxa)

# Log qp data for gaussian family
Y.train.qp <- log(Y.train.qp+0.001) # avoid -Inf in Y data. 

# Species richness across sites.
# hist(rowSums(Y.train))

# names(X.train)

# GIS XFormula set up for the pilot model, plus sampling session:.
# elevation, canopy.ht, yrs.disturb.min

# non.rs <- colnames(X.train)[!grepl("^l_|^mean", colnames(X.train))]
# pairs(X.train[,non.rs])

# non rs,  
preds <- colnames(X.train)
preds

# X.train %>%
#   select(where(is.numeric)) %>%
#   cor() %>%
#   corrplot::corrplot(method = "ellipse", type = "lower")
#   
# X.train %>%
#   select(where(is.numeric)) %>%
#   cor()%>%
#   abs() > 0.5
  
preds <- preds[!preds %in% c("precipitation", "mean.green", "lg_cover2m_4m", "lg_cover2m_max")]
preds

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
  
  m1 = Hmsc(Y = Y.train.pa,
            XData = X.train, XFormula = XFormula,
            phyloTree = P,distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            }
            ),
  
  m2 = Hmsc(Y = Y.train.qp, YScale = TRUE,
            XData = X.train, XFormula = XFormula,
            phyloTree = P, distr = "normal",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            }
            )
)

names(models) <- c("pa_sp_vif", "qp_sp_vif")
modelnames <- c("pa_sp_vif", "qp_sp_vif")

models

rm(S.train, X.train, P, XFormula, rL.site, studyDesign, preds, Y.train.pa, Y.train.qp, xy)
## now in pipeline
# save(models, modelnames, resFolder, file = file.path(modFolder, "unfitted_models.rdata"))

