#### Model set up #### 

## Christian D
## 12/11/2020

## Only local: 
#setwd("J:/UEA/gitRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada/results/RRR_test3")
setwd("~/oregon_ada/results/RRR_test3")
# dir()

library(Hmsc)

# subset of 50 species, biased to higher prevalence
load("data/allData.Rdata") # S,X,Y,P & Tax (where P is based on Tax)

# Check for absent (0) or ubiquitous species (1).
range(colMeans(Y.train > 0))
range(colSums(Y.train > 0))
hist(colSums(Y.train > 0))

# any sites with no species?
sum(rowSums(Y.train) == 0)

# hist(colMeans(Y.train > 0), main = "prevalence")

# Species richness across sites.
# hist(rowSums(Y.train))

names(X.train)

# 4 predictors in XFormula
nrs <- c("elevation_m","insideHJA", "lg_YrsDisturb", "oldGrowthIndex")
# check names
all(nrs %in% colnames(X.train))

XFormula <- as.formula(paste0("~ ", paste0(nrs, collapse = " + ")))
XFormula

# some other Remote Sensed predictors into RRR
rrr.pred <- colnames(X.train)[grepl("mean|l_", colnames(X.train))]

XRRRFormula <- as.formula(paste0("~ ", paste0(rrr.pred, collapse = " + "), "-1"))
XRRRFormula

rm(rrr.pred, nrs)

# Note: X covariates are scaled in hmsc by default.

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
head(xy); tail(xy)

# d <- dist(xy)
# min(d)
# hist(d)

rL.site <- HmscRandomLevel(sData = xy)
str(rL.site)

sum(apply(Y.train, 2, max) > 1) # Y is presence absence only data

## Make models

models <- list(
  
  # basic spatial model, check everything works.
  m1 = Hmsc(Y = Y.train, 
            XData = X.train, XFormula = XFormula,
            phyloTree = P,
            distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            }),
  
  # RRR model, no spatial
  m2 = Hmsc(Y = Y.train,
            XData = X.train, XFormula = XFormula,
            XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 1, # try one axis
            phyloTree = P,
            distr = "probit"),
  
  # RRR model with spatial
  m3 = Hmsc(Y = Y.train,
            XData = X.train, XFormula = XFormula,
            XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 1,
            phyloTree = P,
            distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            }),
  
  # RRR model, no spatial, # 2 axes probably better, looking at ordinationa
  m4 = Hmsc(Y = Y.train,
            XData = X.train, XFormula = XFormula,
            XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 2, 
            phyloTree = P,
            distr = "probit"),
  
  # RRR model, no spatial, # 2 axes probably better, looking at ordinationa
  m5 = Hmsc(Y = Y.train,
            XData = X.train, XFormula = XFormula,
            XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 2, 
            phyloTree = P,
            distr = "probit",
            studyDesign = studyDesign,
            ranLevels = {
              list("site" = rL.site)
            })
  
)


names(models) <- c("pa_sp_simple", "pa_nsp_RRR_nc1", "pa_sp_RRR_nc1", "pa_nsp_RRR_nc2", "pa_sp_RRR_nc2")
modelnames <- c("pa_sp_simple", "pa_nsp_RRR_nc1", "pa_sp_RRR_nc1", "pa_nsp_RRR_nc2", "pa_sp_RRR_nc2")

models

rm(S.train, X.train, SXY.train, xy, Y.train, P, XFormula, XRRRFormula, rL.site, studyDesign)
save(models, modelnames, file = file.path("models", "unfitted_models.rdata"))

