#### Model set up #### 

## Christian D
## 12/11/2020

## Only local: 
# setwd("J:/UEA/Oregon")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/Oregon_winscp")
dir()


library(Hmsc)

load(file = file.path("data", "allData.Rdata")) # S,X,Y,P & Tax (where P is based on Tax)

# Y = 518 Forest Malaise-trapped invertebrate taxa (quasi-probability species)
# with a DNA spike-in standard used to estimate how a species' biomass
# changes from sample to sample).

# Check for absent (0) or ubiquitous species (1).
range(colMeans(Y.train > 0))
min(colSums(Y.train > 0))
# =5.

# Data reduction for the simple pilot model, included species with at least
# observations only.

raretaxa <- which(colSums(Y.train > 0) < 10)
length(raretaxa)

Y.train <- Y.train[, -raretaxa]

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

non.rs <- colnames(X.train)[!grepl("^l_|^mean", colnames(X.train))]
pairs(X.train[,non.rs])

nrs <- c("elevation","canopy.ht","metre.road","metre.stream","yrs.disturb.min","hja") 
# pairs(X.train[,nrs], lower.panel = panel.cor)

# hja - factor binary
# log metre.road and yrs.disturb.min

# plot(X.train[, c("elevation", "canopy.ht", "yrs.disturb.min")])
cor(X.train[, c("elevation", "canopy.ht", "metre.stream")])

XFormula <- ~ elevation + canopy.ht + metre.stream

rs_vars <- colnames(X.train)[grepl("^l_|^mean", colnames(X.train))]
# rs_vars
l_vars <- colnames(X.train)[grepl("^l_.*", colnames(X.train), perl = T)]
l_vars
vi_vars <- colnames(X.train)[grepl("^mean", colnames(X.train))]
vi_vars

# summary(X.train[,rs_vars])
# pairs(X.train[,rs_vars], diag.panel = panel.hist, lower.panel = panel.cor)
vi <- c("mean.EVI", "mean.bright", "mean.wet") # uncorrelated mainly

# pairs(X.train[,l_vars[!grepl("all", l_vars)]], diag.panel = panel.hist, lower.panel = panel.cor)
lv <- c("l_Cover_2m_4m","l_Cover_2m_max","l_Cover_4m_16m","l_rumple")

# pairs(X.train[,lv], diag.panel = panel.hist, lower.panel = panel.cor)
# pairs(X.train[,vi], diag.panel = panel.hist, lower.panel = panel.cor)

XRRRFormula <- as.formula(paste0("~ ", paste0(vi, collapse = " + "), "-1"))
XRRRFormula

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

xy <- data.frame(S.train[match(unique(S.train$SiteName), S.train$SiteName), c("Route_x", "Route_y")]) 
# check for duplicated coordiantes>
sum(duplicated(xy))
rownames(xy) <- unique(S.train$SiteName)
head(xy); tail(xy)

d <- dist(xy)
((87*87) - 87)/2
min(d)
# hist(d)

rL.site <- HmscRandomLevel(sData = xy)
str(rL.site)

sum(apply(Y.train, 2, max) > 1) # Y is presence absence only data

m1 <- Hmsc(
  Y = Y.train, 
  XData = X.train, XFormula = XFormula,
  phyloTree = P,
  distr = "probit"
  )

m2 <- Hmsc(
  Y = Y.train, 
  XData = X.train, XFormula = XFormula,
  phyloTree = P,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = {
    list("site" = rL.site)
  }
)

m3 <- Hmsc(
  Y = Y.train, 
  XData = X.train, XFormula = XFormula,
  XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 1, # try 1 to start with
  phyloTree = P,
  distr = "probit"
  )

m4 <- Hmsc(
  Y = Y.train, 
  XData = X.train, XFormula = XFormula,
  XRRRData = X.train, XRRRFormula = XRRRFormula, ncRRR = 1, # try 1 to start with
  phyloTree = P,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = {
    list("site" = rL.site)
  }
)


m1; m2; m3; m4

models <- list(m1, m2, m3, m4)
names(models) <- c("pa_simp", "pa_simp_sp", "pa_RRR_simp", "pa_RRR_simp_sp")
modelnames <- c("pa_simp", "pa_simp_sp", "pa_RRR_simp", "pa_RRR_simp_sp")

models

save(models, modelnames, file = file.path("models", "unfitted_models.rdata"))

