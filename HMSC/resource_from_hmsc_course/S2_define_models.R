setwd("P:/h572/hmsc_course/li_yuanheng/data")
# OR ON A MAC:
# setwd("/Volumes/group/h572/hmsc_course/li_yuanheng")
localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)

load(file=file.path(DataDir,"allData.R")) #S,X,Y,P & Tax (where P is based on Tax)

# Y = 518 Forest Malaise-trapped invertebrate taxa (quasi-probability species) 
# with a DNA spike-in standard used to estimate how a species' biomass 
# changes from sample to sample).

# Check for absent (0) or ubiquitous species (1).
range(colMeans(Y>0))

min(colSums(Y>0))
# =5.

# Data reduction for the simple pilot model, included species with at least
# observations only.

raretaxa = which(colSums(Y>0)<25)
length(raretaxa)
# 430 of the original 518 taxa occur at < 25 sites. 
# Removed from pilot model so that it runs fast on the course.

Y = Y[,-raretaxa]
# Excluding rare taxa leaves 88 taxa in the analysis.

hist(colMeans(Y>0),main="prevalence")
hist(log(Y[Y>0]),main="log abundance conditional on presence")

# Many species are rare (absent in many samples) - need zero-inflated model. Choice for the
# pilot model is a hurdle model: species presence-absence and log(abundance) separately.

table(rowSums(Y>0))
# Species richness across sites.
hist(rowSums(Y>0))

names(X)

# Mirkka: At Otso's suggestion, I have moved the S variable session into the XData, in order
# to include this as a fixed effect representing possible changes in community structure
# resulting from the effect of the regional fire.

X$session = as.factor(S$session)

# GIS XFormula set up for the pilot model, plus sampling session:.
# elevation, canopy.ht, yrs.disturb.min

X = data.frame(X[, c("elevation" , "canopy.ht" , "yrs.disturb.min", "session")])

plot(X[, c("elevation" , "canopy.ht" , "yrs.disturb.min")])
cor(X[, c("elevation" , "canopy.ht" , "yrs.disturb.min")])

XFormula = ~elevation + canopy.ht + yrs.disturb.min + session

# Note: X covariates are scaled in hmsc by default.

head(S)
# 94 sites with unique spatial coordinates, sampled during two time periods
# (sessions, included in the XFormula) and with either one or two Malaise traps (M1, M2).

# S$site_trap_period is a unique sample code in the StudyDesign matrix, but is currently not included in the models.

studyDesign = data.frame(site = as.factor(S$SiteName), unique_sample = as.factor(S$site_trap_period))

xy = data.frame(S[match(unique(S$SiteName), S$SiteName), 2:3])
rownames(xy) = unique(S$SiteName)
rL.site = HmscRandomLevel(sData = xy)

Ypa = 1*(Y>0)
Yabu = Y
Yabu[Y==0] = NA
Yabu=log(Yabu)

m1 = Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
          phyloTree = P,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site)})

m2 = Hmsc(Y=Yabu, YScale = TRUE,
          XData = X,  XFormula = XFormula,
          phyloTree = P,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site)})

models = list(m1,m2)
modelnames = c("presence_absence","abundance_COP")

save(models,modelnames,file = file.path(ModelDir, "unfitted_models"))
