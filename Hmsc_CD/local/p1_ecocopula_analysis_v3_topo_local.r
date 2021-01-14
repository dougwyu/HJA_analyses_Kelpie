
#### ecocopula analysis - PLOTS

library(mvabund)
library(ecoCopula)
library(doParallel)
library(foreach)

getwd()
wd <- here::here()
setwd(wd)

# get data
source("Hmsc_CD/local/L1_read_data.r")
# otu.pa.csv is Y.train.pa; otu.qp.csv is Y.train.qp
rm(Y.train.qp, Y.train.pa, X.train)

# load new topo vars
load("Hmsc_CD/oregon_ada/data/topo_data.rdata")

# remove NA
sum(is.na(topo.df))
ind <- which(is.na(topo.df$twi))
topo.df <- topo.df[-ind,]

# same with otu
otu.pa.csv <- otu.pa.csv[-ind,]

# same with S.train
S.train <- S.train[-ind,]

## Make list of formulae to loop through, 
colnames(topo.df)
c("be10", "slope", "Ess", "twi", "ht", "cov2_4", "cov4_16", "mTopo") # , "cut.r1k.pt"
c("be10", "slope", "Ess", "twi", "ht.r1k", "cov2_4.r1k", "cov4_16.r1k", "mTopo") # , "cut.r1k.pt"


fm <- list(as.formula(otu.pa ~ be10 + Nss + Ess + ht + cov2_4 +  cov4_16 + mTopo + ht.r1k + cov2_4.r1k + cov4_16.r1k)
           ,as.formula(otu.pa ~ be10 + slope + ht.r1k + cov2_4.r1k + cov4_16.r1k + mTopo + cut.r1k.pt)
           )
fm

## ALternative data specification as list with both response and predictors as list elements 
# (foreach wasn't finding otu.pa)
dataN <- c(list(otu.pa = mvabund(otu.pa.csv)), as.list(topo.df))
str(dataN)


# check all terms are in data object
unique(unlist(sapply(fm, all.vars)))
all(unique(unlist(sapply(fm, all.vars))) %in% names(dataN))

## set up parallel
# nCores <- 16
nCores <- 2

cl <- makeCluster(nCores)
registerDoParallel(cl)

# nBoot <- 10

# i = 1

modList <- foreach(i = seq_along(fm),
                   .combine = list,
                   .multicombine = TRUE, # ,.export = c("otu.pa", "env.vars", "S.train")
                   .noexport = c("otu.pa.csv","otu.qp.csv", "P")
                   ) %dopar% {
                     
    # do glm model
    mod <- mvabund::manyglm(fm[[i]], family = "negative.binomial", data = dataN)
    
    # do ordination
    mod.ord <- ecoCopula::cord(mod)
    
    # Site scores: join ordination axes to data frame of predictors, coordinates
    site_res <- data.frame(mod.ord$scores, topo.df, S.train[,c("UTM_E", "UTM_N")])
    # Species scores
    sp_res <- data.frame(mod.ord$loadings, species = colnames(dataN$otu.pa))
  
    # don't do this locally
    # do anova
    # summ <- mvabund::summary.manyglm(mod, nBoot = nBoot, test = "LR")
    # summ
    
    # list(mod = mod, ord = mod.ord, site = site_res, sp=sp_res, summ = summ)
    list(mod = mod, ord = mod.ord, site = site_res, sp=sp_res)
    
                   }

stopCluster(cl)

cor.preds <- colnames(topo.df)[sapply(topo.df, is.numeric)]

str(modList, max.level =2)
save(modList, fm, cor.preds, file = "Hmsc_CD/local/ecocopula_modList_topo.rdata")

