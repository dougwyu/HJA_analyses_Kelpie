# Jan 18, 2021

## Data set up for models on ADA

# code on sjsdm model run with DNN (for ADA cluster) - from Yuanheng
# to test how long it would take to run one model on the cluster

## I run this locally, save the data, copy it to ADA, then load the saved data file on ADA - Christian
## could also just take the minimap file straight from github online and then follow this script in ADA.

wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"  local
setwd(wd)

# lapply(c('tidyverse','reticulate','sjSDM','glue'), library, character.only=T) # removed here, conflicted

library(glue)
library(dplyr)

# conflict_prefer("filter", "dplyr") # not present
# conflict_prefer("select", "dplyr")# written function in full
# conflict_prefer('colSums', 'base') # not present

# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"

minocc = 5
abund = 'pa' # pa , qp    # !!! change accordingly
trap <- "M1"; period = "S1"
date.model.run = 20210119   # !!! change accordingly

outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('Kelpie_maps/{outputidxstatstabulatefolder}')

sjsdmV = '0.1.3.9000' # package version

# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')

# ..... load data ......

# what file am I loading?
basename(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')))

# when was it last modified?
file.mtime(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')))
# "2021-02-02 11:11:23 GMT"

alldata = read.csv(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
dim(alldata)
names(alldata)[1:10]

# ..... select trap, session .....
trap; period

alldata1 = subset(alldata, trap == 'M1' & period == 'S1')
dim(alldata1)

# select species with minocc
a = alldata1 %>% dplyr::select(contains('__'))
a = a[,which(vegan::specnumber(a, MARGIN=2)>=minocc)]
dim(a)

# join species back to data
alldata1 = cbind(dplyr::select(alldata1, -contains('__')), a)
dim(alldata1)

## Log the 
# lg_DistStream = log(distToStream_m + 0.001),
# lg_DistRoad = log(distToRoad_m + 0.001),
# lg_YrsDisturb = log(YrsSinceDist + 0.001),

alldata1$lg_DistStream <- log(alldata1$distToStream_m + 0.001)
alldata1$lg_DistRoad <- log(alldata1$distToRoad_m + 0.001)
alldata1$lg_YrsDisturb <- log(alldata1$YrsSinceDist + 0.001)

# check data 
# non species names
non.sp.cols <- colnames(alldata1)[!grepl("__", colnames(alldata1))]
395-268
non.sp.cols
head(alldata1[,non.sp.cols])
summary((alldata1[,non.sp.cols]))
# data is not scaled here.. 

# ... 80% of data as training ...
# random selection
num.sample = dim(a)[1]
select.percent = .8
ssdd = 100 		# please keep this value to make results comparable!
set.seed(ssdd)
a = base::sample(1:num.sample, round(num.sample*select.percent))
a
# [1] 74 78 23 70  4 55 85  7 81 83 43 61 12 51 72 18 25  2 75 68 69 52 48 32 21
# [26] 27 39 57 16 11 67 71  6 29 45 30 53 79 86 31 33 49 82 28 47 41 87 42 24 80
# [51]  1  9 20 14 35 40  3 34 84 19 46 63 44 36 26  5 15 22 58 76
otu = dplyr::select(alldata1, contains('__'))
otu.train = otu[a,]
dim(otu.train)
otu.test = otu[-a,]

# old vars
# envnames = c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")

# new vars
cat(paste(non.sp.cols, collapse = '", "'))
newvars <- c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")

## env vars
# all.env.vars <- alldata1[,newvars]
# save(all.env.vars, otu, file = "Hmsc_CD/working_cd/checkData.rdata")

# divide into training and testing
ori.env.train = alldata1[a,newvars]
ori.env.test = alldata1[-a,newvars]

str(ori.env.test)

ori.XY = dplyr::select(alldata1, starts_with('UTM'))
ori.XY.train = ori.XY[a,]
ori.XY.test = ori.XY[-a,]
str(ori.XY.train)

# ... view data ...
par(mfrow=c(1,2))
hist(ori.env.train$be10)
hist(alldata1$be10)

library(ggplot2)
ggplot(ori.XY, aes(UTM_E, UTM_N)) + geom_point() + geom_point(data=ori.XY.train, aes(colour='red')) + scale_colour_manual(labels = c('training'), values = c("red"))

# 0/1
if (abund == 'pa')
{
  otu.train = as.data.frame((otu.train>0)*1)
  otu.test = as.data.frame((otu.test>0)*1)
}

sum(colSums(otu.train) < minocc)
sum(colSums(otu.train) == 0)


# .. env data
colnames(ori.env.train)
scale.env.train.all = dplyr::select(ori.env.train, -'insideHJA') %>% scale() 
str(scale.env.train.all)

scale.env.train = data.frame(scale.env.train.all) %>% tibble::add_column(insideHJA=as.factor(ori.env.train$insideHJA), .before=names(ori.env.train)[2])
str(scale.env.train)

# data frame of means and sd for scaling
dd.env.scaler = data.frame(t(data.frame(env.mean = attr(scale.env.train.all, "scaled:center"), env.sd = attr(scale.env.train.all, "scaled:scale"))))
str(dd.env.scaler)

rm(scale.env.train.all)

scale.env.test = as.data.frame(do.call(rbind, apply(dplyr::select(ori.env.test, -'insideHJA'), 1, function(x){(x-dd.env.scaler['env.mean',])/dd.env.scaler['env.sd',]} ) )) %>% tibble::add_column(insideHJA=as.factor(ori.env.test$insideHJA), .before=names(ori.env.test)[2])
str(scale.env.test)

# .. spatial data
XY.train.all = scale(ori.XY.train)
str(XY.train.all)

XY.train = data.frame(XY.train.all)
str(XY.train)

dd.xy.scaler = data.frame(t(data.frame(sp.mean = attr(XY.train.all, "scaled:center"), sp.sd = attr(XY.train.all, "scaled:scale"))))
str(dd.xy.scaler)
base::rownames(dd.xy.scaler)

rm(XY.train.all)

XY.test = as.data.frame(do.call(rbind, apply(ori.XY.test, 1, function(x){(x-dd.xy.scaler['sp.mean',])/dd.xy.scaler['sp.sd',]} ) ))
str(XY.test)

# ... view data ...
par(mfrow=c(1,2))
hist(ori.env.train[,2],breaks = 10)
hist(ori.env.test[,2], breaks = 10)

hist(XY.test[,2])
hist(XY.train[,2])

rm(dd.env.scaler, dd.xy.scaler)

# make sure input of sjsdm are numeric matrix
s.otu.train = as.matrix(otu.train)
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)

s.otu.test = as.matrix(otu.test)
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)


# save model data
save(s.otu.train,scale.env.train, XY.train,  s.otu.test, scale.env.test, XY.test, abund, 
     ori.XY.train, ori.XY.test, file = paste0("Hmsc_CD/oregon_ada/data/yuanheng_mod_data_newVars_", abund, ".rdata"))

# # one model to see how long it takes in ADA
# load("data/yuanghen_mod_data.rdata")
# # s.otu.train,scale.env.train, XY.train,  s.otu.test, scale.env.test, XY.test
# 
# 
# # set variables
# formula.env = 'envDNN'
# lambda.env = .1
# alpha.env = .5
# lambda.sp = .1
# alpha.sp = .9 
# hidden = c(50L,50L,10L)
# acti.sp = 'relu'
# drop = .3
# 
# model.train = sjSDM(Y = s.otu.train,
#                     env = DNN(data=scale.env.train, formula = ~.,
#                               hidden=hidden, lambda = lambda.env, alpha = alpha.env, activation=acti.sp, dropout=drop, bias = T),
#                     biotic = bioticStruct(lambda=lambda.env, alpha=alpha.env, on_diag=F, inverse = FALSE),
#                     spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp, alpha=alpha.sp),
#                     learning_rate = 0.003, # 0.003 recommended for high species number 
#                     step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L 
# )
# 
# names(model.train)
# 
# # saveRDS(model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS',  glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{date.model.run}_{formula.env}.RDS')) )
# 
# 
# 
