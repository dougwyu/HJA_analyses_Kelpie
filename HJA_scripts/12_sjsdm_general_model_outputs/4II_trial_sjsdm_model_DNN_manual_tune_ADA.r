# Jan 18, 2021
# code on sjsdm model run with DNN (for ADA cluster)
# to test how long it would take to run one model on the cluster



```{r setup}
# setwd('.../HJA_analyses_Kelpie/HJA_scripts/12_sjsdm_general_model_outputs')
	
lapply(c('tidyverse','here','conflicted','reticulate','sjSDM','glue'), library, character.only=T)
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.3.9000’
	

```


```{r set-names}
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
minocc = 5
abund = 'qp' # pa , qp    # !!! change accordingly
trap <- "M1"; period = "S1"
date.model.run = 20210119   # !!! change accordingly
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
```


```{r load-data}
# ..... load data ......
alldata = read.csv(here(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
dim(alldata)
names(alldata)[1:10]
# 103
	
# roughness index
# from "Hmsc_CD/oregon_ada/data/demStats.rdata"
load(here('source','demStats.rdata'))
str(dem_stats)
	

```


```{r select subsets}
# ..... select trap, session .....
trap; period
	
alldata1 = subset(alldata, trap == 'M1' & period == 'S1')
dim(alldata1)
	
names(alldata1)[102:103]
	
a = alldata1 %>% select(contains('__'))
a = a[,which(specnumber(a, MARGIN=2)>=minocc)]
dim(a)
	
alldata1 = cbind(select(alldata1, -contains('__')), a)
dim(alldata1)
	
# ... 80% of data as training ...
# random selection
num.sample = dim(a)[1]
select.percent = .8
ssdd = 100 		# please keep this value to make results comparable!
set.seed(ssdd)
a = base::sample(1:num.sample, round(num.sample*select.percent))
	
# [1] 74 78 23 70  4 55 85  7 81 83 43 61 12 51 72 18 25  2 75 68 69 52 48 32 21
# [26] 27 39 57 16 11 67 71  6 29 45 30 53 79 86 31 33 49 82 28 47 41 87 42 24 80
# [51]  1  9 20 14 35 40  3 34 84 19 46 63 44 36 26  5 15 22 58 76
otu = select(alldata1, contains('__'))
otu.train = otu[a,]
dim(otu.train)
otu.test = otu[-a,]
	
envnames = c("insideHJA", "elevation_m", "canopyHeight_m", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")
	
ori.env = select(left_join(select(alldata1, envnames, "SiteName"), select(dem_stats, 'SiteName', 'tri.pt'), by=c('SiteName'='SiteName')), -'SiteName')
# add 'roughness' -> tri.pt
ori.env.train = ori.env[a,]
ori.env.test = ori.env[-a,]
str(ori.env.test)
	
ori.XY = select(alldata1, starts_with('UTM'))
ori.XY.train = ori.XY[a,]
ori.XY.test = ori.XY[-a,]
str(ori.XY.train)
	
# ... view data ...
par(mfrow=c(1,2))
hist(ori.env.train$elevation_m)
hist(ori.env$elevation_m)
	
ggplot(ori.XY, aes(UTM_E, UTM_N)) + geom_point() + geom_point(data=ori.XY.train, aes(colour='red')) + scale_colour_manual(labels = c('training'), values = c("red"))
	
```


```{r transform-data}
# 0/1
if (abund == 'pa')
{
	otu.train = as.data.frame((otu.train>0)*1)
	otu.test = as.data.frame((otu.test>0)*1)
}
	
# .. env data
scale.env.train.all = select(ori.env.train, -'insideHJA') %>% scale() 
str(scale.env.train.all)
	
scale.env.train = data.frame(scale.env.train.all) %>% add_column(insideHJA=as.factor(ori.env.train$insideHJA), .before=names(ori.env.train)[2])
str(scale.env.train)
	
dd.env.scaler = data.frame(t(data.frame(env.mean = attr(scale.env.train.all, "scaled:center"), env.sd = attr(scale.env.train.all, "scaled:scale"))))
str(dd.env.scaler)
	
rm(scale.env.train.all)
	
scale.env.test = as.data.frame(do.call(rbind, apply(select(ori.env.test, -'insideHJA'), 1, function(x){(x-dd.env.scaler['env.mean',])/dd.env.scaler['env.sd',]} ) )) %>% add_column(insideHJA=as.factor(ori.env.test$insideHJA), .before=names(ori.env.test)[2])
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
par(mfrow=c(2,2))
hist(ori.env.train[,2],xlim=c(1000,5500), breaks = 10)
hist(ori.env.test[,2],xlim=c(1000,5500), breaks = 10)
hist(scale(ori.env.test[,2]))
hist(scale.env.test[,2])
	
par(mfrow=c(1,2))
hist(XY.test[,2])
hist(XY.train[,2])
	 
rm(dd.env.scaler, dd.xy.scaler)
	
```


```{r model-DNN}
# one model to see how long it takes in ADA

# make sure input of sjsdm are numeric matrix
s.otu.train = as.matrix(otu.train)
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
# set variables
formula.env = 'envDNN'
lambda.env = .1
alpha.env = .5
lambda.sp = .1
alpha.sp = .9 
hidden.sp = c(50L,50L,10L)
acti.sp = 'relu'
drop = .3
	
model.train = sjSDM(Y = s.otu.train,
			  env = DNN(data=scale.env.train, formula = ~.,
					  hidden=hidden[[hiddenN]], lambda = lambda.env, alpha = alpha.env, activation=acti.sp, dropout=drop, bias = T),
 			  biotic = bioticStruct(lambda=lambda.env, alpha=alpha.env, on_diag=F, inverse = FALSE),
			  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp, alpha=alpha.sp),
			  learning_rate = 0.003, # 0.003 recommended for high species number 
			  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L 
			 )
	
names(model.train)
	
# saveRDS(model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS',  glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{date.model.run}_{formula.env}.RDS')) )
	
```

