# Feb 22, 2021
# train = 's1m1', test = 's1m2'

#'oldvars') {vars = c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")}
	
#if (formula.env=='newvars') {vars = c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")} 
	
#if (formula.env=='vars2') {vars = c("be10", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")} 
	
#if (formula.env=='vars3') { vars <- c("be10", "slope", "tri", "Nss", "Ess","ndmi_stdDev", "ndvi_p5", "ndvi_p50", "ndvi_p95", "ndmi_p5", "ndmi_p50", "ndmi_p95", "savi_p50", "LC08_045029_20180726_B1", "LC08_045029_20180726_B3", "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7", "LC08_045029_20180726_B10", "ndmi_stdDev_100m", "ndvi_p5_100m", "ndvi_p50_100m", "ndvi_p95_100m", "ndmi_p5_100m", "ndmi_p50_100m", "ndmi_p95_100m", "savi_p50_100m", "LC08_045029_20180726_B1_100m", "LC08_045029_20180726_B3_100m", "LC08_045029_20180726_B4_100m", "LC08_045029_20180726_B5_100m", "LC08_045029_20180726_B7_100m", "LC08_045029_20180726_B10_100m", "tpi250",  "tpi500", "tpi1k" , "ht", "ht.r250", "ht.r1k", "cov2_4.r250", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r1k", "mTopo", "cut.r1k.pt", "insideHJA", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "l_p25", "l_rumple") }

```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm_code_reviewed')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan','pROC', 'gridExtra','ggeffects','corrplot')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.6
	
source(here('source', 'corvif-source.r'))
	

```


```{r set-names}
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
trap <- "M1"; period = "S1"
date.model.run = 20210523  # '20210310' , 20210221, 20210429 20210514
abund = 'pa'		# 'qp','pa'  !!! change accordingly
formula.env = 'vars8'		# 'vars2', vars3, vars6 , vars5 , vars8 vars6x
# vars8 -> from VIF
preocc = ''		# 'Dip'
minocc = 6		# 6 , 20
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.6' # package version  0.1.3.9000
	
# names for graph
sjsdmVfolder = glue( 'sjsdm-{sjsdmV}' )
	
```


```{r load-data}
# ..... load data ......
alldata = read.csv(here(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
dim(alldata)
names(alldata)[110:125]
# 124 Var
	
alldata1 = alldata %>% filter(period == period[1]) #& trap == trap[1]
dim(alldata1); unique(alldata1$period)
	
# ..... delete wrong data points
alldata1[c(115,86),]$SiteName
# [1] "SM-06"    "HOBO-270"
alldata1[c(115,86),]$period; alldata1[c(115,86),]$trap
	
# data are not wrong !!!
dim(alldata1)
	
# remove rare species
alldata2 = alldata1 %>% filter(period == period[1] & trap == trap[1])
a = alldata2 %>% select(contains('__'))
a = a[,which(specnumber(a, MARGIN=2)>=minocc)]
dim(a); rm(alldata2)
	
alldata1 = select(alldata1, -contains('__'), all_of(names(a)))
rm(a); dim(alldata1)
	
if (preocc=='Dip') {
	alldata1 = select(alldata1, -contains('__'),  contains('__Insecta_Diptera'))
	dim(alldata1)
}
	
# ... load from Hmsc_CD ................................
load(here('source','envVars.rdata'))		#; load(here("source",'rndVars.rdata'))
str(allVars)
	
a = sort(alldata1$SiteName)
alldata1 = left_join(select(alldata1, contains('__'), SiteName, trap, period, starts_with('UTM') ), allVars, by='SiteName')
alldata1$meanT_annual[alldata1$SiteName == 'SM-01']
dim(alldata1)
table(a == sort(alldata1$SiteName))
	
envvars = alldata1 %>% dplyr::select(!contains("__"), -starts_with("nor")) %>% 
	mutate(uniqueID = paste(SiteName, trap, period, sep = "_"),
#         elevation_m = elevation_f * 0.3048, ## convert to metres
#         canopyHeight_m = canopyHeight_f * 0.3048,
         lg_DistStream = log(DistStream + 0.001),
         lg_DistRoad = log(DistRoad + 0.001),
#         lg_YrsDisturb = log(YrsSinceDist + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001)) %>% 
         mutate(insideHJA = factor(insideHJA)) 
dim(envvars)
names(envvars)
	
# from VIF (chunk pure VIF)
if (formula.env =='vars8') { load(here(outputpath,'sjsdm_general_outputs', 'sjsdm-0.1.3.9000', 'descriptive','selected-vars.rdata'))
	vars = as.character(vif.count$var)
}
	
if (formula.env =='vars6x') { 
cv = '5CV'
load( here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('vars6_20210429'),'flashlight','sum-flash', glue('data-var-imp-sum_vars6_{abund}_{cv}_{period}_{trap}_min6_20210429')) ) 
vars = ddlist[[3]]
}
	
# ... load from Hmsc_CD ...........................................
# ..... make tables (test->m2, train->m1) .....
if (formula.env =='vars5') { vars <- c("ht30", "gt4_r30", "gt4_500", "cut_r500", "cut_r250", "cut40_r500", "cut40_r250", "be30", "tri30", "slope30", "Nss30", "Ess30", "twi30", "tpi250", "tpi500", "tpi1k", "l_rumple", "insideHJA", "ndmi_stdDev_r100", "nbr_stdDev_r100", "ndvi_p5_r100", "ndvi_p50_r100", "ndvi_p95_r100", "ndmi_p5_r100", "ndmi_p50_r100", "ndmi_p95_r100", "LC08_045029_20180726_B1", "LC08_045029_20180726_B3", "LC08_045029_20180726_B4", "LC08_045029_20180726_B5", "LC08_045029_20180726_B7", "LC08_045029_20180726_B10", "minT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m") }
# cannot! no variables of '`lg_DistStream`, `lg_DistRoad`, `lg_cover2m_max`, `lg_cover2m_4m`, and `lg_cover4m_16m`'
if (formula.env =='vars6') { vars <- c("ht30","gt4_r30","gt4_250", "cut_r250","cut40_r1k","cut40_r250", "be30", "slope30", "Nss30","Ess30","twi30","tpi250","tpi1k", "l_Cover_2m_max", "l_Cover_4m_16m", "l_rumple", "DistStream","DistRoad", "ndmi_stdDev_r100", "ndmi_stdDev_r500","nbr_stdDev_r100", "nbr_stdDev_r500", "ndvi_p5_r100","ndvi_p5_r500", "ndvi_p50_r100","ndvi_p50_r500", "ndmi_p50_r100", "ndmi_p50_r500", "LC08_045029_20180726_B4","LC08_045029_20180726_B5","LC08_045029_20180726_B10", "minT_annual","precipitation_mm") }
if (formula.env =='vars7') { vars <- c("gt4_r30", "gt4_500", "cut_40msk", "cut_r1k", "cut_r250", "cut40_r1k", "cut40_r250", "tri30", "Nss30", "Ess30", "twi30", "tpi250", "tpi1k", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_rumple", "DistStream", "DistRoad", "insideHJA", "ndmi_stdDev_r100", "nbr_stdDev_r250", "ndvi_p5_r100", "ndvi_p95_r500", "ndmi_p95_r100", "ndmi_p95_r500", "LC08_045029_20180726_B3", "LC08_045029_20180726_B5", "LC08_045029_20180726_B10", "minT_annual") }
	
str(select(allVars, tpi250))
str(select(envvars, contains('tpi')))
str(select(alldata1, contains('tpi')))
	
vars
	

```


```{r format-data}
vars; formula.env
ori.env = select(envvars, all_of(vars))
str(ori.env)
	
dim(alldata1); dim(ori.env)
indNA <- complete.cases(ori.env)		# if having missing values
table(indNA)
alldata1 = alldata1[indNA,]
ori.env = ori.env[indNA,]
	
ttindex = which(alldata1$trap %in% 'M1'); vars = 'o'
	
ori.env.train = ori.env[ttindex,]
ori.env.test = ori.env[-ttindex,]
dim(ori.env.train)
	
ori.XY = select(alldata1, starts_with('UTM'))
ori.XY.train = ori.XY[ttindex,]
ori.XY.test = ori.XY[-ttindex,]
str(ori.XY.train)
	
otu = select(alldata1, contains('__'))
otu.train = otu[ttindex,]
dim(otu.train)
otu.test = otu[-ttindex,]
dim(otu.train); dim(otu.test)
	
# .. env data
fenv = names(ori.env.train)[sapply(ori.env.train, is.factor )]; fenv
locfenv = which(sapply(1:ncol(ori.env.train), function(x) {is.factor(ori.env.train[,x])} )==T); locfenv
scale.env.train.all = select(ori.env.train, -all_of(fenv)) %>% scale() 
str(scale.env.train.all)
	
scale.env.train = data.frame(scale.env.train.all) %>% add_column(select(ori.env.train, all_of(fenv))) # , .before=names(ori.env.train)[locfenv+1]
# insideHJA factor !!!
str(scale.env.train)
	
dd.env.scaler = data.frame(t(data.frame(env.mean = attr(scale.env.train.all, "scaled:center"), env.sd = attr(scale.env.train.all, "scaled:scale"))))
str(dd.env.scaler)
	
rm(scale.env.train.all)
	
scale.env.test = as.data.frame(do.call(rbind, apply(select(ori.env.test, -all_of(fenv)), 1, function(x){(x-dd.env.scaler['env.mean',])/dd.env.scaler['env.sd',]} ) )) %>% add_column(select(ori.env.test, all_of(fenv)))
#add_column(insideHJA=as.factor(ori.env.test$insideHJA), .before=names(ori.env.test)[2])
str(scale.env.test)
dim(scale.env.train); dim(scale.env.test)
	
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
dim(XY.train); dim(XY.test)
	

```

```{r save-data}
datapath=glue('../12_sjsdm_general_model_outputs')
save(XY.train,XY.test,scale.env.train,scale.env.test,otu.train,otu.test, file=here(datapath,'source', glue('forada_data_{period}_m1m2_min{minocc}{preocc}_{date.model.run}_{formula.env}.rdata')))
	

```

```{r CD-correlation}
# ... run with only vars from CD
env = allVars
env = subset(env, trap=='M1') %>% select(-SiteName)
	
env1 = subset(envvars, trap=='M1')
allVars$SiteName; env1$SiteName
allVars1 = left_join(select(env1, SiteName), allVars)
names(allVars1)
	
env = left_join(allVars1, select(env1, SiteName, clearcut, YrsSinceDist, elevation_m, canopyHeight_m, EVI_20180726, B_20180726, G_20180726, W_20180726), by = 'SiteName') %>% select(-SiteName)
names(env)
	
fenv = names(env)[sapply(env, is.factor )]; fenv
locfenv = which(sapply(1:ncol(env), function(x) {is.factor(env[,x])} )==T); locfenv
scale.env = select(env, -all_of(fenv)) %>% scale() 
str(scale.env)
	
scale.env = data.frame(scale.env) %>% add_column(select(env, all_of(fenv))) # , .before=names(ori.env.train)[locfenv+1]
str(scale.env)
sort(names(scale.env))
	
for (i in 1:ncol(scale.env)) {
	if (any(is.na(scale.env[,i]))) { 
		print(names(scale.env)[i])
		print(which(is.na(scale.env[,i]))) }
}
# ??? 'cut_r'
	
scale.env = select(scale.env, -cut_r, -aspect30)
nchoose = 40
	
# ............... pure VIF ...............
# loop for better result
dd = scale.env; varrem = 'a'; maxvif = 100
while  (maxvif>=8 ) {
	if (varrem!='a') { dd = select(dd, -all_of(varrem)) }
	
	vif = corvif(dd)
	order = order(vif$GVIF, decreasing=F)
	vif.count = data.frame(var=names(dd)[order], vif=vif[order,])
	maxvif = max(vif.count$vif)
	varrem = as.character(vif.count$var[ncol(dd)])
	
	print(c(ncol(dd), varrem))
	}
	
pairs(select(scale.env,vif.count$var),lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)
	
# save file
# save(vif.count, file=here(outputpath,'sjsdm_general_outputs', sjsdmVfolder, 'descriptive','selected-vars.rdata'))
load(here(outputpath,'sjsdm_general_outputs', sjsdmVfolder, 'descriptive','selected-vars.rdata'))
	
```


```{r explore-correlation}
names(envvars)
env = subset(envvars, trap=='M1')
	
env = select(env, -starts_with('mean.'),-starts_with('lg_'), -ends_with('_20180717'), -ends_with('_20180802'), -ends_with('_20180818'), -ends_with('_f'), -uniqueID, -SiteName, -trap, -period, -lysis_ratio,-COISpike_sum, -starts_with('UTM_'),-l_Cover_2m_4m, -l_Cover_4m_16m, -LC08_045029_20180726_B10_100m, -be10, -G_20180717,-ht,-tpi)
# see 'test linearity', some papams are exactly the same but wit diff names
names(env)
	
fenv = names(env)[sapply(env, is.factor )]; fenv
locfenv = which(sapply(1:ncol(env), function(x) {is.factor(env[,x])} )==T); locfenv
scale.env = select(env, -all_of(fenv)) %>% scale() 
str(scale.env)
	
scale.env = data.frame(scale.env) %>% add_column(select(env, all_of(fenv))) # , .before=names(ori.env.train)[locfenv+1]
str(scale.env)
sort(names(scale.env))
	
# ....................... test linearity ....................... 
a = select(env, slope,tri)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1], ylim=c(min(a[,2]),max(a[,1])))
points(1:nrow(a), a[,2], col='red')
	
a = select(env, ht.r500,ht.r250)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1], ylim=c(min(a[,2]),max(a[,1])))
points(1:nrow(a), a[,2], col='red')
	
a = select(env, B10_20180717, LC08_045029_20180726_B10)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1])
points(1:nrow(a), a[,2], col='red')
	
a = select(env, ht, l_p95, l_p95_all)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1], ylim=c(min(a[,3]),max(a[,1])))
points(1:nrow(a), a[,2], col='red')
points(1:nrow(a), a[,3], col='blue')
	
# !!!
a = select(env, G_20180717, B5_20180717, G_20180717)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1], ylim=c(min(a[,3]),max(a[,3])))
points(1:nrow(a), a[,2], col='red')
points(1:nrow(a), a[,3], col='blue')
	
# !!!!!!!!!!
a = select(env, LC08_045029_20180726_B10_100m, LC08_045029_20180726_B10)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1], ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
	
a = select(scale.env, l_Cover_4m_16m, cov4_16)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1]) #, ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
	
# !!!
a = select(scale.env, l_Cover_2m_4m_all, cov2_4)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1]) #, ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
	
a = select(scale.env, LC08_045029_20180726_B1, B1_20180717)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1]) #, ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
	
# !!!!!!!!!
a = select(scale.env, be10, be500, elevation_m)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1]) #, ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
points(1:nrow(a), a[,3], col='blue')
	
# ....................... test linearity ....................... 
for (i in 1:ncol(scale.env)) {
	if (any(is.na(scale.env[,i]))) { 
		print(names(scale.env)[i])
		print(which(is.na(scale.env[,i]))) }
}
	
scale.env0 = scale.env
scale.env = select(scale.env, -tpi1k, -tpi500)[c(-1,-6),] # ,-'ndvi_p50_100m'
	
nchoose = 40 		# num of variables to choose
# ........ cor-vif (no factorical allowed in 'cor')
env.num = select(scale.env,-fenv)
str(env.num)
a = cor(env.num)
	
order.AOE <- corrMatOrder(a, order = "AOE")
a.AOE <- a[order.AOE,order.AOE]
	
corrplot(a)
corrplot(a.AOE)
	
# ... remove entries with large correlation
nchoose
total = ncol(env.num)*(ncol(env.num)-1)/2; chose = nchoose*(nchoose-1)/2
	
a.cor = a
#a.cor[upper.tri(a, diag=T)] = NA; a.cor
diag(a.cor) = NA
hist(a.cor)
	
sort(abs(a.cor))[1:(chose*2)]
last = tail(sort(abs(a.cor))[1:(chose*2)], n=1)
a.a = a.cor
a.a[abs(a.a)>last] = NA
table(is.na(a.a))
hist(a.a)
corrplot(a.a)
	
# count non NA for each variable
cov.count = data.frame(var=attr(a.a,'dimnames')[[1]], count=rep(0,length.out=nrow(a.a)), sum = rep(0,length.out=nrow(a.a)))
for (i in 1:dim(a.a)[1]) { 
	cov.count$count[i] = length(which(!is.na(a.a[i,])))
	cov.count$sum[i] = sum(abs(a.cor[i,]), na.rm=T) }
str(cov.count)
	
par(mfrow=c(1,2))
hist(cov.count$count); hist(cov.count$sum)
	
# ... sum
cov.sum = cov.count[order(cov.count$sum, decreasing=F),]
sort(cov.sum$var[1:nchoose])
	
# .... II. vif with factorial variables
scale.env2 = select(scale.env, cov.sum$var[1:nchoose], fenv)		# cov.sum , cov.count
str(scale.env2)
	
corvif(scale.env2)
	
vif = corvif(scale.env2)
order = order(vif$GVIF, decreasing=F)
covvif.count = data.frame(var=names(scale.env2)[order], vif=vif[order,])
covvif.count
covvif.count[1:33,1]
#VIF < 30
	
var.select = as.character(covvif.count[1:33,1])
	
# save(var.select, file=here(outputpath,'sjsdm_general_outputs', sjsdmVfolder, 'descriptive','selected-vars.rdata'))
load(here(outputpath,'sjsdm_general_outputs', sjsdmVfolder, 'descriptive','selected-vars.rdata'))
	
scale.env3 = select(scale.env, covvif.count[1:33,1])		# cov.sum , cov.count
str(scale.env3)
	
corvif(scale.env3)
	
vif = corvif(scale.env3)
order = order(vif$GVIF, decreasing=F)
covvif.count = data.frame(var=names(scale.env3)[order], vif=vif[order,])
covvif.count
	
pairs(scale.env3,lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)
	
# ............... pure VIF ...............
# loop for better result
dd = scale.env; varrem = 'a'; maxvif = 100
while  (maxvif>=8 ) {
	if (varrem!='a') { dd = select(dd, -all_of(varrem)) }
	
	vif = corvif(dd)
	order = order(vif$GVIF, decreasing=F)
	vif.count = data.frame(var=names(dd)[order], vif=vif[order,])
	maxvif = max(vif.count$vif)
	varrem = as.character(vif.count$var[ncol(dd)])
	
	print(c(ncol(dd), varrem))
	}
	
pairs(select(scale.env,vif.count$var),lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)
	
vif = corvif(scale.env)
order = order(vif$GVIF, decreasing=F)
vif.count = data.frame(var=names(scale.env)[order], vif=vif[order,])
# only one-time VIF
vif = corvif(scale.env)
order = order(vif$GVIF, decreasing=F)
vif.count = data.frame(var=names(scale.env)[order], vif=vif[order,])
vif.count
	
pairs(select(scale.env,vif.count$var[1:33]),lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)
	
env3 = select(scale.env,vif.count$var[1:33])
vif = corvif(env3)
order = order(vif$GVIF, decreasing=F)
vif.count = data.frame(var=names(env3)[order], vif=vif[order,])
vif.count
	
pairs(env3,lower.panel=panel.smooth2,upper.panel=panel.cor,diag.panel=panel.hist)
	
# ...... plot to see
covvif.count$type='covvif'; vif.count$type='vif'
	
nchoose = 33
compare = rbind(merge(vif.count[1:nchoose,c(1,3)], covvif.count[1:nchoose,c(1,3)], by='var') %>% select(var) %>% add_column(type='both'), anti_join(vif.count[1:nchoose,c(1,3)], covvif.count[1:nchoose,c(1,3)],by='var'), anti_join(covvif.count[1:nchoose,c(1,3)], vif.count[1:nchoose,c(1,3)], by='var'))
	
compare$var = str_replace(compare$var, 'LC08_045029_20180726','LC')
compare$var = str_replace(compare$var, 'l_Cover','l')
	
ss = 6
g1 = ggplot(subset(compare, type=='both'), aes(var)) + geom_bar() + theme(axis.text=element_text(size=ss)) + scale_x_discrete(limits=unique(subset(compare, type=='both')$var))
g2 = ggplot(subset(compare, type!='both'), aes(var,color=type)) + geom_bar() + theme(axis.text=element_text(size=ss), legend.position='bottom') + scale_x_discrete(limits=unique(subset(compare, type!='both')$var))
	
# pdf(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'descriptive', 'choose-var2.pdf'), width=12, height=5)
grid.arrange(g1,g2,nrow=2, heights=c(.44,.56))
	
dev.off()
	
```


```{r explore-correlation-snippet}

# ... count
cov.count = cov.count[order(cov.count$count, decreasing=T),]
sort(cov.count$var[1:nchoose])		#; sort(vif.count$var[1:nchoose])
	
```





```{r explore-graph}
# ..... explore m1,m2 ......
str(alldata1)
otu = as.data.frame((select(alldata1, contains('__'))>0)*1)
dim(otu)
	
exp = select(alldata1, 'UTM_E', 'UTM_N', 'trap', 'SiteName') %>% bind_cols(spprich=rowSums(otu>0))
str(exp)
	
which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'])
which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'] & exp$trap=='M1')
	
# .......... add HJA shape file ...........
# bring in HJA boundary
# https://data-osugisci.opendata.arcgis.com/datasets/74312b6130cb4e9b8c454ae1195f6482_9/data
pacman::p_load(sf)
	 
utm10N = 32610
hja = st_read(here(outputpath, 'HJA_Base_GIS_Layers_2018', "HJA_Boundary.shp"))
	
hja_bound = subset(hja, FP_NAME == "H.J. Andrew Experimental Forest") # just get the boundary
	
hja.utm = st_transform(hja_bound, crs = utm10N) # transform to UTM
class(hja.utm)
plot(st_geometry(hja.utm))
	
hjautm = as.data.frame(st_coordinates(hja.utm))
str(hjautm)
	
stat_sf_coordinates(data=hja.utm, mapping=aes())
	
# pdf(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'descriptive', glue('trap_s1_m1m2_min{minocc}.pdf')), height=6,width=6.5)
	
# pdf(here(outputpath,'prediction_outputs','sjsdm-graph', sjsdmVfolder,glue('coordinate_s1_m1m2_min{minocc}.pdf')), height=10,width=10)
	
# ggplot() + geom_sf(data = st_geometry(hja.utm), colour = "black", fill = NA)
	
hjashape = geom_polygon(data = hjautm, aes(x=X, y=Y), colour = "black", fill = NA) 
p1 = ggplot(exp, aes(UTM_E, UTM_N),shape=trap) + geom_point(aes(colour=trap,shape=trap)) + scale_shape_manual(values=c(19,3))
p2 = ggplot(exp, aes(UTM_E, UTM_N),colour=spprich) + geom_point(aes(colour=spprich)) + geom_point(colour='white') + geom_point(data=filter(exp,trap=='M2'), aes(colour=spprich)) + scale_colour_gradient2(low='black', mid='red',high='blue',midpoint=max(exp$spprich)/2) + ggtitle('M2,test') + theme(legend.position='none')
p3 = ggplot(exp, aes(UTM_E, UTM_N),colour=spprich) + geom_point(aes(colour=spprich)) + geom_point(colour='white') + geom_point(data=exp[which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'] & exp$trap=='M1'),], aes(colour=spprich)) + scale_colour_gradient2(low='black', mid='red',high='blue',midpoint=max(exp$spprich)/2) + ggtitle('M1, explain')
	
grid.arrange(p1,p2,p3, nrow=2, layout_matrix= rbind(c(1), c(2,3)), heights=c(.4,.6),widths=c(.44,.56))
	
p1 + hjashape
	
dev.off()
	

```
