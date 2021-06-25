# Jun 24, 2021
# data 's1', train&test random

```{r setup}
rm(list=ls())
q()
	
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm_general_model_outputs')
	
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
date.model.run = 'ele'  # '20210310' , 20210221, 20210429 20210514 20210523
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
	
# ... load from Hmsc_CD ................................
load(here('source','envVars.rdata'))		#; load(here("source",'rndVars.rdata'))
sort(names(allVars))
# 67
	
```


```{r split-train-test}
# ..... data belong to same sample are assigned to either train/test
# ... 75% of data as training ...
# random selection
num.point = nrow(alldata1)
num.sample = length(unique(alldata1$SiteName))
	
ran.tt = data.frame(table(alldata1$SiteName) )
colnames(ran.tt) = c('SiteName', 'NumTrap')
	
select.percent = .75
ssdd = 88
set.seed(ssdd)
	
sel = 0; num.train = round(num.point*select.percent); num.test = num.point - num.train
while ( sel != num.test ) {
	x1 = (num.train-table(ran.tt$NumTrap)[[1]])/2
	x2 = num.train - table(ran.tt$NumTrap)[[2]]*2
	
	table(ran.tt$NumTrap)
	num1 = base::sample(x1:x2, 1)
	Sseq = base::sample( ran.tt$SiteName, num1)
	# sample.test
	sel = sum(ran.tt$NumTrap[ran.tt$SiteName %in% Sseq])
}
	
Sseq
# [1] HOBO-233 SM-09    HOBO-345 280500   HOBO-055 301004   287338   544292  
# [9] 286789   HOBO-216 647830   124031   150801   HOBO-335 650846   135431  
#[17] HOBO-265 652934   HOBO-357 372502   361840 
	
# add train/test assignment
alldata1$assign = rep('xx', nrow(alldata1))
alldata1$assign[ alldata1$SiteName %in% Sseq ] = 'test'
alldata1$assign[ alldata1$assign=='xx' ] = 'train'
	
# remove rare species
alldata2 = alldata1 %>% filter(period == 'S1' & assign == 'train')
a = alldata2 %>% select(contains('__'))
a = a[,which(specnumber(a, MARGIN=2)>=minocc)]
dim(a); rm(alldata2)
	
alldata1 = select(alldata1, -contains('__'), all_of(names(a)))
rm(a); dim(alldata1)
	
if (preocc=='Dip') {
	alldata1 = select(alldata1, -contains('__'),  contains('__Insecta_Diptera'))
	dim(alldata1)
}
	
dim(select(alldata1, contains('__')))
sort(names(alldata1)[1:126])
names(alldata1)[1:126]
	
a = sort(alldata1$SiteName)
alldata1 = left_join(select(alldata1, SiteName, trap, period, starts_with('UTM'), elevation_f, oldGrowthIndex, YrsSinceDist, contains('__') ), allVars, by='SiteName')
dim(alldata1)
	
table(a == sort(alldata1$SiteName))(
sort(names(alldata1)[1:100]))
sort(names(alldata1)[223:297])
	
envvars = alldata1 %>% dplyr::select(!contains("__"), -starts_with("nor")) %>% 
	mutate(uniqueID = paste(SiteName, trap, period, sep = "_"),
         elevation_m = elevation_f * 0.3048, ## convert to metres
#         canopyHeight_m = canopyHeight_f * 0.3048,
         lg_DistStream = log(DistStream + 0.001),
         lg_DistRoad = log(DistRoad + 0.001),
         lg_YrsDisturb = log(YrsSinceDist + 0.001),
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
	
str(select(alldata1, contains('tpi')))
	
vars
	

```


```{r format-data}
if (formula.env =='vars8' & date.model.run=='ele') {vars = replace(vars, vars=='minT_annual', 'elevation_m')}
vars; formula.env
ori.env = select(envvars, all_of(vars))
str(ori.env)
	
dim(alldata1); dim(ori.env)
indNA <- complete.cases(ori.env)		# if having missing values
table(indNA)
alldata1 = alldata1[indNA,]
ori.env = ori.env[indNA,]
	
dim(envvars); dim(ori.env)
	
ttindex = which(alldata1$trap %in% 'M1'); vars = 'o'
	
rownames(alldata1)
	
env.test = envvars[-ttindex,]
env.test$SiteName
env.train = envvars[ttindex,]
env.train$SiteName
index = which(env.train$SiteName %in% env.test$SiteName)
	
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
save(XY.train,XY.test,scale.env.train,scale.env.test,otu.train,otu.test, file=here(datapath,'source', glue('forada_data_{period}_random_min{minocc}{preocc}_{date.model.run}_{formula.env}.rdata')))
	

```

```{r sjsdm}
if (abund == 'pa')
{
	otu.train = as.data.frame((otu.train>0)*1)
	otu.test = as.data.frame((otu.test>0)*1)
	dim(otu.test)
}
	
str(otu.train) 
names(scale.env.test); formula.env
	
s.otu.train = as.matrix(otu.train)	
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
s.otu.test = as.matrix(otu.test)		
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	

formula.env; abund; minocc
	
lambda.sp = .5
alpha.sp =  .5		# 1.0000000 0.6666667, (pa) 0.5
lambda.bio = .5
alpha.bio =  .5		# 1 0, (pa) 0.5
	
lambda.env = .5
alpha.env =  .5
drop = .2
ind=1
	
itern = 170L; iitt='-i170'
samn = 5000L; sstt=''
	
lrn = 0.002; lltt= '-lr002'
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
	
str(s.otu.train)
str(scale.env.test)
	

model.train = sjSDM(Y = s.otu.train,
	  env = DNN(data=scale.env.train, formula = ~., 
	  hidden=hidden[[ ind ]], lambda = lambda.env[ind], alpha = alpha.env[ind], activation='relu', dropout= drop[ind], bias=T),
	  
	  biotic = bioticStruct(lambda=lambda.bio[ind], alpha=alpha.bio[ind], on_diag=F, inverse = FALSE),
	  
	  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp[ind], alpha=alpha.sp[ind]),
	  
	  device = "cpu", learning_rate = lrn[ind], # 0.003 recommended for high species number 
	  step_size = NULL, iter = itern, family=stats::binomial('probit'), sampling = samn # 150L, 5000L
)  
# saveRDS(model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'random', glue('{formula.env}_{date.model.run}'), glue('s-jSDM{iitt}{sstt}{lltt[ind]}_{period}_{trap}_{abund}_min{minocc}{preocc}_{formula.env}_trial_{date.model.run}.RDS')) )
	
model.train = readRDS( here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'random', glue('{formula.env}_{date.model.run}'), glue('s-jSDM{iitt}{sstt}{lltt}_{period}_{trap}_{abund}_min{minocc}{preocc}_{formula.env}_trial_{date.model.run}.RDS')) )
	
# pdf(here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'random', glue('{formula.env}_{date.model.run}'), glue('model-history{iitt}{sstt}{lltt}_{period}_{trap}_{abund}_min{minocc}{preocc}_{formula.env}_{date.model.run}.pdf')), width=5, height=5)
	
par(cex=.8)
plot(model.train$history)
#mtext(paste0(names(maxdd)[i],': ', round(select(tuning.dd[pp,],names(maxdd)[i])[i,],3)), side=3,line=-1, adj=.95)
#mtext(paste0(names(tuning.dd)[c(2:7,9:10)],': ',select(tuning.dd[pp,],names(tuning.dd)[c(2:7,9:10)])[i,]), side=3,line=seq(-2,-9), adj=.95)
	
model.train$history
	
dev.off()

```


```{r explore-correlation}
names(envvars)
env = subset(envvars, trap=='M1')
	
env = select(env, -starts_with('mean.'),-starts_with('lg_'), -ends_with('_f'), -uniqueID, -SiteName, -trap, -period, -starts_with('UTM_'), -contains('_20180717'), -contains('_20180802'), -contains('_20180818'))
# -LC08_045029_20180726_B10_100m, -tpi, -lysis_ratio,-COISpike_sum, -be10, -G_20180717,, -ht
# see 'test linearity', some papams are exactly the same but wit diff names
sort(names(env))
	
str(select(alldata, oldGrowthIndex))
env$oldGrowthIndex
	
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
a = select(scale.env, be500, elevation_m)
a = a[order(a[,1]),]
plot(1:nrow(a),a[,1]) #, ylim=c(0,max(a[,2])))
points(1:nrow(a), a[,2], col='red')
points(1:nrow(a), a[,3], col='blue')
	
# ....................... test vars containing NA ....................... 
for (i in 1:ncol(scale.env)) {
	if (any(is.na(scale.env[,i]))) { 
		print(names(scale.env)[i])
		print(which(is.na(scale.env[,i]))) }
}
scale.env$cut_r
	
#scale.env = select(scale.env, -tpi1k, -tpi500)[c(-1,-6),] # ,-'ndvi_p50_100m'
scale.env = select(scale.env, -cut_r)
names(scale.env)
	
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
	
dd = select(scale.env, -contains('T_annual'))
names(dd)
	
while  (maxvif>=8 ) {
	if (varrem!='a') { dd = select(dd, -all_of(varrem)) }
	
	vif = corvif(dd)
	order = order(vif$GVIF, decreasing=F)
	vif.count = data.frame(var=names(dd)[order], vif=vif[order,])
	maxvif = max(vif.count$vif)
	varrem = as.character(vif.count$var[ncol(dd)])
	
	print(c(ncol(dd), varrem))
}
	
sort(as.character(vif.count$var)) == sort(vars)
sort(names(scale.env))
	
#allvars1 = names(scale.env); allvars2 = names(select(scale.env, -contains('T_annual'))); select2 = vif.count$var; vars8 = vars
#save(allvars1, allvars2, select2, vars8, file=here('source','vars8-test.rdata'))
#rm(allvars1, allvars2, select2, vars8)
	
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

```{r heat-map}
abund; formula.env; minocc
datapath=glue('../12_sjsdm_general_model_outputs')
	
# ...... load from rdata ........
load(here(datapath,'source', glue('forada_data_{period}_m1m2_min{minocc}_{date.model.run}_{formula.env}.rdata')))
	
# ...... heat map .....
#str(hmp_otus)
#hmp_otus$lineage
	
# .... format data for heat map ....
str(otu.train)
otu.train = as.data.frame((otu.train>0)*1)
otu.test = as.data.frame((otu.test>0)*1)
names(otu.train) == names(otu.test)
	
m1 = bind_rows((otu.train %>% rownames_to_column(var='sample.id') %>% add_column(., id='m1')), 
		(otu.test %>% rownames_to_column(var='sample.id') %>% add_column(., id='m2'))) %>% mutate(idfull = paste0(sample.id,'_', id))
m1 = select(m1,contains('__'), idfull) %>% gather(name, presence, -idfull) %>% spread(idfull, presence) %>% 
		separate(name, c('dummy','taxa'), sep='__', remove=F) %>% separate(taxa, c('taxa','dummy'), sep='_BOLD') %>% separate(taxa, c('class','order','family','genus','species'), sep='_') %>% select(., -dummy)
	
m1$presence
m1[1:5,1:4]
	
head(m1) 
dim(m1)
unique(m1$order)
	
m1$species = na_if(m1$species, 'NA'); m1$genus = na_if(m1$genus, 'NA'); m1$family = na_if(m1$family, 'NA')
# table(is.na(m1$family)); table(is.na(m1$order))
	
# ..... steps for heat map .....
m1h = parse_tax_data(m1, class_cols = c('class','order','family','genus','species'), named_by_rank = TRUE); m1h
m1h %>% taxon_ranks
# 
m1h$data$otu_table <- calc_obs_props(m1h, data = "tax_data", cols = names(select(m1, contains('_m'))) )
m1h$data$tax_table <- calc_taxon_abund(m1h, data = "otu_table", cols = names(select(m1, contains('_m'))) )
	
m1h$data$diff_table = compare_groups(m1h, data='tax_table', cols=names(select(m1, contains('_m'))), groups=as.character(sapply(names(select(m1, contains('_m'))), function(a) str_split(a, '_')[[1]][2])) )
unique(m1h$data$diff_table$taxon_id)
m1h$data$tax_data
	
m1h$data$type_abund <- calc_group_mean(m1h, "tax_table", cols=names(select(m1, contains('_m'))), groups=as.character(sapply(names(select(m1, contains('_m'))), function(a) str_split(a, '_')[[1]][2])))
	
nodraw = c("species", "genus"); tdraw = 'family'
	
m1h %>% filter_taxa(! taxon_ranks %in% nodraw) %>%
heat_tree( node_label = taxon_names,
            node_size = n_obs, 
            node_color = log2_median_ratio, node_color_range = c("cyan", "gray", "tan"),
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio (median)", 
            margin_size = c(0.03, 0.03, 0.03, 0.03)
#           , output_file = here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'descriptive', glue('heat-map_{tdraw}_M1M2_{abund}_min{minocc}{preocc}.pdf'))
)
# layout = "davidson-harel", initial_layout = "re",
	
tdraw = 'noNA'
	
m1h %>% filter_taxa( taxon_names != 'NA') %>%
heat_tree( node_label = taxon_names,
            node_size = n_obs, node_label_size_range = c(0.006,0.02),
            node_color = log2_median_ratio, node_color_range = c("cyan", "gray", "tan"),
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio (median)", 
            margin_size = c(0.03, 0.03, 0.03, 0.03)
	
#           , output_file = here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'descriptive', glue('heat-map_{tdraw}_M1M2_{abund}_min{minocc}{preocc}.pdf'))
)
	
# ... only m2
m1h %>% 
  filter_taxa(nchar(taxon_names)!=0) %>%
  heat_tree(node_label = taxon_names,
            node_size = m2, 
            node_color = m2, 
            layout = "da", initial_layout = "re", 
            title = "Taxa in m2")
				
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
