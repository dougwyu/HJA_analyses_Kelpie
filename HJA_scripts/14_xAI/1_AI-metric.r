# Mar 04, 2021
# last modified: Feb 21, 2021
# play with tuned model (manually tuning in ADA cluster, 20210215-0221)

## explore vip package 
#pred_fun = function(model1) {    
#  apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=NULL, SP=NULL)) , along = -1L), 2:3, mean)
#}
	
#vip::vip(object = model.train, data=data.frame(XY.train,scale.env.train,otu.train[,1]), y = "otu.train...1.", method = "permute", metric = "auc", pred_wrapper = pred_fun, reference_class = 1)
	
#vi(model.train, method='firm')
## only works for certain types of models !!!
## doesn't work


```{r setup}
# set_here()
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/14_xAI')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vip','MetricsWeighted','flashlight')
	
conflict_prefer('vi','vip')
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')                                             
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.3.9000’
	
source(here("source", "xAI-function.r"))
source(here("source", 'sjsdm_function.r'))
	

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
date.model.run = 20210322   # '20210221' '20210310' 20210322
abund = 'qp'		# 'qp','pa'  !!! change accordingly
formula.env = 'vars3'		# 'newVars', 'vars2', 'oldVars' , 'vars3'
preocc = ''		# '' 'm1m2'
minocc = 6
cv = '5CV'		# '10CV'
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
```


```{r load-direct-data}
# ... read data ...
abund; formula.env; minocc; cv
nstep = 2000		# 2000 , 1000
datapath=glue('../12_sjsdm_general_model_outputs')
#nnn = mean_EVAL
	
# ...... from rdata ........
load(here(datapath,'source', glue('forada_data_{period}_m1m2_min{minocc}{preocc}_{date.model.run}_{formula.env}.rdata')))
tuning=read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune',glue('{formula.env}_{date.model.run}'),'best', glue('best_manual_tuning_sjsdm_{cv}_{trap}{period}_mean_AUC_{abund}_min_{minocc}{preocc}_nSteps_{nstep}.csv')), header=T, sep=',') 
	
str(tuning)
	
# 0/1
if (abund == 'pa')
{
	otu.train = as.data.frame((otu.train>0)*1)
	otu.test = as.data.frame((otu.test>0)*1)
	str(otu.test)
}
	
str(otu.train); abund
names(scale.env.test); formula.env
	
s.otu.train = as.matrix(otu.train)	
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
s.otu.test = as.matrix(otu.test)		
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	

```

```{r load-model}
lambda.env = tuning[,'lambda.env']; alpha.env = tuning[,'alpha.env']; lambda.sp = tuning[,'lambda.sp']
alpha.sp =  tuning[,'alpha.sp']; lambda.bio = tuning[,'lambda.bio']; alpha.bio =  tuning[,'alpha.bio']
hidden1 = tuning[,'hidden.ind']; acti = as.character(tuning[,'acti.sp']); drop = tuning[,'drop']    
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
i=1; metric='AUC.test_mean'
if (cv=='5CV' & formula.env=='newVars' & minocc==6) {i=4; metric='plr.test_mean'}
if (cv=='10CV' & abund=='qp' & formula.env=='newVars' & minocc==6) {i=4; metric='plr.test_mean'}
if (cv=='5CV' & abund=='qp' & formula.env=='vars2' & minocc==6 & preocc=='m1m2') {i=3; metric='plr.test_mean'}
if (cv=='5CV' & abund=='pa' & formula.env=='vars3' & minocc==6 & preocc=='') {i=3; metric='plr.test_mean'}
metric; i; select(tuning, metric)
	
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS',glue('{formula.env}_{date.model.run}'), glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_{cv}_min{minocc}{preocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{lambda.bio[i]}_{drop[i]}_{date.model.run}.RDS')) )
	
plot(model.train$history)
	

```

```{r flashlight-coord} 
# flashlight package
old=Sys.time(); newt = old
# for (i in 1:100 ) {
# for (i in 101:200 ) {
# for (i in 201:ncol(otu.train) ) {
for (i in 1:ncol(otu.train) ) {
	old=Sys.time()
	
	print(i)
	custom_predict <- function(model1,new_data) {
	newdd = select(new_data, -'otu.train...i.', -'UTM_E',-'UTM_N')
	spdd = select(new_data, 'UTM_E','UTM_N')
#		print(c(dim(newdd),dim(spdd)))
	apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=spdd)) , along = -1L), 2:3, mean)[,i]
}
	
	fl = flashlight::flashlight(model = model.train, data=data.frame(scale.env.train,XY.train, otu.train[,i]), y = "otu.train...i.", label = "oregon", predict_function = custom_predict)
	
	imp <- flashlight::light_importance(fl, m_repetitions = 3,type = "permutation")
	save(imp, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','coord', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
	vv = imp$data$variable[sapply(1:10, function(i) which(abs(imp$data$value)==sort(abs(imp$data$value),decreasing=T)[i]) ) ]
#	int <- flashlight::light_interaction(fl, pairwise=F, type='H', v=vv, grid_size = round(nrow(otu.train)*.8), n_max = nrow(otu.train), seed = 42)
	int <- flashlight::light_interaction(fl, pairwise=F, type='H', v=vv, grid_size = round(nrow(otu.train)*.5), n_max = round(nrow(otu.train)*.9), seed = 42)
	newt = Sys.time() - old
	print(c(i,newt),digit=3)
	
	save(int, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2','coord',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	}
	

```

```{r flashlight} 
# flashlight package
old=Sys.time(); newt = old
# for (i in 1:100 ) {
# for (i in 101:200 ) {
# for (i in 201:ncol(otu.train) ) {
for (i in 1:ncol(otu.train) ) {
	old=Sys.time()
	
	print(i)
	custom_predict <- function(model1,new_data) {
		newdd = select(new_data, -'otu.train...i.', -'UTM_E',-'UTM_N')
		spdd = select(new_data, 'UTM_E','UTM_N')
#		print(c(dim(newdd),dim(spdd)))
		apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=spdd)) , along = -1L), 2:3, mean)[,i]
}
	
	fl = flashlight::flashlight(model = model.train, data=data.frame(scale.env.train,XY.train, otu.train[,i]), y = "otu.train...i.", label = "", predict_function = custom_predict)
	
	imp <- flashlight::light_importance(fl, m_repetitions = 3,type = "permutation", v=names(scale.env.train))
	save(imp, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
	vv = imp$data$variable[sapply(1:10, function(i) which(abs(imp$data$value)==sort(abs(imp$data$value),decreasing=T)[i]) ) ]
	int <- flashlight::light_interaction(fl, pairwise=F, type='H', v=vv, grid_size = round(nrow(otu.train)*.5), n_max = round(nrow(otu.train)*.9), seed = 42)
	save(int, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
	newt = Sys.time() - old
	print(c(i,newt),digit=3)
	rm(imp, int)
#	load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{formula.env}_spp{i}.rdata' ) ))
	
#	plot(imp, fill = "darkred")
#	most_important(imp, 2)
	}
	
# ................... interaction .....................
old=Sys.time(); newt = old
# for (i in 1:100 ) {
# for (i in 101:200 ) {
# for (i in 201:ncol(otu.train) ) {
for (i in 1:ncol(otu.train) ) {
	old=Sys.time()
	
	custom_predict <- function(model1,new_data) {
		newdd = select(new_data, -'otu.train...i.', -'UTM_E',-'UTM_N')
		spdd = select(new_data, 'UTM_E','UTM_N')
#		print(c(dim(newdd),dim(spdd)))
		apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=spdd)) , along = -1L), 2:3, mean)[,i]
}
	
	fl = flashlight::flashlight(model = model.train, data=data.frame(scale.env.train,XY.train, otu.train[,i]), y = "otu.train...i.", label = "oregon", predict_function = custom_predict)
	load( here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
	vv = imp$data$variable[sapply(1:10, function(i) which(abs(imp$data$value)==sort(abs(imp$data$value),decreasing=T)[i]) ) ]
	int <- flashlight::light_interaction(fl, pairwise=F, type='H', v=vv, grid_size = round(nrow(otu.train)*.5), n_max = round(nrow(otu.train)*.9),seed = 42)
	? normalize = F, 
#	grid_size = 30, n_max = 50,
	newt = Sys.time() - old
	print(c(i,newt),digit=3)
	
#	plot(light_ice(fl,v="B5_20180717"))
#	plot(int)
	save(int, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	rm(int)
}
	

```                                       


```{r plot-individual}
mm = c('vint','vind'); i=1; mm=mm[i]; mm
	
# ... choose species ...  
# prevalent spp
taxadd = data.frame(sum=colSums(otu.train>0), otu=names(otu.train))
taxadd = taxadd[order(taxadd$sum, decreasing=T),]
taxadd$sum.seq = 1:nrow(taxadd)
str(taxadd)
	
spp.list = taxadd[c(sample(1:10,1),sample.int((nrow(taxadd)*.45):(nrow(taxadd)*.55),1),sample.int((nrow(taxadd)*.9):nrow(taxadd),1 )), ]  
spp.list
	
plot.list=list()
for (j in 1:3) {
	pp=ggplot()
	ind = which(names(otu.train)==spp.list$otu[j])
	if (mm == 'vint') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{ind}.rdata' ) ) )
		pp = plot(int, fill = "darkred")}
	if (mm == 'vind') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{ind}.rdata' ) ))
		pp = plot(imp, fill = "darkred")}
	
	name = toString(strsplit(strsplit(as.character(spp.list$otu[j]),'__')[[1]][2],'_')[[1]][2:6])
	plot.list[[j]] = pp + ggtitle(paste0(spp.list$sum[j],', ', name)) + 
    theme(plot.title = element_text(size = 9))
}
	 
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}{preocc}_{date.model.run}.pdf')), width=12, height=7)
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('fl-overall-int_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}{preocc}_{date.model.run}.pdf')), width=12, height=5)
	
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]], nrow=1)
	
dev.off()
	
```


```{r data-allspp-imp}
wxy = ''		# '' , 'coord'
mm = c('vint','vind'); name='a'; i=2; mm=mm[i]; mm
	
varx = data.frame(vardd = c(names(scale.env.train),names(XY.train)), ind = 1:(ncol(scale.env.train)+ncol(XY.train)))
impdd = data.frame(out=character(), var.abs=character(), varx.abs=numeric(), value.abs=numeric(), var.pos=character(),varx.pos=numeric(), value.pos=numeric())
	
for (i in 1:ncol(otu.train)) {
	print(i)
#	tryCatch({
	xpos='a'; xabs='a'; idd='a'
	if (mm == 'vint') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
		loadd = int; rm(int)
		
		xabs = most_important(loadd, 1); n1=1
		if (wxy=='' & startsWith(xabs, 'UTM_')) { while (startsWith(xabs, 'UTM_')) {n1=n1+1; xabs=most_important(loadd, n1)[n1]; print(xabs)} }
		n1=n1+1; xpos = most_important(loadd, n1)[n1]
		if (wxy=='' & startsWith(xpos, 'UTM_')) { while (startsWith(xpos, 'UTM_')) {n1=n1+1; xpos=most_important(loadd, n1)[n1]; print(xpos)} }
		
		idd = data.frame(otu=names(otu.train)[i], var.abs=xabs, varx.abs=varx$ind[varx$vardd==xabs], value.abs=loadd$data$value[loadd$data$variable==xabs], var.pos=xpos,varx.pos=varx$ind[varx$vardd==xpos], value.pos=loadd$data$value[loadd$data$variable==xpos])
		name = glue( 'overall-int-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
		}
	if (mm == 'vind') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ))
		loadd = imp; rm(imp)
		
		n1=1; xabs = which(abs(loadd$data$value)==max(abs(loadd$data$value))); vabs = loadd$data$variable[xabs]
		
		if (wxy=='' & startsWith(vabs, 'UTM_')) { while (startsWith(vabs, 'UTM_')) {n1=n1+1; xabs=which(abs(loadd$data$value)==sort(abs(loadd$data$value),decreasing=T)[n1]); vabs = loadd$data$variable[xabs]; print(vabs) } }
		n1=1; xpos = which(loadd$data$value==max(loadd$data$value)); vpos = loadd$data$variable[xpos]
		if (wxy=='' & startsWith(vpos, 'UTM_')) { while (startsWith(vpos, 'UTM_')) {n1=n1+1; xpos=which(loadd$data$value==sort(loadd$data$value,decreasing=T)[n1]); vpos = loadd$data$variable[xpos]; print(vpos) } }
		
		idd = data.frame(otu=names(otu.train)[i], var.abs=loadd$data$variable[xabs], varx.abs=varx$ind[varx$vardd==loadd$data$variable[xabs]], value.abs=loadd$data$value[xabs], var.pos=loadd$data$variable[xpos],varx.pos=varx$ind[varx$vardd==loadd$data$variable[xpos]], value.pos=loadd$data$value[xpos])
		name = glue( 'imp-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
		}
	
	impdd = rbind(impdd, idd)
	rm(idd)
	
	
	write.table(impdd, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name ), row.names=F,sep=',')
#	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
	
str(impdd)
	
name = glue( 'imp-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
if (mm == 'vint') { name = glue( 'overall-int-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' ) }
name
	
impdd = read.table(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name), header=T,sep=',')
	

```


```{r plot-circular}
table(impdd$otu == names(otu.train)); wxy; mm
# variables index & value which has biggest coefficient
effect_comb = data.frame(max_effects=impdd$varx.abs, V2=impdd$value.abs)
if (nrow(effect_comb)<ncol(otu.train)) { print('caution') }
str(effect_comb)
	
effect_comb2 = data.frame(max_effects=impdd$varx.pos, V2=impdd$value.pos)
str(effect_comb2)
	
otu.text = glue("incidence ({abund})")
version.text = ''
names(scale.env.train)
	
evnames = c("be10",'HJA',"tri","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple')
if (formula.env=='oldVars') {
	evnames = c('HJA','ele','canopy','minT','preci','road','stream','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','2m','2_4m','4_16m','p25','p95','rumple')}
if (formula.env=='vars2') {
	evnames = c("be10","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple','HJA')}
if (formula.env=='vars3') {
	evnames = c("be10", "slope","tri", "Nss", "Ess",'ndmi.std',"ndviP5", "ndviP50", "ndviP95", "ndmiP5", "ndmiP50", "ndmiP95","saviP50",'lcB1','lcB3','lcB4','lcB5','lcB7','lcB10','ndmi.std100m',"ndviP5.100m","ndviP50.100m","ndviP95.100m","ndmiP5.100m","ndmiP50.100m","ndmiP95.100m","saviP50.100m",'lcB1.100m','lcB3.100m','lcB4.100m','lcB5.100m','lcB7.100m','lcB10.100m',"tpi250",  "tpi500", "tpi1k" , "ht", "ht.250", "ht.1k", "2_4.250", "2_4.1k", "4_16", "4_16.250", "4_16.1k", "mTopo","cut1k",'stream', 'road','disturb','p25','rumple','HJA' )
	}
	
if (wxy=='coord') { evnames=append(evnames, c('UTM_E','UTM_N')) }; evnames
	

# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI',glue('{formula.env}'), glue('var{wxy}-imp-pos-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
	
cov.circle.env(version.text, evnames, otu.text, effect_comb=effect_comb2, otu.tbl=otu.train) 
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('var{wxy}-imp-overall-int-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('var{wxy}-imp-abs-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
	
cov.circle.env(version.text, evnames, otu.text, effect_comb=effect_comb, otu.tbl=otu.train)
	
dev.off()
	

```


```{r extract-max-vars}
wxy; mm; ddlist=list(); numvar=10
# wxy = ''; mm = 'vind'
	
for (i in 1:2) {
	aaa = c('vind','vint'); mm=aaa[i]
	if (mm=='vind') {numvar=25}; if (mm=='vint') {numvar=10}
	name = glue( 'imp-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
	if (mm == 'vint') { name = glue( 'overall-int2-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' ) }
	print(c(numvar, name)) 
	
	impdd = read.table(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name), header=T,sep=',')
	
	effect_comb = data.frame(max_effects=impdd$varx.abs, value=impdd$value.abs, var=impdd$var.abs, otu=impdd$otu)
	if (nrow(effect_comb)<ncol(otu.train)) { print('caution') }
#	str(effect_comb)
	
	dd = left_join(effect_comb, (effect_comb %>% count(max_effects)), by=c('max_effects'='max_effects')) %>% add_column(., type=mm)
	aaa = effect_comb %>% count(max_effects); aaa = aaa[order(aaa$n, decreasing=T),]; aaa = aaa$max_effects[1:numvar]
	aaa = unlist(sapply(1:numvar, function(x) {which(dd$max_effects==aaa[x])})); dd = dd[aaa,]
	dd$var2 = evnames[dd$max_effects]
	print(c(mm, i, length(unique(dd$max_effects)))); ddlist[[i]]=dd
	
} 
	
str(ddlist)
unique(ddlist[[1]]$var); unique(ddlist[[2]]$var)
ddlist[[3]] = unique(c(as.vector(unique(ddlist[[1]]$var)), as.vector(unique(ddlist[[2]]$var)))); ddlist[[3]] 
save(ddlist, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', glue('data-var{wxy}-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}')))
	
g1=ggplot(ddlist[[1]], aes(x=var2)) + geom_bar(stat="count") + scale_x_discrete(limits=unique(ddlist[[1]]$var2)) + xlab('individual') + theme(axis.text=element_text(size=6.5))
g2=ggplot(ddlist[[2]], aes(x=var2)) + geom_bar(stat="count") + scale_x_discrete(limits=unique(ddlist[[2]]$var2)) + xlab('interaction')
# , y=n  stat="identity"
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('var{wxy}-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=10)
grid.arrange(g1,g2, nrow=2)			# , widths=c(.44,.56)
	
dev.off()
	

```

```{r extract-3-maxV}
evnames
wxy = ''		# '' , 'coord'
mm = c('vint','vind'); name='a'; i=2; mm=mm[i]; mm
	
varx = data.frame(vardd = c(names(scale.env.train),names(XY.train)), ind = 1:(ncol(scale.env.train)+ncol(XY.train)))
impdd = data.frame(out=character(), var.abs=character(), varx.abs=numeric(), value.abs=numeric(), order=numeric())
	
for (i in 1:ncol(otu.train)) {
	print(i)
#	tryCatch({
	xpos='a'; xabs='a'; idd='a'
	if (mm == 'vint') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction2',glue( 'fl-overall-int_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
		loadd = int; rm(int)
		
		xabs = most_important(loadd, 1); n1=1
		if (wxy=='' & startsWith(xabs, 'UTM_')) { while (startsWith(xabs, 'UTM_')) {n1=n1+1; xabs=most_important(loadd, n1)[n1]; print(xabs)} }
		n1=n1+1; xpos = most_important(loadd, n1)[n1]
		if (wxy=='' & startsWith(xpos, 'UTM_')) { while (startsWith(xpos, 'UTM_')) {n1=n1+1; xpos=most_important(loadd, n1)[n1]; print(xpos)} }
		
		idd = data.frame(otu=names(otu.train)[i], var.abs=xabs, varx.abs=varx$ind[varx$vardd==xabs], value.abs=loadd$data$value[loadd$data$variable==xabs], var.pos=xpos,varx.pos=varx$ind[varx$vardd==xpos], value.pos=loadd$data$value[loadd$data$variable==xpos])
		name = glue( 'overall-int-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
		}
	if (mm == 'vind') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ))
		loadd = imp; rm(imp)
		
		values=sort(abs(loadd$data$value),decreasing=T)[1:5]
		vars = vector(); n1=1 
		idd=data.frame()
		for (j in 1:5) {
			xabs = which(abs(loadd$data$value)==values[j]); vabs = loadd$data$variable[xabs]
			idd=rbind(idd, data.frame(otu=names(otu.train)[i], var.abs=vabs, varx.abs=varx$ind[varx$vardd==vabs], value.abs=loadd$data$value[xabs], order=j))
			vars[j]=vabs
		}
		utm = startsWith(vars, 'UTM_')
		ccc=vector()
		for (j in 1:5) {
			if (utm[j]) { print(idd$var.abs[j]); ; ccc=append(ccc,j) }
		}
		if (length(ccc)>0) {idd=idd[-as.numeric(ccc),]}; rm(ccc)
		if (nrow(idd)!=3) { idd = idd[1:3,] }
#		most_important(loadd, 2)
		
		name = glue( 'imp-var{wxy}3-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
	}
	
	impdd = rbind(impdd, idd)
	rm(idd)
	
	write.table(impdd, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name ), row.names=F,sep=',')
#	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
str(impdd)
	
wxy; mm; ddlist=list(); numvar=10
	
for (i in 1:3) {
	mm='vind'
	if (mm=='vind') {numvar=10}; if (mm=='vint') {numvar=10}
	name = glue( 'imp-var{wxy}3-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' )
	if (mm == 'vint') { name = glue( 'overall-int2-var{wxy}-fl_min{minocc}{preocc}_{abund}_{cv}_{formula.env}.csv' ) }
	print(c(numvar, name)) 
	
	impdd = read.table(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name), header=T,sep=',')
	
	impdd=subset(impdd, order==i)
	effect_comb = data.frame(max_effects=impdd$varx.abs, value=impdd$value.abs, var=impdd$var.abs, otu=impdd$otu, order=impdd$order)
	if (nrow(effect_comb)<ncol(otu.train)) { print('caution') }
#	str(effect_comb)
	
	dd = left_join(effect_comb, (effect_comb %>% count(max_effects)), by=c('max_effects'='max_effects')) %>% add_column(., type=mm)
	aaa = effect_comb %>% count(max_effects); aaa = aaa[order(aaa$n, decreasing=T),]; aaa = aaa$max_effects[1:numvar]
	aaa = unlist(sapply(1:numvar, function(x) {which(dd$max_effects==aaa[x])})); dd = dd[aaa,]
	dd$var2 = evnames[dd$max_effects]
	print(c(mm, i, length(unique(dd$max_effects)))); ddlist[[i]]=dd
	ddlist[[i]] = dd
	rm(dd,impdd,effect_comb,aaa)
} 
	
str(ddlist)
unique(ddlist[[1]]$var); unique(ddlist[[2]]$var); unique(ddlist[[3]]$var)
ddlist[[4]] = unique(c(as.vector(unique(ddlist[[1]]$var)), as.vector(unique(ddlist[[2]]$var)), as.vector(unique(ddlist[[3]]$var)))); ddlist[[4]]
save(ddlist, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', glue('data-var{wxy}3-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}')))
	
dd = rbind(select(ddlist[[1]], var2, type, order, otu), select(ddlist[[2]], var2, type, order, otu), select(ddlist[[3]], var2, type, order, otu)); dd$order=as.character(dd$order); str(dd)
	
load( here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', glue('data-var{wxy}-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}')) )
dint = select(ddlist[[2]], var2, type, otu) %>% add_column(order=ddlist[[2]]$type, .before='otu')
dd = rbind(dd, dint)
	
g1=ggplot(dd[dd$order=='vint',], aes(x=var2)) + geom_bar(stat='count') + scale_x_discrete(limits=unique(dd$var2[dd$order=='vint'])) + xlab('interaction') #+ theme(axis.text=element_text(size=6.5))
	
g2 = ggplot(dd[dd$order!='vint',], aes(x=var2, fill=order)) + geom_bar(position="dodge") + xlab('individual') + scale_x_discrete(limits=unique(dd$var2[dd$order!='vint'])) + theme(axis.text=element_text(size=6), legend.position=c(.9,.8)) 
# , y=n  stat="identity" 
	
a = dd %>% count(var2); a = a[order(a$n,decreasing=T),] 
g3 = ggplot(dd, aes(x=var2, colour=order)) + geom_bar(stat="count") + xlab('individual') + theme(axis.text=element_text(size=6), legend.position=c(.9,.8)) + scale_x_discrete(limits=a$var2)
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('var{wxy}3-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=12, width=10)
grid.arrange(g2,g1,g3, nrow=3)			# , widths=c(.44,.56)
	
dev.off()
	

```


```{r table-max-vars}
dind=data.frame(); dint=data.frame(); dsum=data.frame()
for (i in c('qp','pa')) {
	jabund=i
	load( here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', glue('data-var{wxy}-imp-sum_{formula.env}_{jabund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}')) )
	dind=rbind(dind,ddlist[[1]] %>% add_column(dd=i))
	dint=rbind(dint,ddlist[[2]] %>% add_column(dd=i))
	dsum=rbind(dsum,data.frame(ddlist[[3]]) %>% add_column(dd=i, minocc=minocc, CV=cv, Svars=formula.env) %>% rename(var=1))
}
str(dind); str(dint); dsum
	
write.table(dsum, file=here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),'selected-var.csv'), row.names=F, sep=',')
	
g1 = ggplot(dind, aes(x=var2, colour=dd)) + geom_bar(stat="count") + scale_x_discrete(limits=unique(dind$var2)) + xlab('individual') + theme(axis.text=element_text(size=6.2))
g2 = ggplot(dint, aes(x=var2, colour=dd)) + geom_bar(stat="count") + scale_x_discrete(limits=unique(dint$var2)) + xlab('interaction') + theme(axis.text=element_text(size=8))
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('{formula.env}'),glue('var{wxy}-imp-sum_{formula.env}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=14)
grid.arrange(g1,g2, nrow=2)			# , widths=c(.44,.56)
	
dev.off()
	
names(ddlist)
formula.env = c('oldVars','newVars','vars2')
abund=c('pa','qp');i=2;abund=abund[i]
date.model.run = c('20210221','20210310'); i=2; date.model.run = date.model.run[i]; 
	
evnames = c("be10",'HJA',"tri","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple')
if (formula.env=='oldVars') {
	evnames = c('HJA','ele','canopy','minT','preci','road','stream','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','2m','2_4m','4_16m','p25','p95','rumple')}
if (formula.env=='vars2') {
	evnames = c("be10","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple','HJA')}
	
load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', glue('data-var{wxy}-imp-sum_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}')))
	
c(abund, cv, minocc, formula.env)
evnames
ddlist[[3]]
	


```


```{r data-for-boxplot}
numvar = 10
chs = c('abs','pos'); i=1
choice = chs[i]
# variables index & value which has biggest coefficient
effect_comb=data.frame()
if (choice=='abs') {effect_comb = data.frame(max_effects=impdd$varx.abs, V2=impdd$value.abs, var=impdd$var.abs, otu=impdd$otu)}
if (choice=='pos') {effect_comb = data.frame(max_effects=impdd$varx.pos, V2=impdd$value.pos, var=impdd$var.pos, otu=impdd$otu)}
str(effect_comb)
	
dd = left_join(effect_comb, (effect_comb %>% count(max_effects)), by=c('max_effects'='max_effects'))
aaa = effect_comb %>% count(max_effects); aaa = aaa[order(aaa$n, decreasing=T),]; aaa = aaa$n[numvar]
dd = dd %>% filter(n>=aaa); length(unique(dd$max_effects))
if(length(unique(dd$max_effects))!=numvar) {
	i=length(unique(dd$max_effects)) - numvar
	for (ii in 1:i) { dd = dd %>% filter(max_effects!=unique(dd$max_effects[dd$n==min(dd$n)])[ii]) }
} 
str(dd)
	
# .... prevalent spp
taxadd = data.frame(sum=colSums(otu.train>0), otu=names(otu.train))
taxadd = taxadd[order(taxadd$sum, decreasing=T),]
taxadd$sum.seq = 1:nrow(taxadd)
str(taxadd)
# 352, 3 -> min5 ; 303, 3 -> min6
	
# ... taxonomy 
taxadd$order = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
taxadd$class = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
taxadd$family = sapply(strsplit(sapply(str_split(taxadd$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
taxadd = left_join(taxadd, (taxadd %>% count(order)), by=c('order'='order')) %>% rename(nOrder=n)
str(taxadd)
	
dd = left_join(dd, taxadd, by=c('otu'='otu'))
str(dd)
	

# ................... load-auc-data ..........................
set = 'test' # 'explain', test
roc.dd = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, glue('{formula.env}_{date.model.run}'), set, glue('roc_result_{set}_{period}_{trap}_{abund}_{cv}_min{minocc}_{formula.env}_{metric}_{date.model.run}.RDS')))
names(roc.dd)
dim(roc.dd$otu)
	
# ... make long table 
auc.all = data.frame(otu=as.character(names(otu.train)), auc.test=rep(.1,length=ncol(otu.train)), auc.exp=rep(.1,length=ncol(otu.train)))
str(auc.all)
	
formula = paste0(abund,'.',formula.env,'-',strsplit(metric,'.test')[[1]][1]) 
	
b.t = data.frame( auc.test = roc.dd$roc.allS, otu.t = names(roc.dd$otu) )
str(b.t)
	
auc.all = left_join(auc.all, b.t, by=c('otu'='otu.t'), suffix=c('', glue('.{formula}')), copy=T)
str(auc.all)
	
auc.all = dplyr::select(auc.all, -'auc.test',-'auc.exp')
	
# ... extract taxonomy info
auc.all = right_join(auc.all, dd, by=c('otu'='otu'))
abc = data.frame(seq.var=letters[1:length(unique(dd$n))], var=sort(unique(dd$n),decreasing=T))
auc.all$oVar = sapply(1:nrow(auc.all), function(x) paste(abc$seq.var[abc$var==auc.all$n[x]],'.',auc.all$var[x],'.',auc.all$n[x], sep=''))
names(auc.all)
	
```


```{r plot-box}
abund; minocc; formula.env; cv
plot.date=20210309
choice; numvar
dd = auc.all %>% rename(auc.test=2)
names(dd)
	
cc = dd %>% count(var)
cc = cc[order(cc$n, decreasing=T),]
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('box_overall-int_max{numvar}Vars_{period}_{trap}_{abund}_{cv}_min{minocc}_{formula.env}_{metric}_{date.model.run}.pdf')), width=12, height=6)
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('box_{choice}_max{numvar}Vars_{period}_{trap}_{abund}_{cv}_min{minocc}_{formula.env}_{metric}_{date.model.run}.pdf')), width=12, height=6)
# glue('{formula.env}_{date.model.run}'),
	
ggplot(dd, aes(x=var, y=V2)) + geom_jitter(width = 0.35, aes(colour = order, size=auc.test)) + geom_boxplot(width=0.1) + ylab('max_effect') + scale_x_discrete(limit=cc$var) + annotate(geom="text", y=min(dd$V2)*1.1, x=1:length(cc$n), label=as.character(cc$n), col='red') + theme_minimal() + annotate(geom="text", x=numvar,y=max(dd$V2), label=glue('AUC.test: {round(mean(dd$auc.test,na.rm=T),3)}'))
#  + geom_hline(yintercept=.75, col='gray') + theme(legend.position = c(.95,0.3) ) + guides(color = FALSE)
	
dev.off()
	

```



