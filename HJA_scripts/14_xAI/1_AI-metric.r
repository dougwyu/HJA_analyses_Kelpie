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
date.model.run = 20210221   # according to tuning result
abund = 'pa'		# 'qp','pa'  !!! change accordingly
formula.env = 'newVars'		# 'envDNN.newVars', 'changed.envDNN'
minocc = 6
cv = '10CV'
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
```


```{r load-direct-data}
# ... read data ...
abund; formula.env; minocc; cv
nstep=1000
datapath=glue('../12_sjsdm_general_model_outputs')
	
# ...... from rdata ........
load(here(datapath,'source', glue('forada_data_{period}_m1m2_min{minocc}_{date.model.run}_{formula.env}.rdata')))
tuning=read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune',glue('{formula.env}_{date.model.run}'),'best', glue('best_manual_tuning_sjsdm_{cv}_{trap}{period}_mean_AUC_{abund}_min_{minocc}_nSteps_{nstep}.csv')), header=T, sep=',') 
	
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
metric; names(tuning)
	
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS',glue('{formula.env}_{date.model.run}'), glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_{cv}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{lambda.bio[i]}_{drop[i]}_{date.model.run}.RDS')) )
	
plot(model.train$history)
	

```


```{r demo-metrics} 
# flashlight package
# for (i in 1:100 ) {
# for (i in 101:200 ) {
# for (i in 201:ncol(otu.train) ) {
for (i in 1:ncol(otu.train) ) {
	print(i)
	custom_predict <- function(model1,new_data) {
		newdd = select(new_data, -'otu.train...i.')
		apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=XY.train)) , along = -1L), 2:3, mean)[,i]
}
	
	fl = flashlight::flashlight(model = model.train, data=data.frame(scale.env.train,otu.train[,i]), y = "otu.train...i.", label = "oregon", predict_function = custom_predict)
	
	imp <- flashlight::light_importance(fl, m_repetitions = 3,type = "permutation")
	save(imp, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
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
	load( here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
	vv = imp$data$variable[sapply(1:10, function(i) which(abs(imp$data$value)==sort(abs(imp$data$value),decreasing=T)[i]) ) ]
	int <- flashlight::light_interaction(fl, pairwise=F, type='H', v=vv, grid_size = 30, n_max = 50, seed = 42)
	newt = Sys.time() - old
	print(c(i,newt),digit=3)
	
#	plot(light_ice(fl,v="B5_20180717"))
#	plot(int)
	save(int, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction',glue( 'fl-overall-int_min{minocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
	
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
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction',glue( 'fl-overall-int_min{minocc}_{abund}_{cv}_{formula.env}_spp{ind}.rdata' ) ) )
		pp = plot(int, fill = "darkred")}
	if (mm == 'vind') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{cv}_{formula.env}_spp{ind}.rdata' ) ))
		pp = plot(imp, fill = "darkred")}
	
	name = toString(strsplit(strsplit(as.character(spp.list$otu[j]),'__')[[1]][2],'_')[[1]][2:6])
	plot.list[[j]] = pp + ggtitle(paste0(spp.list$sum[j],', ', name)) + 
    theme(plot.title = element_text(size = 9))
}
	 
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), width=12, height=5)
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('fl-overall-int_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), width=12, height=5)
	
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]], nrow=1)
	
dev.off()
	
```


```{r data-allspp-imp}
mm = c('vint','vind'); name='a'; i=1; mm=mm[i]; mm
	
varx = data.frame(vardd = names(scale.env.train), ind = 1:ncol(scale.env.train))
impdd = data.frame(out=character(), var.abs=character(), varx.abs=numeric(), value.abs=numeric(), var.pos=character(),varx.pos=numeric(), value.pos=numeric())
	
for (i in 1:ncol(otu.train)) {
	print(i)
#	tryCatch({
	xpos='a'; xabs='a'; idd='a'
	if (mm == 'vint') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', 'interaction',glue( 'fl-overall-int_min{minocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ) )
		loadd = int; rm(int)
		xabs = most_important(loadd, 1); xpos = most_important(loadd, 2)[2]
		idd = data.frame(otu=names(otu.train)[i], var.abs=xabs, varx.abs=varx$ind[varx$vardd==xabs], value.abs=loadd$data$value[loadd$data$variable==xabs], var.pos=xpos,varx.pos=varx$ind[varx$vardd==xpos], value.pos=loadd$data$value[loadd$data$variable==xpos])
		name = glue( 'overall-int-var-fl_min{minocc}_{abund}_{cv}_{formula.env}.csv' )}
	if (mm == 'vind') {
		load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{cv}_{formula.env}_spp{i}.rdata' ) ))
		loadd = imp; rm(imp)
		xabs = which(abs(loadd$data$value)==max(abs(loadd$data$value))); xpos = which(loadd$data$value==max(oadd$data$value))
		idd = data.frame(otu=names(otu.train)[i], var.abs=imp$data$variable[xabs], varx.abs=varx$ind[varx$vardd==imp$data$variable[xabs]], value.abs=imp$data$value[xabs], var.pos=imp$data$variable[xpos],varx.pos=varx$ind[varx$vardd==imp$data$variable[xpos]], value.pos=imp$data$value[xpos])
		name = glue( 'imp-var-fl_min{minocc}_{abund}_{cv}_{formula.env}.csv' )}
	
	impdd = rbind(impdd, idd)
	rm(idd)
	
	
	write.table(impdd, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name ), row.names=F,sep=',')
#	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
str(impdd)
	
impdd = read.table(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight','sum-flash', name), header=T,sep=',')
	

```
#> model.old=sjSDM(Y = s.otu.train,
#+   env = linear(data=scale.env.train, formula = ~.,
#+   lambda = lambda.env, alpha = alpha.env),
#+    biotic = bioticStruct(lambda=lambda.env, alpha=alpha.env, on_diag=F, inverse = FALSE),
#+ spatial = linear(data=XY.train, ~0+UTM_E*UTM_N,lambda = lambda.sp , alpha = alpha.sp),
#+ learning_rate = 0.003, # 0.003 recommended for high species number 
#+   step_size = NULL, iter =1,family=stats::binomial('probit'), sampling =1)
#model = model.old
#result = list( beta = coef(model), sigma = getCov(model),loss=NULL,history=model$history, p=NULL, logLik=logLik(model))

```{r plot-circular}
impdd$otu == names(otu.train)
# variables index & value which has biggest coefficient
effect_comb = data.frame(max_effects=impdd$varx.abs, V2=impdd$value.abs)
if (nrow(effect_comb)<ncol(otu.train)) {}
str(effect_comb)
	
effect_comb2 = data.frame(max_effects=impdd$varx.pos, V2=impdd$value.pos)
str(effect_comb2)
	
otu.text = glue("incidence ({abund})")
version.text = ''
names(scale.env.train)
	
evnames = c("be10",'HJA',"tri","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple')
if (formula.env=='oldVars') {evnames = c('HJA','ele','canopy','minT','preci','road','stream','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','2m','2_4m','4_16m','p25','p95','rumple')}
	

# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('var-imp-pos-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
	
cov.circle.env(version.text, evnames, otu.text, effect_comb=effect_comb2, otu.tbl=otu.train) 
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('var-imp-overall-int-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('var-imp-abs-flashlight-train_{formula.env}_{abund}_{cv}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
	
cov.circle.env(version.text, evnames, otu.text, effect_comb=effect_comb, otu.tbl=otu.train)
	
dev.off()
	

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



