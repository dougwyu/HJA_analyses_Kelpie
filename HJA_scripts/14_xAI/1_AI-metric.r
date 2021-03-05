# Mar 04, 2021
# last modified: Feb 21, 2021
# play with tuned model (manually tuning in ADA cluster, 20210215-0221)



```{r setup}
# setwd('/media/yuanheng/SD-64g3/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/14_xAI')
# set_here()
	
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
abund = 'qp'		# 'qp','pa'  !!! change accordingly
formula.env = 'newVars'		# 'envDNN.newVars', 'changed.envDNN'
minocc = 6
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
```


```{r load-direct-data}
# ... read data ...
abund; formula.env; minocc
nstep=1000
datapath=glue('../12_sjsdm_general_model_outputs')
	
# ...... from rdata ........
if (formula.env == "oldVars") {load(here(datapath,'source','yuanheng_mod_data.rdata'))}
if (formula.env == "newVars") {
	load(here(datapath,'source',glue('forada_data_{period}_m1m2_min{minocc}_{date.model.run}_{formula.env}.rdata')))
	tuning=read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune',glue('{formula.env}_{date.model.run}'), glue('best_manual_tuning_sjsdm_5CV_{trap}{period}_mean_AUC_{abund}_min_{minocc}_nSteps_{nstep}.csv')), header=T, sep=',')}
	
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
i=2
	
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS',glue('{formula.env}_{date.model.run}'), glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{lambda.bio[i]}_{drop[i]}_{date.model.run}.RDS')) )  
plot(model.train$history)
	

```

## explore vip package 
#pred_fun = function(model1) {    
#  apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=NULL, SP=NULL)) , along = -1L), 2:3, mean)
#}
	
#vip::vip(object = model.train, data=data.frame(XY.train,scale.env.train,otu.train[,1]), y = "otu.train...1.", method = "permute", metric = "auc", pred_wrapper = pred_fun, reference_class = 1)
	
#vi(model.train, method='firm')
## only works for certain types of models !!!
## doesn't work

```{r demo-metrics} 
# flashlight package
# 1:100; 101:200; 201:ncol(otu.train)
# for (i in 
for (i in 1:ncol(otu.train) ) {
	print(i)
	custom_predict <- function(model1,new_data) {
		newdd = select(new_data, -'otu.train...i.')
		apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=XY.train)) , along = -1L), 2:3, mean)[,i]
}
	
	fl = flashlight::flashlight(model = model.train, data=data.frame(scale.env.train,otu.train[,i]), y = "otu.train...i.", label = "oregon", predict_function = custom_predict)
	
	imp <- flashlight::light_importance(fl, m_repetitions = 3,type = "permutation")
	save(imp, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{formula.env}_spp{i}.rdata' ) ) )
	
#	load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{formula.env}_spp{i}.rdata' ) ))
	
#	plot(imp, fill = "darkred")
#	most_important(imp, 2)
	}
	

```                                       


```{r plot-individual}
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
	ind = which(names(otu.train)==spp.list$otu[j])
	load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{formula.env}_spp{ind}.rdata' ) ))
	pp = plot(imp, fill = "darkred")
	name = toString(strsplit(strsplit(as.character(spp.list$otu[j]),'__')[[1]][2],'_')[[1]][2:6])
	plot.list[[j]] = pp + ggtitle(paste0(spp.list$sum[j],', ', name)) + 
    theme(plot.title = element_text(size = 9))
}
	 
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('flashlight-train_{formula.env}_{abund}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), width=12, height=5)
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]], nrow=1)
	
dev.off()
	
```


```{r data-allspp-imp}
varx = data.frame(vardd = names(scale.env.train), ind = 1:ncol(scale.env.train))
impdd = data.frame(out=character(), var.abs=character(), varx.abs=numeric(), value.abs=numeric(), var.pos=character(),varx.pos=numeric(), value.pos=numeric())
	
for (i in 1:ncol(otu.train)) {
	tryCatch({
	load(here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'fl_min{minocc}_{abund}_{formula.env}_spp{i}.rdata' ) ))
	xabs = which(abs(imp$data$value)==max(abs(imp$data$value)))
	xpos = which(imp$data$value==max(imp$data$value))
#	most_important(imp, 1)
	idd = data.frame(otu=names(otu.train)[i], var.abs=imp$data$variable[xabs], varx.abs=varx$ind[varx$vardd==imp$data$variable[xabs]], value.abs=imp$data$value[xabs], var.pos=imp$data$variable[xpos],varx.pos=varx$ind[varx$vardd==imp$data$variable[xpos]], value.pos=imp$data$value[xpos])
	
	impdd = rbind(impdd, idd)
	rm(idd)
	
	write.table(impdd, file=here(outputpath,'prediction_outputs','sjsdm-model-RDS',sjsdmVfolder,glue('{formula.env}_{date.model.run}'),'flashlight', glue( 'imp-var-fl_min{minocc}_{abund}_{formula.env}.csv' ) ), row.names=F,sep=',')
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
str(impdd)
	

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
effect_comb = data.frame(max_effects=impdd$varx.abs, V2=impdd$value.abs)
# variables index & value which has biggest coefficient
str(effect_comb)
	
effect_comb2 = data.frame(max_effects=impdd$varx.pos, V2=impdd$value.pos)
	
otu.text = glue("incidence ({abund})")
version.text = ''
names(scale.env.train)
	
evnames = c("be10",'HJA',"tri","slope","Nss","Ess","ht","ht.250","ht.500","ht.1k","2_4","2_4.250","2_4.500","2_4.1k","4_16","4_16.250", "4_16.500","4_16.1k","be500","mTopo","cut1k",'minT', 'maxT','preci','stream', 'road','disturb','B1','B2','B3','B4','B5','B6','B7','B10','B11','NDVI','EVI','B','G','W','p25','rumple')
	

# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'xAI', glue('var-imp-pos-flashlight-train_{formula.env}_{abund}_{period}_{trap}_min{minocc}_{date.model.run}.pdf')), height=8, width=8)
	
cov.circle.env(version.text, evnames, otu.text, effect_comb2, otu.tbl=otu.train) 
	
dev.off()
	

```









