# Feb 02, 2021
# code on sjsdm model run with DNN (on laptop)



```{r setup}
# setwd('/media/yuanheng/SD-64g2/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm_general_model_outputs')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan','pROC','MLmetrics', 'gridExtra')
	
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
abund = 'qp'		# 'qp','pa'  !!! change accordingly
trap <- "M1"; period = "S1"
date.model.run = 20210202   # !!! change accordingly
	
formula.env = 'envDNN-added'			# 'envDNN-newgis'
	
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
# 124 Var
	
# roughness index
# from "Hmsc_CD/oregon_ada/data/demStats.rdata"
load(here('..','10_eo_data','topo_data.rdata'))
str(topo.df)
# cov2_4.r500, Ess, cut.r1k.pt, ht
	

```


```{r select subsets-1}
# ..... select trap, session .....
trap; period
	
alldata1 = subset(alldata, trap == "M1" & period == "S1")
dim(alldata1)
	
names(alldata1)[102:125]
	
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
	
```

```{r select subsets-2 }
envnames = c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist")
	
# envnames = c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")
# ori.env = select(left_join(select(alldata1, envnames, "SiteName"), select(topo.df, 'tri','cov2_4.r500', 'Ess', 'cut.r1k.pt', 'cov4_16.r250','siteName'), by=c('SiteName'='siteName')), -'SiteName') 
# add 'roughness' -> tri; cov2_4.r500, Ess, cut.r1k.pt, cov4_16.r250
	
ori.env = select(left_join(select(alldata1, envnames, "SiteName"), select(topo.df,-'twi'), by=c('SiteName'='siteName')), -'SiteName')
	
str(ori.env)
	
ori.env.train = ori.env[a,]
ori.env.test = ori.env[-a,]
str(ori.env.test)
	
ori.XY = select(alldata1, starts_with('UTM'))
ori.XY.train = ori.XY[a,]
ori.XY.test = ori.XY[-a,]
str(ori.XY.train)
	
# ... view data ...
par(mfrow=c(2,2))
hist(ori.env.train$elevation_f,xlim=c(1000,5500))
hist(ori.env$elevation_f)
	
hist(ori.env.train$cov2_4.r500,xlim=c(.01,.13))
hist(ori.env$cov2_4.r500)
	
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
# insideHJA factor !!!
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


```{r save-data}
s.otu.test = as.matrix(otu.test)
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	
s.otu.train = as.matrix(otu.train)
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
write.table(s.otu.train, file=here('source','new-gis-data', glue('otu.train.{abund}.csv')), row.names=F, col.names=F, sep=',')
write.table(s.otu.test, file=here('source','new-gis-data', glue('otu.test.{abund}.csv')), row.names=F, col.names=F, sep=',')
	
write.table(scale.env.train, file=here('source','new-gis-data', 'scale.env.train.csv'), row.names=F, sep=',')
write.table(scale.env.test, file=here('source','new-gis-data', 'scale.env.test.csv'), row.names=F, sep=',')
	
write.table(XY.train, file=here('source','new-gis-data', 'scale.XY.train.csv'), row.names=F, sep=',')
write.table(XY.test, file=here('source','new-gis-data', 'scale.XY.test.csv'), row.names=F, sep=',')
	
```


```{r read-data}
# ...... from rdata ........
if (formula.env == "envDNN-added") {load(here('source','yuanheng_mod_data_envDNN-added.rdata'))}
if (formula.env == "envDNN.newVars") {load(here('source','yuanheng_mod_data_newVars.rdata'))}
	
# save(s.otu.test, s.otu.train, scale.env.test, scale.env.train, XY.test, XY.train, file=here(datapath,'source','yuanheng_mod_data_envDNN-added.rdata'))
	
if (abund == 'pa')
{
	s.otu.train = (s.otu.train>0)*1
	s.otu.test = (s.otu.test>0)*1
}
	
str(scale.env.test); formula.env
str(s.otu.train); abund
	
otu.name = read.table(here('source','otu-names.csv'),header=T)
	
```


```{r load-data-newgis}
# ... read data ...
s.otu.train = as.matrix(read.table(here('source','new-gis-data', glue('otu.train.{abund}.csv')), header=F, sep=','))
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
scale.env.train = read.table(here('source','new-gis-data', 'scale.env.train.csv'), header=T, sep=',')
str(scale.env.train)
	
XY.train =read.table(here('source','new-gis-data', 'scale.XY.train.csv'), header=T, sep=',')
str(XY.train)
	
# ... test data
s.otu.test = as.matrix(read.table(here('source','new-gis-data',glue('otu.test.{abund}.csv')), header=F, sep=','))
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	
scale.env.test = read.table(here('source','new-gis-data','scale.env.test.csv'), header=T, sep=',')
str(scale.env.test)
	
XY.test =read.table(here('source','new-gis-data','scale.XY.test.csv'), header=T, sep=',')
str(XY.test)
	

```


```{r model-DNN-newgis}
date.model.run = 20210204 
# one model to see how long it takes in ADA

# make sure input of sjsdm are numeric matrix
s.otu.train = as.matrix(otu.train)
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
str(scale.env.train)
	
# ..... for test run .....
# set variables
formula.env = 'envDNN-newgis'
lambda.env = .1
alpha.env = .5
lambda.sp = .1
alpha.sp = .9 
hidden.sp = c(50L,50L,10L)
acti.sp = 'relu'
drop = .3
	
model.train = sjSDM(Y = s.otu.train,
			  env = DNN(data=scale.env.train, formula = ~.,
					  hidden=hidden.sp, lambda = lambda.env, alpha = alpha.env, activation=acti.sp, dropout=drop, bias = T),
 			  biotic = bioticStruct(lambda=lambda.sp, alpha=alpha.sp, on_diag=F, inverse = FALSE),
			  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp, alpha=alpha.sp),
			  learning_rate = 0.003, # 0.003 recommended for high species number 
			  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L 
			 )
	
names(model.train)
	
# saveRDS(model=model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS',  glue('s-jSDM_model_{period}_{trap}_{abund}_min{minocc}_20210204_{formula.env}.RDS')) )
	
# ..... for test run .....
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS', glue('s-jSDM_model_{period}_{trap}_{abund}_min{minocc}_20210204_{formula.env}.RDS')) )
	
# pdf(here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'DNN_tune', glue('model-history_model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.pdf')), width=5, height=5)
	
plot(model.train$history)
	
dev.off()
	
```


```{r model-DNN.env}
# make sure input of sjsdm are numeric matrix
str(s.otu.train)
str(scale.env.train)
	
# set variables
formula.env = 'envDNN-added'
lambda.env = seq(.01,.3, length.out=3)	# .1
alpha.env = seq(.8,1, length.out=3)		# .9
lambda.sp = seq(0.01,.3, length.out=3)	# .1
alpha.sp =  seq(0.4,.6, length.out=3)	# .5 
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.2,.4, length.out=3) # .3
sample.bio = seq(0,1,length.out=11)
	

lambda.envN = 2		# sub in ubuntu 1,2 (qp); (pa) 2
	
for (alpha.envN in 1:3) {
	for (lambda.spN in 1:3) {
		for (alpha.spN in 1:3) {
			for (dropN in 1:3) {
				for (hiddenN in 1:2) {
					
					lambda.bioN = sample(1:11,1)
					alpha.bioN = sample(1:11,1)
					print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN, hiddenN))
					
					model.train = sjSDM(Y = s.otu.train,
					  env = DNN(data=scale.env.train, formula = ~.,
					  hidden=hidden[[hiddenN]], lambda = lambda.env[lambda.envN], alpha = alpha.env[alpha.envN], activation=acti.sp, dropout=drop[dropN], bias=T),
					  
					  biotic = bioticStruct(lambda=sample.bio[lambda.bioN], alpha=sample.bio[alpha.bioN], on_diag=F, inverse = FALSE),
					  
					  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp[lambda.spN], alpha=alpha.sp[alpha.spN]),
					  
					  learning_rate = 0.003, # 0.003 recommended for high species number 
					  step_size = NULL, iter = 140L, family=stats::binomial('probit'), sampling = 900L # 150L, 5000L
					 )
					 saveRDS(list(model=model.train, random=data.frame('lambda.bioN'=lambda.bioN, 'alpha.bioN'=alpha.bioN)), here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS','added-gis', glue('s-jSDM_tuning.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
	
				}
			}
		}
	}
}
	

```


```{r manually-tune}
# 20210202 run on laptop
# set variables (from 'model-...' chunk)
formula.env = 'envDNN-added'
lambda.env = seq(.01,.3, length.out=3)	# .1
alpha.env = seq(.8,1, length.out=3)		# .9
lambda.sp = seq(0.01,.3, length.out=3)	# .1
alpha.sp =  seq(0.4,.6, length.out=3)	# .5 
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.2,.4, length.out=3) # .3
sample.bio = seq(0,1,length.out=11)
	
tuning.dd = data.frame(lambda.env = numeric(), alpha.env = numeric(),lambda.sp = numeric(), alpha.sp = numeric(), lambda.bio = numeric(), alpha.bio = numeric(), drop = numeric(), hidden = character(), loglike = numeric(), loss= numeric(), AUC.explain=numeric(), AUC.test=numeric(), prAUC.explain=numeric(), prAUC.test=numeric(), ll.explain=numeric(), ll.test=numeric(), nagel.explain=numeric(), nagel.test=numeric())
	
# ccc = 0
# hiddenN=1; lambda.envN=2; alpha.envN=2; lambda.spN=2; alpha.spN=2;dropN=2
	
for (hiddenN in 1:2) {
	for (lambda.envN in 1:3) {
		for (alpha.envN in 1:3) {
			for (lambda.spN in 1:3) {
				for (alpha.spN in 1:3) {
					for (dropN in 1:3) {
						
		tryCatch({
			model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS','added-gis', glue('s-jSDM_tuning.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
#			plot(model.train$model$history)

			lambda.bioN = model.train$random$lambda.bioN
			alpha.bioN = model.train$random$alpha.bioN
			
			tdd = data.frame(lambda.env = lambda.env[lambda.envN], alpha.env = alpha.env[alpha.envN], lambda.sp = lambda.sp[lambda.spN], alpha.sp = alpha.sp[alpha.spN], lambda.bio = sample.bio[lambda.bioN], alpha.bio = sample.bio[alpha.bioN], drop = drop[dropN], hidden = as.character(hiddenN), loglike = logLik(model.train$model), loss= model.train$model$history[length(model.train$model$history)], AUC.explain=.1, AUC.test=.1, prAUC.explain=.1, prAUC.test=.1,ll.explain=.1, ll.test=.1,nagel.explain=.1, nagel.test=.1)
			
			for (pred in 1:2) {
				# 1 -> 'test'
				newdd = scale.env.test ; newsp = XY.test; otudd = s.otu.test
				if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train }
				
				pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train$model, newdata=newdd, SP=newsp)) , along = -1L), 2:3, mean)
				attr(pred.dd, 'dimnames') = NULL
				
				if (pred==1) {
					otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
					pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
					otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
				}
				
				otudd.pa = (otudd>0)*1
				auc.mean = mean(sapply(1:ncol(otudd), function(i) as.numeric(pROC::roc(otudd.pa[,i], pred.dd[,i],direction='<',quiet=T)$auc)))
				str(pred)
				prauc=mean(sapply(1:ncol(otudd), function(i) MLmetrics::PRAUC(pred.dd[,i],otudd.pa[,i])))
#				pred = ROCR::prediction(pred.dd[,i],otudd.pa[,i],label.ordering = c(0,1))
#				precision = pred@tp[[1]]/(pred@tp[[1]]+pred@fp[[1]])
#				recall = pred@tp[[1]]/(pred@tp[[1]]+pred@fn[[1]])
#				plot(recall, precision)
				if (pred==2) {tdd$AUC.explain=auc.mean; tdd$prAUC.explain=prauc}
				if (pred==1) {tdd$AUC.test=auc.mean; tdd$prAUC.test=prauc}
				
				# ll & nagel spp in sample
				rsq = data.frame(ll=rep(.1, length.out=dim(pred.dd)[1]), nagel=rep(.1, length.out=dim(pred.dd)[1]))
				for (j in 1:dim(pred.dd)[1]) {
#					print(j)
					p = pred.dd[j,]
					y = otudd.pa[j,]
					loglikP = sum( log( p*y + (1-p)*(1-y) ) )
					loglikN = sum( log( mean(p)*y + (1-mean(p))*(1-y) ) )		 
			#		rsq$logrsq[j] = (loglikP-loglikN)/(1-loglikN)
					rsq$nagel[j] = (1-exp(2/length(p)*(loglikN-loglikP))) / (1-exp(2/length(p)*loglikN))
					rsq$ll[j] = loglikP
				}
				if (pred==2) {tdd$ll.explain=mean(rsq$ll); tdd$nagel.explain=mean(rsq$nagel)}
				if (pred==1) {tdd$ll.test=mean(rsq$ll); tdd$nagel.test=mean(rsq$nagel)}
				
			}
			
		tuning.dd = rbind(tuning.dd, tdd)
		print(c(hiddenN, lambda.envN, alpha.envN, lambda.spN, alpha.spN, lambda.bioN, alpha.bioN, dropN))
		
		rm(model.train, tdd)
#		ccc = ccc+1
		
		write.table(tuning.dd, file=here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('manualII_tuning.sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')), row.names=F, sep=',')
	
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		
} }}}}}
	
ccc
str(tuning.dd)
	

```


```{r analyze-tune}
formula.env; abund 
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
	
tuning.dd = read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('manualII_tuning.sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')), header=T, sep=',')
str(tuning.dd)
	
pp = c(which(tuning.dd$AUC.test==max(tuning.dd$AUC.test)), which(tuning.dd$ll.test==sort(tuning.dd$ll.test,decreasing=T)[1]), which(tuning.dd$prAUC.test==max(tuning.dd$prAUC.test)))
tuning.dd[pp,]
# nagel, ll gives same result
	
p = which(tuning.dd$ll.test>tuning.dd$ll.explain)
tuning.dd[p,]
	
# pdf(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('plot_tuning.sjsdm.laptop_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.pdf')), width=11, height=4.5)
	
#pp = which(tuning.dd$AUC.test==max(tuning.dd$AUC.test))
#tuning.dd[pp,]
	
par(mfrow=c(1,2), cex=.8)
plot(y=tuning.dd$AUC.test, x=tuning.dd$loglike, pch=3, col='blue',ylim=c(min(tuning.dd$AUC.test,tuning.dd$AUC.explain),max(tuning.dd$AUC.test,tuning.dd$AUC.explain)),xlab='log-Likelihood',ylab='AUC')
points(y=tuning.dd$AUC.explain, x=tuning.dd$loglike, pch=8)
points(y=tuning.dd$AUC.test[pp], x=tuning.dd$loglike[pp], pch=19, col='red')
points(y=tuning.dd$AUC.explain[pp], x=tuning.dd$loglike[pp], pch=8, col='red')
abline(v=tuning.dd$loglike[pp], lty=2, col='red')
legend('topright',pch=c(3,8),col=c('blue','black'),c(paste('auc.test,',round(tuning.dd[pp,]$AUC.test,3)),paste('auc.train,',round(tuning.dd[pp,]$AUC.explain,3))),bty='n')
	
plot(y=tuning.dd$ll.test, x=tuning.dd$AUC.test, pch=3, col='blue',ylim=c(min(tuning.dd$ll.test,tuning.dd$ll.explain),max(tuning.dd$ll.test,tuning.dd$ll.explain)), xlab='AUC.test',ylab='ll metrics')
points(y=tuning.dd$ll.explain, x=tuning.dd$AUC.test, pch=8)
points(y=tuning.dd$ll.test[pp2], x=tuning.dd$AUC.test[pp2], pch=3, col='red')
points(y=tuning.dd$ll.explain[pp2], x=tuning.dd$AUC.test[pp2], pch=8, col='red')
abline(v=tuning.dd$AUC.test[pp2], lty=2, col='red')
	
points(y=tuning.dd$ll.test[pp3], x=tuning.dd$AUC.test[pp3], pch=3, col='brown')
points(y=tuning.dd$ll.explain[pp3], x=tuning.dd$AUC.test[pp3], pch=8, col='brown')
abline(v=tuning.dd$AUC.test[pp3], lty=2, col='brown')
	
points(y=tuning.dd$ll.test[pp], x=tuning.dd$AUC.test[pp], pch=3, col='brown')
points(y=tuning.dd$ll.explain[pp], x=tuning.dd$AUC.test[pp], pch=8, col='brown')
abline(v=tuning.dd$AUC.test[pp], lty=2, col='brown')
	
dev.off()
	
# model with optimal parameters
lambda.env = tuning.dd[pp,'lambda.env']
alpha.env = tuning.dd[pp,'alpha.env']
lambda.sp = tuning.dd[pp,'lambda.sp']
alpha.sp =  tuning.dd[pp,'alpha.sp']
lambda.bio = tuning.dd[pp,'lambda.bio']
alpha.bio =  tuning.dd[pp,'alpha.bio']
hidden1 = tuning.dd[pp,'hidden']
acti.sp = 'relu'
drop = tuning.dd[pp,'drop']
	
str(s.otu.train)
str(scale.env.train)
	
i = 1
	
for (i in 2:3) {
model.train = sjSDM(Y = s.otu.train,
	  env = DNN(data=scale.env.train, formula = ~.,
	  hidden= hidden[[hidden1[i]]], lambda = lambda.env[i], alpha = alpha.env[i], activation=acti.sp, dropout=drop[i], bias = T),
	  
	  biotic = bioticStruct(lambda=lambda.bio[i], alpha=alpha.bio[i], on_diag=F, inverse = FALSE),
	  
	  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp[i], alpha=alpha.sp[i]),
	  
	  learning_rate = 0.003, # 0.003 recommended for high species number 
	  step_size = NULL, iter = 140L, family=stats::binomial('probit'), sampling = 2000L # 150L, 5000L
)
	 
# saveRDS(model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS', glue('s-jSDM_tuned.model.laptop_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS')) )
	
}
	
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS', glue('s-jSDM_tuned.model.laptop_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS')) )
	
# pdf(here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'DNN_tune', glue('model-history_tuned.model.laptop_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env}_{alpha.env}_{lambda.sp}_{alpha.bio}_{drop}_{date.model.run}.pdf')), width=5, height=5)
	
plot(model.train$history)
	
dev.off()
	
```



```{r prediction}
model1 = model.train
formula = formula.env  # 
	
# predictive AUC
newdd = scale.env.test
newsp = XY.test
otudd = s.otu.test
set = 'test'
	
# explanatory AUC
 newdd = NULL
 newsp = NULL
 otudd = s.otu.train
 set = 'explain'
	
pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata=newdd, SP=newsp)) , along = -1L), 2:3, mean)
	
str(pred.dd)
attr(pred.dd, 'dimnames') = NULL
	
otudd = as.data.frame(otudd)
names(otudd)=otu.name[,1]
names(otudd)[1:5]
	
# .. delete otus that only 0/1
# ..... (needed for testing) .....
dim(otudd)
	
otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
which(otudd[dim(otudd)[1],] == 0)
	
pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
str(pred.dd)
	
otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
dim(otudd)
	
# ..... (needed for testing) .....
otudd.pa = (otudd>0)*1
table(otudd.pa==otudd)
	
# calculate AUC
roc.dd = sapply(1:ncol(otudd), function(i) as.numeric(pROC::roc(otudd.pa[,i], pred.dd[,i],direction='<')$auc))
	
head(roc.dd)
	
auc.mean = mean(roc.dd)
formula; abund; auc.mean
	
# saveRDS(list(pred.Y=pred.dd, otu=otudd, roc.allS=roc.dd, auc.mean=auc.mean), here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, 'prediction', glue('roc_result_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS')))
	
```


```{r roc_result-for-regression}
# .. set names
# predictive
set = 'test' # 'explain', test
formula
	
explain = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, 'prediction', glue('roc_result_explain_S1_M1_pa_min5_envDNN_20210125.RDS'))) 
names(explain)
# needed for selecting data
	
# .... prevalent spp
a = data.frame(sum=colSums(explain$otu), otu=attr(colSums(explain$otu),'names'))
a = a[order(a$sum, decreasing=T),]
a$sum.seq = 1:dim(a)[1]
str(a)
# 268, 3
	
# ... taxonomy 
a$order = sapply(strsplit(sapply(str_split(a$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[2])
a = left_join(a, (a %>% count(order)), by=c('order'='order'))
unique(a$n)
str(a)
	
a$class = sapply(strsplit(sapply(str_split(a$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[1])
a$family = sapply(strsplit(sapply(str_split(a$otu, '__'), function(aa) aa[2]), '_'), function(aa) aa[3])
	
# ... load roc result ...
roc.dd = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, 'prediction', glue('roc_result_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS'))) 
names(roc.dd)
	
# ... make long table 
roc.dd1 = as.data.frame(bind_cols(roc.dd$otu %>% pivot_longer(everything(), names_to='otu.name', values_to = 'obs'), pred=as.vector(t(roc.dd$pred.Y)), auc=rep(roc.dd$roc.allS, nrow(roc.dd$otu))))
roc.dd1 = left_join(roc.dd1, a, by = c('otu.name'='otu'), copy=T)
str(roc.dd1)
# 18*243 = 4374
# 70*268 = 18760
	
roc.dd1$otu.name = as.factor(roc.dd1$otu.name)
	
# check on dataset
roc.dd$otu[1,1:6]; roc.dd$pred.Y[1,1:6]; roc.dd1[1:6,2:3]
	
roc.dd50 = roc.dd1[which(roc.dd1$otu.name %in% as.character(a$otu[1:50])),]
sort(unique(roc.dd50$otu.name)) == sort(a$otu[1:50])
	
as.character(sort(unique(roc.dd50$otu.name))) == as.character(sort(a$otu[1:50]))
	
```


```{r correlate-auc}
# .... (define) ....
auc.all = data.frame(otu=as.character(otu.name$otu), auc.test=rep(.1,length=nrow(otu.name)), auc.exp=rep(.1,length=nrow(otu.name)))
str(auc.all)
# dim[1] == 268!!!
	
formula = 'envDNN-added'		# 'envDNN'
i=3
	
roc.dd.t = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, 'prediction', glue('roc_result_test_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS'))) 
names(roc.dd.t)
	
roc.dd.e = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, 'prediction', glue('roc_result_explain_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS'))) 
names(roc.dd.e)
	
tuning.dd[pp,]
if (i==1) {formula = paste0('envDNN-added-','roc')} # roc , ll , pr
if (i==2) {formula = paste0('envDNN-added-','ll')}; if (i==3) {formula = paste0('envDNN-added-','pr')}
	
# ... make long table 
b.t = data.frame( auc.test = roc.dd.t$roc.allS, otu.t = names(roc.dd.t$otu) )
str(b.t)
	
b.e = data.frame( auc.exp = roc.dd.e$roc.allS, otu.e = names(roc.dd.e$otu) )
str(b.e)
	
auc.te = inner_join(b.t, b.e, by=c('otu.t'='otu.e'))
str(auc.te)
	
auc.all = left_join(auc.all, auc.te, by=c('otu'='otu.t'), suffix=c('', glue('.{formula}')), copy=T)
str(auc.all)
	
# .. after all variables are added
auc.all = dplyr::select(auc.all, -'auc.test', -'auc.exp')
	
# ... extract taxonomy info
auc.all = left_join(auc.all, select(a, 'otu','order','class','family','sum','n'), by=c('otu'='otu'))
abc = data.frame(seq.order=letters[1:length(unique(a$n))], order=sort(unique(a$n),decreasing=T))
auc.all$Oorder = sapply(1:dim(auc.all)[1], function(x) paste(abc$seq.order[abc$order==auc.all$n[x]],'.',auc.all$order[x],'.',auc.all$n[x], sep=''))
str(auc.all)
	
# ....... regression .......
set  # 'test' , 'explain'
formula ='envDNN-added'
	
# ... loop ...
setS = c('roc','ll', 'pr')  
	
dd.ggplot = vector(mode = "list", length = 3)
str(dd.ggplot)
	
plot.list=list(ggplot(),ggplot(),ggplot())
	
for (j in 1:3) {
	set = setS[j]  
#	formula = formulaS[i]  # 'envDNN-spAll' , 'spDNN' , 'noDNN'
	ii=j*2; jj=ii+1
	auc.1 = select(auc.all, ii, jj, 'sum', 'order', 'class', 'family', 'Oorder') %>% rename(auc.test=1, auc.exp=2, incident=3)
	auc.1 = na.omit(auc.1)
	str(auc.1)
	dd.ggplot[[j]] = auc.1
	
	gp = ggplot(dd.ggplot[[j]], aes(auc.exp, auc.test)) + geom_point(aes(colour=factor(Oorder), size=incident))+ scale_size(range = c(1, 7)) + scale_colour_manual(values=colorRampPalette(c('dodgerblue3','firebrick2','yellow'))(12)) + geom_smooth(method='lm', se=F, colour='gray') + geom_abline(slope=1,intercept=0, linetype = "dashed", colour='gray', size=1.5) + theme(panel.background=element_rect(fill='snow')) + ggtitle(glue('{set}, {formula}')) + geom_hline(yintercept=0.5, linetype = "dashed", colour='red') + geom_vline(xintercept=0.5, linetype = "dashed", colour='red') + xlim(0,1) +ylim(0,1)
	if (j==2) {plot.list[[j]] = gp + theme(legend.position='left')} else {plot.list[[j]] = gp + theme(legend.position='none')}
	
}
	
#m1 = lm(auc.qp ~ auc.pa, data=auc.1)
#summary(m1)
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction', glue('auc-correlation_{period}_{trap}_min{minocc}_{formula}-laptop_tuned_{date.model.run}.pdf')), width=18, height=6)
	
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]], nrow=1, widths=c(.3,.4,.3))   #, layout_matrix= rbind(c(1), c(2,3), c(4,5)), heights=c(.4,5,5)
	
dev.off()
	

```


```{r violin-auc}
dd = roc.dd1
names(dd)
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction', glue('{envvar}'), glue('violin_order_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_{date.cross.validation}.pdf')), width=9, height=4)
	
cc = a %>% count(order)
cc = cc[order(cc$n, decreasing=T),]
	
ggplot(dd, aes(x=order, y=auc)) + geom_hline(yintercept=.75, col='gray') + geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + geom_boxplot(width=0.1) + theme_minimal() + scale_x_discrete(limit=cc$order) + annotate(geom="text", y=-.43, x=1:length(cc$n), label=as.character(cc$n), col='red')
	
dev.off()
	
# .. choose prevalent group
sort(unique(dd$sum), decreasing=T)
sort(unique(dd$sum.seq), decreasing=F)
	
dd = subset(dd, sum.seq<=sort(unique(dd$sum.seq))[25])
str(dd)
	
g1.1 = ggplot(dd, aes(x=as.factor(sum.seq), y=c(pred, obs))) + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1) + geom_boxplot(width=0.1) + theme_minimal() + scale_x_discrete(limit=as.character(sort(unique(dd$sum.seq))), labels=a$sum[sort(unique(dd$sum.seq))] )
#xlim( as.character(dd$sum))
	
dd$auc[dd$sum.seq==7]
	
```


```{r AUC-regression}
set = 'explain' # 'test'
dd = roc.dd1
	
#set = 'test.50spp'
#dd = roc.dd50
	
str(dd)
	
par(mfrow=c(2,2))
hist(dd$sum)
hist(dd$sum)
hist(dd$auc)
	
hist(dd$obs)
hist(log1p(dd$obs[dd$obs!=0]))
	
plot(dd$sum, dd$obs)
	
#/family /order  class/order/family
lme1 = lme(auc~obs + sum, random= ~1 | class/order/family, method='ML', data = dd, control = lmeControl(opt='optim', msMaxIter = 200, msMaxEval = 500, msVerbose = T))
	
summary(lme1)
	
# pdf(here(outputpath,'prediction_outputs','sjsdm-graph', sjsdmVfolder, 'prediction', glue('lme-auc_{set}_{datafolder}_{abund}_{formula}.pdf')), height=12,width=17)
# _spDNN
	
par(mfrow=c(2,3),cex=1.05)
	
hist(dd$obs, cex.lab=1.2)
hist(dd$sum, cex.lab=1.2)
hist(dd$auc, cex.lab=1.2)
	
obs.seq = seq(min(dd$obs),max(dd$obs), length=100)
sum.seq = seq(min(dd$sum),max(dd$sum), length=100)
fit = summary(lme1)$tTable[[1]] + summary(lme1)$tTable[[2]]*obs.seq + summary(lme1)$tTable[[3]]*sum.seq
	
#plot(dd$sum, dd$auc, main='lme', cex.lab=1.2, pch=20)
#lines(sum.seq, fit)
	
dd1 = dd[seq(1, dim(dd)[1],by=2),]
	
plot(dd1$obs, dd1$auc, main='lme', cex.lab=1.2, pch=20)
lines(obs.seq, fit)
	
par.resid.obs = resid(lme1)[seq(1, dim(dd)[1],by=2)] + summary(lme1)$tTable[[2]]*dd1$obs
	
plot(dd1$obs, par.resid.obs, main='lme', cex.lab=1.2, ylab='AUC(partial.resid.obs)', pch=20)
lines(obs.seq, summary(lme1)$tTable[[2]]*obs.seq)
	
par.resid.sum = resid(lme1)[seq(1, dim(dd)[1],by=2)] + summary(lme1)$tTable[[3]]*dd1$sum 
	
plot(dd1$sum, par.resid.sum, cex.lab=1.2, ylab='AUC(partial.resid.sum)', main='auc~obs + sum, ~1 | class/order/family', pch=20)
lines(sum.seq, summary(lme1)$tTable[[3]]*sum.seq)
	
dev.off()
	
# obs~pred
dd = subset(dd, obs!=0)
str(dd)
	
lme2 = lme(obs~pred + sum, random= ~1 | class/order/family, method='ML', data = dd, control = lmeControl(opt='optim', msMaxIter = 200, msMaxEval = 500, msVerbose = T))
	
summary(lme2)
	
# pdf(here(outputpath,'prediction_outputs','sjsdm-graph', sjsdmVfolder, 'prediction', glue('lme-obs_{set}_{datafolder}_{abund}_{formula}.pdf')), height=12,width=17)
	
par(mfrow=c(2,3),cex=1.05)
	
hist(dd$pred, cex.lab=1.2)
hist(dd$sum, cex.lab=1.2)
hist(dd$obs, cex.lab=1.2)
	
pred.seq = seq(min(dd$pred),max(dd$pred), length=100)
sum.seq = seq(min(dd$sum),max(dd$sum), length=100)
fit = summary(lme2)$tTable[[1]] + summary(lme2)$tTable[[2]]*pred.seq + summary(lme2)$tTable[[3]]*sum.seq 
	
plot(dd$pred, dd$obs, main='lme', cex.lab=1.2, pch=20)
lines(pred.seq, fit)
	
par.resid.pred = resid(lme2) + summary(lme2)$tTable[[2]]*dd$pred
	
plot(dd$pred, par.resid.pred, main='lme', cex.lab=1.2, ylab='AUC(partial.resid.pred)', pch=20)
lines(pred.seq, summary(lme2)$tTable[[2]]*pred.seq)
	
par.resid.sum = resid(lme2) + summary(lme2)$tTable[[3]]*dd$sum
	
plot(dd$sum, par.resid.sum, cex.lab=1.2, ylab='AUC(partial.resid.sum)', main='obs~pred + sum, ~1 | class/order/family', pch=20)
lines(sum.seq, summary(lme2)$tTable[[3]]*sum.seq)
	
dev.off()
	
# ..... obs(0/1)~pred .....
#predict.lme.x1 <- logit.lme[1] + predicted.df.x1$x1*logit.lme[2] + predicted.df.x1$x2*logit.lme[3]
#	predict.lme.x1 <- exp(predict.lme.x1)/(1+exp(predict.lme.x1))
#predicted.df.x2 <- data.frame(x2 = seq(min(df$x2), max(df$x2), length.out = 1000), x1 = mean(df$x1))
#predicted.df.x2 <- as.data.frame(do.call(rbind, replicate(3, predicted.df.x2, simplify = F)))
#predicted.df.x2$x3 <- rep(c('a', 'b', 'c'), each = 1000)
#plot(df$x1, predict(test.glmer, type = 'response'), pch = 16, cex = 0.4, col = 'blue', xlab = 'x1', ylab = 'Predicted Prob.', main = 'glmer Function')
#lines(predicted.df.x1$x1[1:1000], predict.glmer.x1[1:1000], col = 'black', lwd = 2)

#pred.seq = seq(min(dd$pred),max(dd$pred), length.out=100)
#sum.seq = seq(min(dd$sum),max(dd$sum), length.out=100)
#fit <- predict(glmer1, newdata=data.frame(pred=pred.seq, sum=sum.seq), type = 'response')
	
glmer1 = glmer(obs ~ pred + sum + (1|class/order/family), data = dd, family = binomial,  control = glmerControl(optimizer = "bobyqa"), verbose = 1)
isSingular(glmer1)
	
summary(glmer1)$coefficients
	
# pdf(here(outputpath,'prediction_outputs','sjsdm-graph', sjsdmVfolder, 'prediction', glue('lme-obs_{set}_{datafolder}_{abund}_{formula}.pdf')), height=12,width=17)
	
par(mfrow=c(2,3),cex=1.05)
	
hist(dd$pred, cex.lab=1.2)
hist(dd$sum, cex.lab=1.2)
hist(dd$obs, cex.lab=1.2)
	
dd1 = dd[seq(1,dim(dd)[1],by=2),]
	
plot(dd1$pred, dd1$obs, pch=20, main='lme')
points(dd1$pred, predict(glmer1, dd1, type = 'response'), pch=4, col='lightblue')
	
plot(dd1$sum, dd1$obs, pch=20, main='obs~pred + sum, ~1 | class/order/family')
points(dd1$sum,predict(glmer1, dd1,type = 'response'), pch=4, col='lightblue')
	
dev.off()
	

```



