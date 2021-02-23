# Feb 04, 2021
# last modified: Feb 21, 2021
# play with tuned model (manually tuning in ADA cluster, 20210215-0221)



```{r setup}
# setwd('/media/yuanheng/SD-64g2/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/13_sjsdm_prediction')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan','pROC', 'gridExtra','ggeffects','cowplot')
	
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
	
trap <- "M1"; period = "S1"
date.model.run = 20210221   # according to tuning result
abund = 'pa'		# 'qp','pa'  !!! change accordingly
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
datapath=glue('../12_sjsdm_general_model_outputs')
	
# ...... from rdata ........
if (formula.env == "envDNN") {load(here(datapath,'source','yuanheng_mod_data.rdata'))}
if (formula.env == "newVars") {load(here(datapath,'source',glue('forada_data_{period}_m1m2_min{minocc}_{date.model.run}_{formula.env}.rdata')))}
	
#save(s.otu.test, s.otu.train, scale.env.test, scale.env.train, XY.test, XY.train, file=here(datapath,'source','yuanheng_mod_data_newVars.rdata'))
	
# 0/1
if (abund == 'pa')
{
	otu.train = as.data.frame((otu.train>0)*1)
	otu.test = as.data.frame((otu.test>0)*1)
	str(otu.test)
}
	
str(otu.train); abund
str(scale.env.test); formula.env
	
s.otu.train = as.matrix(otu.train)	
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
s.otu.test = as.matrix(otu.test)		
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	

```


```{r analyze-tune}
nstep=1000
plot.date=20210222
formula.env; abund; minocc 		
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
	
# analyze tuning data from ada (3rd week of 2021Feb)
tuning.ddI = read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune',glue('{formula.env}_{date.model.run}'), glue('manual_tuning_sjsdm_5CV_{trap}{period}_{abund}_min_{minocc}_nSteps_{nstep}.csv')), header=T, sep=',')
str(tuning.ddI)
tuning.ddS = read.table(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune',glue('{formula.env}_{date.model.run}'), glue('manual_tuning_sjsdm_5CV_{trap}{period}_mean_AUC_{abund}_min_{minocc}_nSteps_{nstep}.csv')), header=T, sep=',')
str(tuning.ddS)
	
tuning.dd = tuning.ddS
tuning.dd = dplyr::rename(tuning.dd, 'diff.train_mean'='tjur.train_mean', 'diff.test_mean'='tjur.test_mean')
	
which.max(tuning.dd$AUC.test_mean); which.max(tuning.dd$cor.test_mean)
which.max(tuning.dd$nagel.test_mean); which.max(tuning.dd$ll.test_mean)
	
maxdd = select(tuning.dd, ends_with('.test_mean'))
pp = c(sapply(1:(ncol(maxdd)-1), function(i) which.max(maxdd[,i])),which.min(maxdd$auc.lt5.test_mean))
	
# !!! if auc==cor; nagel==ll !!!
maxdd = select(tuning.dd, ends_with('.test_mean'), -starts_with('ll'), -'cor.test_mean')
pp = c(sapply(1:(ncol(maxdd)-1), function(i) which.max(maxdd[,i])),which.min(maxdd$auc.lt5.test_mean))
	
# pdf(here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('{formula.env}_{date.model.run}'),glue('plot_tuning_sjsdm_{trap}{period}_{abund}_min{minocc}_{formula.env}_{date.model.run}_{plot.date}.pdf')), width=6, height=4.5)
	
par(cex=.7)
names(tuning.dd)
ycol = 18; xcol = 12; xann = 'log-likelihood'; yann = 'AUC'; yrange = c(0.5,1)
ycol.train = 11
	
plot(y=tuning.dd[,ycol], x=tuning.dd[,xcol], pch=3, col='blue',ylim=yrange,xlab=xann,ylab=yann)
points(y=tuning.dd[,ycol.train], x=tuning.dd[,xcol], pch=8)
points(y=tuning.dd[pp,ycol], x=tuning.dd[pp,xcol], pch=19, col='red')
text(x=tuning.dd[pp,xcol], y=seq(.5,.55,length.out=length(pp)), strsplit(names(maxdd),'.test_mean'))
abline(v=tuning.dd[pp,xcol], lty=2, col='red')
	
legend('topright',pch=c(3,8),col=c('blue','black'),c('auc.test','auc.train'),bty='n')
#legend('topright',pch=c(3,3,8,8),col=c('blue','blue','black','black'),c(paste('auc.test,',round(tuning.dd[pp,]$AUC.test,3)),paste('auc.train,',round(tuning.dd[pp,]$AUC.explain,3))),bty='n')
	
dev.off()
	
tuning.dd[pp,]
	
# model with optimal parameters
lambda.env = tuning.dd[pp,'lambda.env']		# 0.0 0.2, (pa) 0.2
alpha.env = tuning.dd[pp,'alpha.env']		# 0.7 1.0, (pa) 0.9
lambda.sp = tuning.dd[pp,'lambda.sp']		# 0.8333333 0.1666667, (pa) 0.6666667
alpha.sp =  tuning.dd[pp,'alpha.sp']		# 1.0000000 0.6666667, (pa) 0.5
lambda.bio = tuning.dd[pp,'lambda.bio']		# 0.9 0.3, (pa) 1
alpha.bio =  tuning.dd[pp,'alpha.bio']		# 1 0, (pa) 0.5
hidden1 = tuning.dd[pp,'hidden.ind']			# 2 2, (pa) 2
acti = as.character(tuning.dd[pp,'acti.sp'])
drop = tuning.dd[pp,'drop']					# 0.1 0.1, (pa) 0.1
	
str(s.otu.train)
str(scale.env.test)
	
for (i in 1:length(pp)) {
model.train = sjSDM(Y = s.otu.train,
	  env = DNN(data=scale.env.train, formula = ~.,
	  hidden=hidden[[hidden1[i]]], lambda = lambda.env[i], alpha = alpha.env[i], activation=acti[i], dropout=drop[i], bias=T),
	  
	  biotic = bioticStruct(lambda=lambda.bio[i], alpha=alpha.bio[i], on_diag=F, inverse = FALSE),
	  
	  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp[i], alpha=alpha.sp[i]),
	  
	  learning_rate = 0.003, # 0.003 recommended for high species number 
	  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L # 150L, 5000L
)
	 
# saveRDS(model.train, here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS', glue('{formula.env}_{date.model.run}'), glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{lambda.bio[i]}_{drop[i]}_{date.model.run}.RDS')) )
}
	
model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder, 'sjsdm-model-RDS',glue('{formula.env}_{date.model.run}'), glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{lambda.bio[i]}_{drop[i]}_{date.model.run}.RDS')) )
	
# pdf(here(outputpath, 'sjsdm_general_outputs',sjsdmVfolder, 'DNN_tune', glue('{formula.env}_{date.model.run}'),glue('model-history_tuned.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{names(maxdd)[i]}_{date.model.run}.pdf')), width=5, height=5)
	
par(cex=.8)
plot(model.train$history)
mtext(paste0(names(maxdd)[i],': ', round(select(tuning.dd[pp,],names(maxdd)[i])[i,],3)), side=3,line=-1, adj=.95)
mtext(paste0(names(tuning.dd)[c(2:7,9:10)],': ',select(tuning.dd[pp,],names(tuning.dd)[c(2:7,9:10)])[i,]), side=3,line=seq(-2,-9), adj=.95)
	
dev.off()
	

```


```{r prediction}
model1 = model.train
formula = formula.env 
	
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
	
#table(names(otu.test)==names(otu.train))
otudd = as.data.frame(otudd)
names(otudd)=names(otu.train)
str(otudd)
	
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
#roc.dd2 = sapply(1:ncol(otudd), function(i) Metrics::auc(otudd.pa[,i], pred.dd[,i]))
	
#table(round(roc.dd,6) == round(roc.dd2,6))
	
auc.mean = mean(roc.dd)
formula; abund; names(maxdd)[i]; auc.mean
tuning.dd[pp[i],]
	
# saveRDS(list(pred.Y=pred.dd, otu=otudd, roc.allS=roc.dd, auc.mean=auc.mean), here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, glue('{formula.env}_{date.model.run}'), glue('roc_result_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_{names(maxdd)[i]}_{date.model.run}.RDS')))
	
```


```{r roc_result-for-regression}
# .. set names
# predictive
set = 'test' # 'explain', test
formula 
	
# ....... needed for selecting data ..........
explain = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder,glue('{formula}_{date.model.run}'), glue('roc_result_explain_S1_M1_pa_min5_newVars_AUC.test_mean_20210221.RDS'))) 
names(explain)
	
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
roc.dd = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, glue('roc_result_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env[i]}_{alpha.env[i]}_{lambda.sp[i]}_{drop[i]}_{date.model.run}.RDS'))) 
names(roc.dd)
	
names(roc.dd$otu)
	
# ... make long table 
roc.dd1 = as.data.frame(bind_cols(roc.dd$otu %>% pivot_longer(everything(), names_to='otu.name', values_to = 'obs'), pred=as.vector(t(roc.dd$pred.Y)), auc=rep(roc.dd$roc.allS, nrow(roc.dd$otu))))
roc.dd1 = left_join(roc.dd1, a, by = c('otu.name'='otu'), copy=T)
str(roc.dd1)
# 18*230 = 4140
# 70*268 = 18760
	
# ......... check on dataset .........
table(roc.dd1$sum.seq[1:243])
roc.dd1$otu.name = as.factor(roc.dd1$otu.name)
	
roc.dd$otu[1,1:6]; roc.dd$pred.Y[1,1:6]; roc.dd1[1:6,2:3]
unique(roc.dd1$otu.name[roc.dd1$auc==0]) ; unique(roc.dd1$otu.name[roc.dd1$auc==1])
# envDNN,qp -> unique(roc.dd1$otu.name[roc.dd1$auc==0])[1:2]; unique(roc.dd1$otu.name[roc.dd1$auc==1])[c(2,6)]
```


```{r violin-auc}
plot.date=20210210
dd = roc.dd1
names(dd)
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction', glue('box_order_{set}_{period}_{trap}_{abund}_min{minocc}_{formula}_{lambda.env[i]}_{date.model.run}_{plot.date}.pdf')), width=12, height=6)
	
cc = a %>% count(order)
cc = cc[order(cc$n, decreasing=T),]
	
dd1 = dd[sapply(1:length(unique(dd$otu.name)), function(x) {which(dd$otu.name==unique(dd$otu.name)[x])[1]} ),]
str(dd1)
	
ttt=1
if (set=='explain') {ttt=round(tuning.dd$AUC.explain[pp[i]],3)}
if (set=='test') {ttt=round(tuning.dd$AUC.test[pp[i]],3)}
	
ggplot(dd1, aes(x=order, y=auc)) + geom_hline(yintercept=.75, col='gray') + geom_jitter(width = 0.35, aes(colour = order, size=sum)) + geom_boxplot(width=0.1) + scale_x_discrete(limit=cc$order) + annotate(geom="text", y=-.43, x=1:length(cc$n), label=as.character(cc$n), col='red') + theme_minimal() + theme(legend.position = c(.95,0.3) ) + guides(color = FALSE) + annotate(geom="text", x=2, y=-.3, label=glue('AUC.{set}: {ttt}, {round(roc.dd$auc.mean,3)}'))
	
dev.off()
	
```


```{r correlate-auc}
# .... (define) ....
minocc
auc.all = data.frame(otu=as.character(names(otu.train)), auc.test=rep(.1,length=ncol(otu.train)), auc.exp=rep(.1,length=ncol(otu.train)))
str(auc.all)
# dim[1] == 352, min5; 303, min6;
	
set = c('test','explain')  # 1'test' , 2'explain'
formula1 = formula.env 
	
names(maxdd)
	
for (i in 1:ncol(maxdd)) {
	
	roc.dd.t = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder,glue('{formula1}_{date.model.run}'), glue('roc_result_{set[1]}_{period}_{trap}_{abund}_min{minocc}_{formula1}_{names(maxdd)[i]}_{date.model.run}.RDS'))) 
	names(roc.dd.t)
	roc.dd.e = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder,glue('{formula1}_{date.model.run}'), glue('roc_result_{set[2]}_{period}_{trap}_{abund}_min{minocc}_{formula1}_{names(maxdd)[i]}_{date.model.run}.RDS'))) 
	names(roc.dd.e)
	formula = paste0('newVars-',strsplit(names(maxdd)[i],'.test')[[1]][1]) 
	# ... make long table
	b.t = data.frame( auc.test = roc.dd.t$roc.allS, otu.t = names(roc.dd.t$otu) )
	str(b.t)
	
	b.e = data.frame( auc.exp = roc.dd.e$roc.allS, otu.e = names(roc.dd.e$otu) )
	str(b.e)
	
	auc.te = inner_join(b.t, b.e, by=c('otu.t'='otu.e'))
	str(auc.te)
	
	auc.all = left_join(auc.all, auc.te, by=c('otu'='otu.t'), suffix=c('', glue('.{formula}')), copy=T)
	str(auc.all)
	
}
	 
# .. after all variables are added
auc.all = dplyr::select(auc.all, -'auc.test',-'auc.exp')
	
# ... extract taxonomy info
auc.all = left_join(auc.all, select(a, 'otu','order','class','family','sum','n'), by=c('otu'='otu'))
abc = data.frame(seq.order=letters[1:length(unique(a$n))], order=sort(unique(a$n),decreasing=T))
auc.all$Oorder = sapply(1:dim(auc.all)[1], function(x) paste(abc$seq.order[abc$order==auc.all$n[x]],'.',auc.all$order[x],'.',auc.all$n[x], sep=''))
str(auc.all)
	
# ....... regression .......
formula = formula1 
	
# ... loop ...
setS = sapply(1:ncol(maxdd), function(i) strsplit(names(maxdd)[i],'.test')[[1]][1])
	
dd.ggplot = vector(mode = "list", length = ncol(maxdd))
plot.list=list(ggplot(),ggplot(),ggplot(),ggplot(),ggplot(),ggplot(),ggplot())
#plot.list=list(ggplot(),ggplot(),ggplot(),ggplot(),ggplot())
	
# test vs. train
for (j in 1:ncol(maxdd)) {
	set = setS[j]
	ii=j*2; jj=ii+1
	print(c(names(auc.all)[ii],names(auc.all)[jj]))
	
	ab = strsplit(names(auc.all)[ii],'[.]')[[1]][2]
	auc.1 = select(auc.all, ii, jj, 'sum', 'order', 'class', 'family', 'Oorder') %>% rename(auc.test=1, auc.train=2, incidence=3)
	auc.1 = na.omit(auc.1)
	dd.ggplot[[j]] = auc.1
	
	gp = ggplot(dd.ggplot[[j]], aes(auc.train, auc.test)) + geom_point(aes(colour=factor(Oorder), size=incidence))+ scale_size(range = c(1, 7)) + scale_colour_manual(values=colorRampPalette(c('dodgerblue3','firebrick2','yellow'))(length(unique(auc.1$Oorder)))) + geom_smooth(method='lm', se=F, colour='gray') + geom_abline(slope=1,intercept=0, linetype = "dashed", colour='gray', size=1.5) + theme(panel.background=element_rect(fill='snow')) + ggtitle(glue('{set}, {formula}, min{minocc}')) + geom_hline(yintercept=0.5, linetype = "dashed", colour='red') + geom_vline(xintercept=0.5, linetype = "dashed", colour='red') + xlim(0,1) +ylim(0,1)+ annotate(geom="text", y=.05, x=.15, label=glue('auc.train: {round(mean(auc.1$auc.train),3)}, auc.test: {round(mean(auc.1$auc.test),3)}'))
	if (j==1 |j==3 | j==5 | j==7) {plot.list[[j]] = gp + theme(legend.position='none')} else {plot.list[[j]] = gp + theme(legend.position='left')}
	
	
}
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction', glue('{formula1}_{date.model.run}'),glue('auc-test-train_{formula1}_{period}_{trap}_min{minocc}_tuned_{date.model.run}.pdf')), width=14, height=20)
	
#grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]], nrow=3, widths=c(.44,.56))  
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],plot.list[[7]],widths=c(.44,.56), nrow=4)#layout_matrix= rbind(c(1,2), c(3,4), c(5,6,7)),
	
dev.off()
	

# ............. lme model ................
auc.all$order =  as.factor(auc.all$order)
names(auc.all)[1:5]
i = 2
abund='qp'; formula='envDNN'
# auc.qp.test.envDNN.newVars
	
dd = auc.all[which(auc.all[,i]!='NA'),]
mod1 <- lme4::lmer(dd[,i] ~ sum + (1 | order), data = dd, REML = FALSE)
mod2 <- lme4::lmer(dd[,i] ~ 1 + (1 | order), data = dd, REML = FALSE)
anova(mod1, mod2)
	
mod1 <- lme4::lmer(dd[,i] ~ sum + (1 | order), data = dd, REML = T)
r2=MuMIn::r.squaredGLMM(mod1)
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction', glue('lme_{abund}_{period}_{trap}_min{minocc}_{formula}_tuned_{date.model.run}.pdf')), width=6, height=8.5)
	
#plot(dd$sum, dd[,i], xlab='incidence',ylab=paste0('auc.',abund,'.test.',formula))
## auc.all$auc.qp.test.envDNN.newVars
#abline(a=coef(mod1)$order[,1], b=coef(mod1)$order[,2])
	
p1 = ggplot(dd, aes(x = sum, y = dd[,i], group = order, colour = order)) + geom_point() + geom_line(data = cbind(dd, pred = predict(mod1)), aes(y = pred), size = 1) + labs(x = "incidence", y = paste0('AUC.',abund,'.test.',formula)) + annotate(geom='text',x=(max(dd$sum)*.9), y=.1, label=paste0('R^2: ',round(r2[[1]],3))) + scale_colour_viridis_d(option = "cividis") + theme(legend.position='bottom')
p2 = ggplot(dd, aes(x = sum, y = dd[,i], group = order, colour = order)) + geom_point() + geom_line(data = cbind(dd, pred = predict(mod1)), aes(y = pred), size = 1) + labs(x = "incidence", y = paste0('AUC.',abund,'.test.',formula)) + theme_cowplot() + facet_wrap(vars(order)) + scale_colour_viridis_d(option = "cividis") + theme(legend.position='none')
grid.arrange(p1,p2, nrow=2, heights=c(.65,.35))
	
dev.off()
	
```


```{r I-map-plot}
# load data 
set = c('test','explain')  
abund; formula  
	
dd.test=data.frame(); dd.explain=data.frame()
for (i in 1:2) {
	if (formula=='envDNN' & abund=='qp') {lambda.env1=0.2; alpha.env1=1; lambda.sp1=0.166666666666667; drop1=0.1}		# envDNN, qp
	if (formula=='envDNN.newVars' & abund=='qp') {lambda.env1=0; alpha.env1=1; lambda.sp1=0.666666666666667; drop1=0.1}		# newVars, qp
	if (formula=='envDNN' & abund=='pa') {lambda.env1=0.2; alpha.env1=0.9; lambda.sp1=0.666666666666667; drop1=0.1}		# envDNN, pa
	if (formula=='envDNN.newVars' & abund=='pa') {lambda.env1=0.1; alpha.env1=1; lambda.sp1=0.5; drop1=0.1}						# newVars, pa
	
	dd = readRDS(here(outputpath,'prediction_outputs','sjsdm-model-RDS', sjsdmVfolder, glue('roc_result_{set[i]}_{period}_{trap}_{abund}_min{minocc}_{formula}_lambdaE{lambda.env1}_{alpha.env1}_{lambda.sp1}_{drop1}_{date.model.run}.RDS'))) 
	dd$pred.Y = as.data.frame(dd$pred.Y); names(dd$pred.Y)=names(dd$otu)
	dd$roc.allS = data.frame(auc=dd$roc.allS, name=names(dd$otu))
	if (set[i]=='test') {dd.test = dd}; if (set[i]=='explain') {dd.explain = dd}
}
	
# .. select OTU 
names(dd.test); str(dd.test$pred.Y)
	
dd.map = bind_rows( bind_cols(select(dd.test$otu, dd.test$roc.allS$name[dd.test$roc.allS$auc<.06][1:2], dd.test$roc.allS$name[dd.test$roc.allS$auc>.99][1:2]), XY.test, select(dd.test$pred.Y, dd.test$roc.allS$name[dd.test$roc.allS$auc<.06][1:2], dd.test$roc.allS$name[dd.test$roc.allS$auc>.99][1:2]) %>% rename_with(~ tolower(gsub("R", "pred.", .x, fixed = TRUE))), type='test'), bind_cols(select(dd.explain$otu, dd.test$roc.allS$name[dd.test$roc.allS$auc<.06][1:2], dd.test$roc.allS$name[dd.test$roc.allS$auc>.99][1:2]), XY.train, select(dd.explain$pred.Y, dd.test$roc.allS$name[dd.test$roc.allS$auc<.06][1:2], dd.test$roc.allS$name[dd.test$roc.allS$auc>.99][1:2]) %>% rename_with(~ tolower(gsub("R", "pred.", .x, fixed = TRUE))), type='train'))
str(dd.map)
	
dd.map[(nrow(dd.map)+1),1:(ncol(dd.map)-1)] = c(0,0,0,0,-2,-2,0,0,0,0); dd.map[nrow(dd.map),ncol(dd.map)] = 'z'
dd.map[(nrow(dd.map)+1),1:(ncol(dd.map)-1)] = c(1,1,1,1,-2,-2,1,1,1,1); dd.map[nrow(dd.map),ncol(dd.map)] = 'z'
dd.map$type = factor(dd.map$type,levels=c('train','test','z'))
	

```

```{r II-map-plot}
plotlist=list(ggplot(),ggplot(),ggplot(),ggplot())
	
# envDNN, qp
for(i in 1:4) {
	print(i)
	auctest=round(dd.test$roc.allS$auc[dd.test$roc.allS$name==names(dd.map)[i]],3); auctrain=round(dd.explain$roc.allS$auc[dd.explain$roc.allS$name==names(dd.map)[i]],3)
#	names(dd.map)[i]; names(dd.map)[i+6]; 
	if (i==1) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R9119__Insecta_Diptera_Muscidae_Phaonia_nigricauda_BOLD_AAP6480_size.1261, colour=pred.9119__insecta_diptera_muscidae_phaonia_nigricauda_bold_aap6480_size.1261))}
	if (i==2) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R2343.18__Insecta_Lepidoptera_Noctuidae_Diarsia_esurialis_BOLD_ABX6710_size.304, colour=pred.2343.18__insecta_lepidoptera_noctuidae_diarsia_esurialis_bold_abx6710_size.304))}
	if (i==3) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R2197.71__Insecta_Diptera_Sciaridae_Bradysia_placida_BOLD_ACR4350_size.42, colour=pred.2197.71__insecta_diptera_sciaridae_bradysia_placida_bold_acpred.4350_size.42))}
	if (i==4) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=type, size=R5811.2__Insecta_Raphidioptera_Raphidiidae_Agulla_unicolor_BOLD_ACA6995_size.789, colour=pred.5811.2__insecta_pred.aphidioptera_pred.aphidiidae_agulla_unicolor_bold_aca6995_size.789))}
	
	plotlist[[i]] = g1st + geom_point() + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9), panel.background=element_rect(fill='gray98')) + labs(shape='data type', size='obs', colour='pred') + ggtitle(paste0(strsplit(strsplit(names(dd.map)[i],'__')[[1]][2],'D_')[[1]][1],', auc.test: ',auctest,', auc.train: ',auctrain)) + scale_colour_gradient2(low='black', mid='red',high='lavender',midpoint=.5)
	
} 
	
# envDNN, pa
for(i in 1:4) {
	print(i)
	auctest=round(dd.test$roc.allS$auc[dd.test$roc.allS$name==names(dd.map)[i]],3); auctrain=round(dd.explain$roc.allS$auc[dd.explain$roc.allS$name==names(dd.map)[i]],3)
#	names(dd.map)[i]; names(dd.map)[i+6]; 
	if (i==1) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R2339.2__Insecta_Diptera_Tabanidae_Chrysops_NA_BOLD_AAP7670_size.2174, colour=pred.2339.2__insecta_diptera_tabanidae_chrysops_na_bold_aap7670_size.2174))}
	if (i==2) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R5499__Insecta_Hymenoptera_Ichneumonidae_NA_NA_BOLD_AAU8873_size.147, colour=pred.5499__insecta_hymenoptera_ichneumonidae_na_na_bold_aau8873_size.147))}
	if (i==3) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R4091.6__Insecta_Diptera_Asilidae_Nevadasilus_NA_BOLD_AAH2298_size.3235, colour=pred.4091.6__insecta_diptera_asilidae_nevadasilus_na_bold_aah2298_size.3235))}
	if (i==4) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=type, size=R10098__Insecta_Diptera_Mycetophilidae_Mycomya_NA_BOLD_ACJ6978_size.215, colour=pred.10098__insecta_diptera_mycetophilidae_mycomya_na_bold_acj6978_size.215))}
	
	plotlist[[i]] = g1st + geom_point() + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9), panel.background=element_rect(fill='gray98')) + labs(shape='data type', size='obs', colour='pred') + ggtitle(paste0(strsplit(strsplit(names(dd.map)[i],'__')[[1]][2],'D_')[[1]][1],', auc.test: ',auctest,', auc.train: ',auctrain)) + scale_colour_gradient2(low='black', mid='red',high='lavender',midpoint=.5)
	
} 
	
# envDNN.newVars, qp
for(i in 1:4) {
	print(i)
	auctest=round(dd.test$roc.allS$auc[dd.test$roc.allS$name==names(dd.map)[i]],3); auctrain=round(dd.explain$roc.allS$auc[dd.explain$roc.allS$name==names(dd.map)[i]],3)
#	names(dd.map)[i]; names(dd.map)[i+6] 
	if (i==1) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R10103.4__Insecta_Diptera_Mycetophilidae_Syntemna_NA_BOLD_ACG6555_size.20, colour=pred.10103.4__insecta_diptera_mycetophilidae_syntemna_na_bold_acg6555_size.20))}
	if (i==2) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R3045.20__Insecta_Diptera_Piophilidae_Actenoptera_NA_BOLD_AAG0471_size.87, colour=pred.3045.20__insecta_diptera_piophilidae_actenoptera_na_bold_aag0471_size.87))}
	if (i==3) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R10178.4__Insecta_Diptera_Cecidomyiidae_NA_NA_BOLD_ACC2335_size.24, colour=pred.10178.4__insecta_diptera_cecidomyiidae_na_na_bold_acc2335_size.24))}
	if (i==4) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=type, size=R1519.2__Insecta_Diptera_NA_NA_NA_BOLD_ACX8657_size.2595, colour=pred.1519.2__insecta_diptera_na_na_na_bold_acx8657_size.2595))}
	
	plotlist[[i]] = g1st + geom_point() + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9), panel.background=element_rect(fill='gray98')) + labs(shape='data type', size='obs', colour='pred') + ggtitle(paste0(strsplit(strsplit(names(dd.map)[i],'__')[[1]][2],'D_')[[1]][1],', auc.test: ',auctest,', auc.train: ',auctrain)) + scale_colour_gradient2(low='black', mid='red',high='lavender',midpoint=.5)
	
} 
	
# envDNN.newVars, pa
for(i in 1:4) {
	print(i)
	auctest=round(dd.test$roc.allS$auc[dd.test$roc.allS$name==names(dd.map)[i]],3); auctrain=round(dd.explain$roc.allS$auc[dd.explain$roc.allS$name==names(dd.map)[i]],3)
#	names(dd.map)[i]; names(dd.map)[i+6]
	if (i==1) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R2339.2__Insecta_Diptera_Tabanidae_Chrysops_NA_BOLD_AAP7670_size.2174, colour=pred.2339.2__insecta_diptera_tabanidae_chrysops_na_bold_aap7670_size.2174))}
	if (i==2) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R5928.2__Insecta_Lepidoptera_Copromorphidae_Ellabella_editha_BOLD_ABX3483_size.558, colour=pred.5928.2__insecta_lepidoptera_copromorphidae_ellabella_editha_bold_abx3483_size.558))}
	if (i==4) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=as.factor(type), size=R1519.2__Insecta_Diptera_NA_NA_NA_BOLD_ACX8657_size.2595, colour=pred.1519.2__insecta_diptera_na_na_na_bold_acx8657_size.2595))}
	if (i==3) {g1st = ggplot(dd.map, aes(UTM_E, UTM_N, shape=type, size=R10098__Insecta_Diptera_Mycetophilidae_Mycomya_NA_BOLD_ACJ6978_size.215, colour=pred.10098__insecta_diptera_mycetophilidae_mycomya_na_bold_acj6978_size.215))}
	
	plotlist[[i]] = g1st + geom_point() + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9), panel.background=element_rect(fill='gray98')) + labs(shape='data type', size='obs', colour='pred') + ggtitle(paste0(strsplit(strsplit(names(dd.map)[i],'__')[[1]][2],'D_')[[1]][1],', auc.test: ',auctest,', auc.train: ',auctrain)) + scale_colour_gradient2(low='black', mid='red',high='lavender',midpoint=.5)
	
} 
	
# pdf(here(outputpath,'prediction_outputs', 'sjsdm-graph', sjsdmVfolder, 'prediction','descriptive', paste0('auc-map_',abund,'_',period,'_',trap,'_min',minocc,'_',formula,'_tuned_',date.model.run,'_minmax.pdf')), width=17, height=13)
	
grid.arrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], nrow=2, layout_matrix= rbind(c(1,2), c(3,4)), heights=c(.5,.5))   # g.1, 
	
dev.off()
	

```





