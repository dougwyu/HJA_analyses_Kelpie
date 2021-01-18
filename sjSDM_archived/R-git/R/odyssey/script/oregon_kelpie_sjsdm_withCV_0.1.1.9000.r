# created: Oct 29, 2020
# last modified:

```{r setup}
# setwd('/media/yuanheng/SD-64g2/Downloads/backup2/HJA_analyses_Kelpie/sjSDM/R-git')
	
lapply(c("ggplot2",'vegan', 'tidyverse','gridBase', 'grid','gridExtra', 'ggcorrplot','here', 'conflicted','reticulate','sjSDM','mgsub','vcd','RColorBrewer', 'reshape2','glue'), library, character.only=T)
# 'scatterplot3d', 'labdsv',
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("importance", "sjSDM")
	
dr_here()
packageVersion('sjSDM')
#[1] ‘0.1.1.9000’
	
source(here("R", "source", "sjsdm_function.r"))
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
```

```{r set-names}
# load data from 'crossvalidation' folder 
# they're for M1S1 ???
date.cross.validation = 20201030
minocc = 5
envvar = 'gismslidar'
abund = 'qp'
session = 's1'; malaise = 'm1'
	
datafolder = glue("data_{date.cross.validation}_{minocc}minocc_{envvar}")
	
cvfolder = glue("results_{date.cross.validation}_{minocc}minocc_{envvar}_{abund}_loocv")
	
otudata = glue('otu.{abund}.csv')
	
sjsdmV = '0.1.1.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
describefolder1 = glue('sjsdm_{session}_{malaise}_{abund}_min{minocc}')
describefolder2 = glue('model_{envvar}_{date.cross.validation}')
	
```

```{r load files}
# .. otu data
otu = read.csv(here('results', 'crossvalidation', cvfolder, datafolder, otudata), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(otu[,40:46])
	
dim(otu)
# [1]  88 268 
	
# .. covariates
scale.env = read.csv(here('results', 'crossvalidation', cvfolder, datafolder, 'scale.env.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(scale.env)
	
# .. spatial data
XY = read.csv(here('results', 'crossvalidation', cvfolder, datafolder, 'XY.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(XY)
	
dim(XY); dim(scale.env); dim(otu)
# 2+27+268 = 297
	
```

```{r basic-check}
par(mfrow=c(1,2))
hist(XY$UTM_E)
hist(XY$UTM_N)
# spatial data are scaled
	
par(mfrow=c(1,2))
hist(scale.env$elevation_m)
hist(scale.env$B1_20180717)
	
par(mfrow=c(1,2))
hist(scale.env$NDVI_20180717)
hist(scale.env$EVI_20180717)
# spans of NDVI & EVI are large!!!
	
```

```{r sjsdm-model}
# ... session 1, Malaise I, (>4) ...
# ... with all environmental data
# . sjSDM-version = ‘0.1.1.9000’
# . specified algorithm-related parameters ???

# load tune result (only best, adagpu)
best = readRDS(here('results', 'crossvalidation', cvfolder, glue('sjsdm_tune_results_HJA_{date.cross.validation}_bestonly.RDS')))
best
	
# run sjSDM
model = sjSDM(Y = as.matrix(otu),
			  env = linear(data=as.matrix(scale.env), 
			  lambda=best[['lambda_coef']], alpha= best[['alpha_coef']]),
			  biotic = bioticStruct(lambda=best[['lambda_cov']], alpha=best[['alpha_cov']], on_diag=F, inverse = FALSE),
			  spatial = linear(data=as.matrix(XY), ~0+UTM_E:UTM_N, lambda =best[['lambda_spatial']] , alpha = best[['alpha_spatial']]),
			  learning_rate = 0.003, # 0.003 recommended for high species number 
			  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L
			 )
	
loss=unlist(model$logLik)
history=model$history
p = getSe(model)
result = list( beta = coef(model), sigma = getCov(model),loss=loss,history=model$history, p=p, logLik=logLik(model))
	
# saveRDS(result, file = here('results','sjsdm-model-RDS', glue('sjsdm-{sjsdmV}'), glue('s-jSDM_result_s1_m1_{abund}_min{minocc}_{envvar}_{date.cross.validation}.RDS')) )
# saveRDS(model, here('results','sjsdm-model-RDS', glue('sjsdm-{sjsdmV}'), glue('s-jSDM_model_s1_m1_{abund}_min{minocc}_{envvar}_{date.cross.validation}.RDS')) )
	
```

```{r load-model}
result = readRDS(here('results','crossvalidation', cvfolder, glue('sjsdm_result_HJA_{date.cross.validation}.RDS')))
names(result)
	
model = readRDS(here('results', 'crossvalidation', cvfolder, glue('sjsdm_model_HJA_{date.cross.validation}.RDS')))
names(model)
	
model["cl"]
model["settings"]
	
```

```{r graph-history.anova}
# .... model history
# pdf(here('results','sjsdm-graph', sjsdmVfolder, describefolder1,  describefolder2, glue('model-history_{envvar}_{abund}_min{minocc}_{date.cross.validation}_{sjsdmV}.pdf')), width=5, height=5)
	
plot(result$history)
	
dev.off()
	
# ... anova plot
# model1 = readRDS(here('results','sjsdm-model-RDS', sjsdmVfolder, glue('s-jSDM_model_s1_m1_{abund}_min{minocc}_{envvar}_{date.cross.validation}.RDS')) ) # run on 'cpu'
# an = anova(model1, cv = F)
	
an = readRDS(here('results','sjsdm-model-RDS', sjsdmVfolder,glue('anova_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}.RDS')))
names(an)
	
# pdf(here('results','sjsdm-graph',sjsdmVfolder, describefolder1,  describefolder2, glue('anova_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}_{sjsdmV}.pdf')), width=5, height=5)
	
plot(an)
title(paste('analysis of variance, ',ncol(result$sigma), ' OTUs', ', full-R2: ', round(an$result[9,5], 3), sep=''))
	
dev.off()
	
```


```{r graph-ternary}
# ... ternary plots ...
# imp =importance(model)
imp = readRDS(here('results', 'crossvalidation', cvfolder, glue('sjsdm_imp_HJA_{date.cross.validation}.RDS')))
	
# pdf(here('results','sjsdm-graph', sjsdmVfolder, describefolder1,  describefolder2, glue('vari.par_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}_{sjsdmV}.pdf')), width=5, height=5)
	
plot(imp, cex=.8) #from s-jSDM package
title(paste('variation partition, ',ncol(result$sigma), ' OTUs', sep=''))
	
dev.off()
	
# ... coloured ternary plot differentiate order & family
impD = data.frame(OTU=imp$names, imp$res$total)
str(impD)
	
impD$order = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[2])
table(impD$order)
	
impD$family = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[3])
sort(table(impD$family))
	
a = as.data.frame(impD %>% count(order))
a = subset(a, n>9)
impD1 = inner_join(impD, a, by=c('order'='order'))
str(impD1)
	
impD1$order = mgsub(impD1$order, c('Diptera', 'Coleoptera', 'Lepidoptera', 'Hemiptera', 'Hymenoptera'), c('Fly', 'Beetle & Weevil', 'Butterfly & Moth', 'Bug', 'Wasp & Bee'))
impD1$family = mgsub(impD1$family, c('Geometridae', 'Cecidomyiidae', 'Ichneumonidae'), c('Moth', 'Gall midge', 'Ichneumon wasp'))
# 'Muscidae', 'Fly'
str(impD1)
	
a = as.data.frame(impD1 %>% count(family))
a = subset(a, n>9)
impD2 = inner_join(impD1, a, by=c('family'='family'))
impD2 = subset(impD2, family!='BOLD' & family != 'NA')
str(impD2)
	
table(impD2$order)
table(impD2$family)
	
# pdf(here('results','sjsdm-graph',sjsdmVfolder, describefolder1, describefolder2, glue('vari.order_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}_{sjsdmV}.pdf')), width=8, height=5)
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Order (occurrence>9)'
vp.plot(data=impD1, ind=5, x1=2, x2=4, textM=titleText, textL=legendText)
	
dev.off()
	
#summary(cut(impD2$spatial, breaks = seq(0,1,length.out = 12)))
#summary(cut(impD2$biotic, breaks = seq(0,1,length.out = 12)))
#summary(cut(impD2$env, breaks = seq(0,1,length.out = 12)))
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Family (occurrence>9)'
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
# pdf(here('results','sjsdm-graph', sjsdmVfolder, describefolder1, describefolder2, glue('vari.family_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}_{sjsdmV}.pdf')), width=8, height=5)
	
vp.plot(data=impD2, ind=6, x1=2, x2=4, textM=titleText, textL=legendText)
	
dev.off()
	
# plot ranges of env, biotic, spatial
a = data.frame(range=names(summary(cut(impD$biotic, breaks = seq(0,1,length.out = 12)))), env=summary(cut(impD$env, breaks = seq(0,1,length.out = 12))), spatial=summary(cut(impD$spatial, breaks = seq(0,1,length.out = 12))), biotic=summary(cut(impD$biotic, breaks = seq(0,1,length.out = 12))), row.names=NULL)
a = melt(a, id="range", measure.vars=names(a)[2:4])
	
# pdf(here('results','sjsdm-graph',sjsdmVfolder, describefolder1, describefolder2, glue('vari.barplot_sjSDM_{session}_{malaise}_{abund}_min{minocc}_{envvar}_{date.cross.validation}_{sjsdmV}.pdf')), width=11, height=3)
	
ggplot(a, aes(range, value, fill = variable)) + geom_bar(stat="identity", position = "dodge") + scale_fill_brewer(palette = "Set1") #+ ggtitle(legendText)
	 
dev.off()
	
```

```{r graph-correlation}
# 27 variables
str(scale.env)
	
# ... correlation plot ...
dim(result$sigma)# cov <- result$sigma 
# [1] 268 268
otu.table = otu
str(otu.table[,1:5])
	
co.env.spp <- cov2cor(result$sigma)
rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
range(co.env.spp)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
rm(cut.t)
	
# ... species correlation plot
# pdf(here('results','sjsdm-graph',sjsdmVfolder, describefolder1, describefolder2, paste('species_correlation_', 'coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']],'.pdf', sep='')), height=12, width=12)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',sjsdmV,', tune',date.cross.validation,', coef-alpha: ', best[['alpha_coef']],', lambda: ',best[['lambda_coef']], sep=''))
	
dev.off()
	
# ..... circular Drawing .....
# . species association .
otu.tbl = otu
number=10
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
# pdf(here('results','sjsdm-graph', sjsdmVfolder, describefolder1, describefolder2, paste('species_covariance_circular_', 'coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']],'.pdf', sep='')), height=8, width=8)
	
version.text = paste('sjSDM version: ', sjsdmV,', tune',date.cross.validation,', coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']], sep='')
otu.text = "sum of quasiP"
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl, result=result)
	
dev.off()
	
# ... graph of circle with environmental factors
beta = as.matrix(data.frame(result$beta$env, result$beta$spatial)[,2:29])
effects = apply(beta, 1, function(o) sum(abs(o)))
str(effects)
	
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
	
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
str(effect_comb)
# variables index & value which has biggest coefficient
	
version.text = paste('sjSDM: ', sjsdmV,', tune',date.cross.validation,', coef-alpha: ', best[['alpha_coef']], ', lambda: ', best[['lambda_coef']], sep='')
otu.text = "spp. quasiP"
	
# names(scale.env)
evnames =c('HJA', 'ele', 'canopy','min.T', 'preci', 'road','stream', 'disturb', 'B1', 'B2', 'B3','B4','B5','B6','B7','B10','B11', 'NDVI','EVI','bright','green','wet','2m.max','2-4m','4-16m','p25','rumple', 'space')
	
# . plot of circle with max envir factor
# pdf(here('results','sjsdm-graph', sjsdmVfolder, describefolder1, describefolder2,'max_environ_spp_cov.pdf'), height=8, width=8)
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle.env(version.text=version.text, evnames=evnames, otu.text=otu.text,result=result, effect_comb=effect_comb, otu.tbl=otu.tbl)
	
dev.off()
	
```

```{r barplot-subset}
# continue with previous chalk
# names(scale.env)
sppnum = 40 # define how many spp to plot
	
str(effect_comb)
	
str(otu.tbl)
a = data.frame(sum=colSums(otu.tbl), otu=names(otu.tbl))
a = a[order(a$sum, decreasing=T),]
#a1 = a[c(1:10, (nrow(a)-9):nrow(a)), ]  
a1 = a[c(1:sppnum), ]
a1$rank=c(1:sppnum) 
a1$otu = as.character(a1$otu)
str(a1)
	
effect_comb$otu=names(otu.tbl)
effect_ind = dplyr::inner_join(effect_comb, a1, by=c('otu'='otu'))
effect_ind$name = paste(sapply(strsplit(sapply(str_split(effect_ind$otu, '__'), function(a) a[2]), '_'), function(a) paste(a[2],'_',a[3],'_',a[4], sep='')),'_', sapply(str_split(effect_ind$otu, '__'), function(a) a[1]), sep='')
sort(table(effect_ind$name))
	
beta1 = data.frame(result$beta$env, result$beta$spatial)[,2:29]
names(beta1) = evnames
beta1$otu = names(otu.tbl)
effect_ind = inner_join(effect_ind, beta1, by=c('otu'='otu'))
rm(beta1)
	
effect_ind$max_effects = as.character(effect_ind$max_effects)
a=unique(effect_ind$max_effects)
for (i in 1:length(a)) {effect_ind$max_effects = replace(effect_ind$max_effects, which(effect_ind$max_effects==a[i]), evnames[as.numeric(a[i])])}
str(effect_ind)
	
effect_ind = as.data.frame(pivot_longer(effect_ind, evnames, names_to='factors',values_to='beta'))
effect_ind$factors = factor(effect_ind$factors)
effect_ind$max = NA
effect_ind$max[effect_ind$max_effects==effect_ind$factors] = 'max'
	
effect_ind = effect_ind[order(effect_ind$rank),]
str(effect_ind)
	
col = sample(viridis::viridis(length(unique(effect_ind$max_effects))))
	
# pdf(here('results','sjsdm-graph',sjsdmVfolder, describefolder1, describefolder2, glue('max_envir_max{sppnum}spp_{envvar}_{abund}_min{minocc}_beta.pdf')), height=24, width=15)
	
bar.coef(effect_ind, .8)
	
dev.off()
	
```

