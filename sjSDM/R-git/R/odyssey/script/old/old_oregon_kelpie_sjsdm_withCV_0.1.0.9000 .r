# created: Sep 10, 2020
# last modified:
# store old scripts not used currently


```{r sjsdm-analysis-correlation-lidar (0/1)}
# . sjSDM-version = ‘0.1.0.9000’
# 0/1, >4
# c('elevation.scale', 'l_Cover_2m_max.scale', 'l_Cover_2m_4m.scale', 'l_Cover_4m_16m.scale', 'l_p25.scale', 'l_p95.scale', 'l_rumple.scale')]
str(scale.env1)
	
# . load result
model = readRDS(here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_model_s1_m1_present_no4_lidar_0.1.0.RDS'))
result = readRDS(here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_result_s1_m1_present_no4_lidar_0.1.0.RDS') )
summary(result)
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', 'model_history.pdf'), height=5, width=5)
	
plot(result$history)
	
best = readRDS(here('results', 'crossvalidation', 'results_20200904_5minocc_lidar_pa_loocv', "sjsdm_tune_results_HJA_20200904_bestonly.RDS"))
best
	
# ... correlation plot ...
dim(result$sigma)# cov <- result$sigma 
# [1] 274 274
otu.table = as.data.frame(dplyr::select(all.data, contains('__')))
str(otu.table[,1:5])
	
co.env.spp <- cov2cor(result$sigma)
rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
range(co.env.spp)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# ... species correlation plot
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', paste('species_correlation_', 'coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']],'.pdf', sep='')), height=12, width=12)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', result$regulation[2],', lambda: ', result$regulation[1],sep=''))
	
#dev.off()
	

# ..... Polygon Drawing .....
# . species association .
otu.tbl = as.data.frame(dplyr::select(all.data, contains('__')))
number=10
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
#pdf(here( 'results','sjsdm-graph', 'sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', paste('species_covariance_circular_', 'coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']],'.pdf', sep='')), height=8, width=8)
	
version.text = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200723, ', 'coef-alpha', best[['alpha_coef']], '_lambda', best[['lambda_coef']], sep='')
otu.text = "sum of presence"
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl, result=result)
	
#dev.off()
	

# . graph of polygon with environmental factors
# c('elevation.scale', 'l_Cover_2m_max.scale', 'l_Cover_2m_4m.scale', 'l_Cover_4m_16m.scale', 'l_p25.scale', 'l_p95.scale', 'l_rumple.scale')
beta = as.matrix(data.frame(result$beta$env, result$beta$spatial)[,2:9])
effects = apply(beta, 1, function(o) sum(abs(o)))
str(effects)
	
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
	
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
str(effect_comb)
# variables index & value which has biggest coefficient
	
version.text = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200723, coef-alpha: ', best[['alpha_coef']], ', lambda: ', best[['lambda_coef']], result$regulation[1], sep='')
otu.text = "presence spp."
	
# c('elevation.scale', 'l_Cover_2m_max.scale', 'l_Cover_2m_4m.scale', 'l_Cover_4m_16m.scale', 'l_p25.scale', 'l_p95.scale', 'l_rumple.scale')
evnames =c('ele', '2m_max','2m_4m', '4m_16m', 'p25', 'p95','rumple','space')
	
# . plot of polygon with max envir factor
#pdf(here('results','sjsdm-graph', 'sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model','max_environ_spp_cov.pdf'), height=8, width=8)
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle.env(version.text=version.text, evnames=evnames, otu.text=otu.text,result=result, effect_comb=effect_comb, otu.tbl=otu.tbl)
	
#dev.off()
	
```

```{r sjsdm-analysis-barplot-lidar (0/1)}
# continue with previous chalk
str(effect_comb)
	
str(otu.tbl)
a = data.frame(sum=colSums(otu.tbl), otu=names(otu.tbl))
a = a[order(a$sum, decreasing=T),]
#a1 = a[c(1:10, (nrow(a)-9):nrow(a)), ]  
#a1$rank=c(1:10,-10:-1) 
a1 = a[c(1:40), ]
a1$rank=c(1:40) 
a1$otu = as.character(a1$otu)
str(a1)
	
effect_comb$otu=names(otu.tbl)
effect_ind = dplyr::inner_join(effect_comb, a1, by=c('otu'='otu'))
effect_ind$name = sapply(strsplit(sapply(str_split(effect_ind$otu, '__'), function(a) a[2]), '_'), function(a) paste(a[2],'_',a[3],'_',a[4],'_',a[5],'_',a[6], sep=''))
table(effect_ind$name)
	
# change names
effect_ind$order = sapply(strsplit(sapply(str_split(effect_ind$name, '_'), function(a) a[1]), '_'), function(a) a[1])
table(effect_ind$order)
	
effect_ind$family = sapply(strsplit(sapply(str_split(effect_ind$name, '_'), function(a) a[2]), '_'), function(a) a[1])
table(effect_ind$family)
effect_ind$family[effect_ind$family=='BOLD'] = 'NA'
	
effect_ind$genus = sapply(strsplit(sapply(str_split(effect_ind$name, '_'), function(a) a[3]), '_'), function(a) a[1])
table(effect_ind$genus)
effect_ind$genus[effect_ind$genus=='ACF6271' | effect_ind$genus=='ACI9871'|effect_ind$genus=='ADC1012'] ='NA'
	
effect_ind$sp = sapply(strsplit(sapply(str_split(effect_ind$name, '_'), function(a) a[5]), '_'), function(a) a[1])
table(effect_ind$sp)
effect_ind$sp[effect_ind$sp=='BOLD'] = 'NA'
	
effect_ind$name = paste(effect_ind$order,'_', effect_ind$family, '_',effect_ind$genus,'_', effect_ind$sp, sep='')
effect_ind$name
	
effect_ind = as.data.frame(select(effect_ind, -'order',-'family',-'genus',-'sp'))
	
# done change name
beta1 = data.frame(result$beta$env, result$beta$spatial)[,2:9]
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
	
#pdf(here('results','sjsdm-graph', 'sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model','max_envir_max40spp_beta.pdf'), height=8, width=15)
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
bar.coef(effect_ind=effect_ind, t_turn=.9)
	
dev.off()
	
```


```{r sjsdm-model-lidar (quasiP)}
# ... session 1, Malaise I, (>4) ...
# ... with lidar environmental data
# . sjSDM-version = ‘0.1.0.9000’

all.data = dataI.1.quasiP5
otu.data = as.data.frame(select(all.data, contains('__'))) 
	
env.data = read.csv(here('results', 'crossvalidation', 'results_20200904_5minocc_lidar_qp_loocv', 'data_20200904_5minocc_lidar', 'scale.env1.csv'), header=T, sep=',')
names(env.data)
	
scale.env1 = as.data.frame(select(all.data, names(env.data)))
str(scale.env1)
	
XY = as.data.frame(select(all.data, starts_with('UTM_')) %>% scale)
	
par(mfrow=c(1,2))
hist(XY$UTM_E)
hist(XY$UTM_N)
	
# load tune result with only best (adagpu)
best = readRDS(here('results', 'crossvalidation', 'results_20200904_5minocc_lidar_qp_loocv', "sjsdm_tune_results_HJA_20200904_bestonly.RDS"))
best
	
# run sjSDM
model = sjSDM( Y = as.matrix(otu.data),
		   env = linear(data=as.matrix(scale.env1), 
#		   formula = ~ elevation.scale + l_Cover_2m_max.scale + l_Cover_2m_4m.scale + l_Cover_4m_16m.scale + l_p25.scale + l_p95.scale + l_rumple.scale, 
		   lambda=best[['lambda_coef']], alpha= best[['alpha_coef']]),
		   biotic = bioticStruct(lambda=best[['lambda_cov']], alpha=best[['alpha_cov']], on_diag=F, inverse = FALSE),
		   spatial = linear(data=as.matrix(XY), ~0+UTM_E:UTM_N, 
					   lambda =best[['lambda_spatial']] , alpha = best[['alpha_spatial']]),
		   learning_rate = 0.01, step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L
		   )
	
loss=unlist(model$logLik)
history=model$history
p = getSe(model)
result = list( beta = coef(model), sigma = getCov(model),loss=loss,history=model$history, p=p, logLik=logLik(model))
	
# saveRDS(result, file = here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_result_s1_m1_quasiP_no4_lidar_0.1.0.RDS') )
# saveRDS(model, here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_model_s1_m1_quasiP_no4_lidar_0.1.0.RDS') )
	
model = readRDS( here('results', 'sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_model_s1_m1_quasiP_no4_lidar_0.1.0.RDS') )
result = readRDS( here('results', 'sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_result_s1_m1_quasiP_no4_lidar_0.1.0.RDS'))
	

# . plot variation partition
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000', 'sjsdm_s1_m1_qp_no4',  'model_lidar','vari.par_sjSDM_s1_m1_qp_no4_lidar_0.1.0.pdf'), width=5, height=5)
	
imp =importance(model)
names(imp)
	
plot(imp, cex=.8) #from s-jSDM package
title(paste('variation partition, ',ncol(result$sigma), ' OTUs', sep=''))
	
dev.off()
	
# differentiate order & family
impD = data.frame(OTU=imp$names, imp$res$total)
str(impD)
	
impD$order = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[2])
table(impD$order)
	
impD$family = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[3])
sort(table(impD$family))
	
a = as.data.frame(impD %>% count(order))
a = subset(a, n>9)
impD = inner_join(impD, a, by=c('order'='order'))
str(impD)
	
impD$order = mgsub(impD$order, c('Diptera', 'Coleoptera', 'Lepidoptera', 'Hemiptera', 'Hymenoptera'), c('Fly', 'Beetle & Weevil', 'Butterfly & Moth', 'Bug', 'Wasp & Bee'))
impD$family = mgsub(impD$family, c('Muscidae', 'Geometridae', 'Cecidomyiidae', 'Ichneumonidae'), c('Fly', 'Moth', 'Gall midge', 'Ichneumon wasp'))
str(impD)
	
a = as.data.frame(impD %>% count(family))
a = subset(a, n>9)
impD2 = inner_join(impD, a, by=c('family'='family'))
impD2 = subset(impD2, family!='BOLD' & family != 'NA')
str(impD2)
	
table(impD2$order)
table(impD2$family)
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000', 'sjsdm_s1_m1_qp_no4', 'model_lidar', 'vari.family_sjSDM_s1_m1_qp_no4_lidar_0.1.0.pdf'), width=8, height=5)
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Family (occurrence>9)'
vp.plot(data=impD2, ind=6, textM=titleText, textL=legendText)
	
dev.off()
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000', 'sjsdm_s1_m1_qp_no4',  'model_lidar', 'vari.order_sjSDM_s1_m1_qp_no4_lidar_0.1.0.pdf'), width=8, height=5)
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Order (occurrence>9)'
vp.plot(data=impD2, ind=5, textM=titleText, textL=legendText)
	
dev.off()
	
summary(cut(impD2$env, breaks = seq(0,1,length.out = 12)))
	

# anova plot
an = anova(model,cv = F)
#saveRDS(an, file = here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','anova_sjSDM_s1_m1_quasiP_no4_lidar_0.1.0.RDS'))
an = readRDS(here('results','sjsdm-model-RDS','sjsdm-0.1.0.9000','anova_sjSDM_s1_m1_quasiP_no4_lidar_0.1.0.RDS'))
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000', 'sjsdm_s1_m1_qp_no4',  'model_lidar','anova_sjSDM_s1_m1_qp_no4_lidar_0.1.0.pdf'), width=5, height=5)
	
plot(an)
title(paste('analysis of variance, ',ncol(result$sigma), ' OTUs', sep=''))
	
dev.off()
	

```



```{r sjsdm-model-lidar (0/1)}
# . sjSDM-version = ‘0.1.0.9000’
# 0/1, >=5
all.data = dataI.1.present5
otu.data = as.data.frame(select(all.data, contains('__'))) 
	
env.data = read.csv(here('results', 'crossvalidation', 'results_20200904_5minocc_lidar_pa_loocv', 'data_20200904_5minocc_lidar', 'scale.env1.csv'), header=T, sep=',')
names(env.data)
	
scale.env1 = as.data.frame(select(all.data, names(env.data)))
str(scale.env1)
	
XY = as.data.frame(select(all.data, starts_with('UTM_')) %>% scale)
#XY = as.data.frame(select(all.data, starts_with('UTM_')) )
	
par(mfrow=c(2,2))
	
hist(all.data$UTM_E)
hist(all.data$UTM_N)
hist(XY$UTM_E)
hist(XY$UTM_N)
	
# load tune result with only best (adagpu)
best = readRDS(here('results', 'crossvalidation', 'results_20200904_5minocc_lidar_pa_loocv', "sjsdm_tune_results_HJA_20200904_bestonly.RDS"))
best
	
# run sjSDM
model = sjSDM( Y = as.matrix(otu.data),
		   env = linear(data=as.matrix(scale.env1), 
#		   formula = ~ elevation.scale + l_Cover_2m_max.scale + l_Cover_2m_4m.scale + l_Cover_4m_16m.scale + l_p25.scale + l_p95.scale + l_rumple.scale, 
		   lambda=best[['lambda_coef']], alpha= best[['alpha_coef']]),
		   biotic = bioticStruct(lambda=best[['lambda_cov']], alpha=best[['alpha_cov']], on_diag=F, inverse = FALSE),
		   spatial = linear(data=as.matrix(XY), ~0+UTM_E:UTM_N, 
					   lambda =best[['lambda_spatial']] , alpha = best[['alpha_spatial']]),
		   learning_rate = 0.01, step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L
		   )
	
loss=unlist(model$logLik)
history=model$history
p = getSe(model)
result = list(beta = coef(model), sigma = getCov(model),loss=loss,history=model$history, p=p, logLik=logLik(model))
	
#saveRDS(result, file = here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_result_s1_m1_present_no4_lidar_0.1.0.RDS') )
#saveRDS(model, file = here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_model_s1_m1_present_no4_lidar_0.1.0.RDS') )
	
rm(model)
	
model = readRDS(here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_model_s1_m1_present_no4_lidar_0.1.0.RDS'))
result = readRDS(here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','s-jSDM_result_s1_m1_present_no4_lidar_0.1.0.RDS') )
str(result)
	
summary(model)
	
summary.p=summary(result$p)
str(summary.p)
	
# .. model with space
# . plot variation partition
imp =importance(model)
names(imp)
names(imp$res)
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', 'vari.par_sjSDM_s1_m1_pa_no4_0.1.0.pdf'), width=5, height=5)
	
plot(imp, cex=.8) #from s-jSDM package
title(paste('variation partition, ',ncol(result$sigma), ' OTUs', sep=''))
	
dev.off()
	
# differentiate order & family
impD = data.frame(OTU=imp$names, imp$res$total)
str(impD)
	
impD$order = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[2])
table(impD$order)
# Diptera (flies), Coleoptera (beetles, weevils), Lepidoptera (butterflies, moths), Hymenoptera (wasp,  bees, ants), Hemiptera (bugs)
impD$family = sapply(strsplit(sapply(str_split(impD$OTU, '__'), function(a) a[2]), '_'), function(a) a[3])
sort(table(impD$family))
# 24 not defined
# Muscidae (house flies), Geometridae (geometer moths), Cecidomyiidae (gall midges), Ichneumonidae (ichneumon wasps)
	
a = as.data.frame(impD %>% count(order))
a = subset(a, n>9)
impD = inner_join(impD, a, by=c('order'='order'))
str(impD)
	
impD$order = mgsub(impD$order, c('Diptera', 'Coleoptera', 'Lepidoptera', 'Hemiptera', 'Hymenoptera'), c('Fly', 'Beetle & Weevil', 'Butterfly & Moth', 'Bug', 'Wasp & Bee'))
impD$family = mgsub(impD$family, c('Muscidae', 'Geometridae', 'Cecidomyiidae', 'Ichneumonidae'), c('Fly', 'Moth', 'Gall midge', 'Ichneumon wasp'))
str(impD)
	
a = as.data.frame(impD %>% count(family))
a = subset(a, n>9)
impD2 = inner_join(impD, a, by=c('family'='family'))
impD2 = subset(impD2, family!='BOLD' & family != 'NA')
str(impD2)
	
table(impD2$order)
table(impD2$family)
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', 'vari.family_sjSDM_s1_m1_pa_no4_0.1.0.pdf'), width=8, height=5)
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Family (occurrence>9)'
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
#summary(cut(impD$spatial, breaks = seq(0,1,length.out = 12)))
#summary(cut(impD$biotic, breaks = seq(0,1,length.out = 12)))
#summary(cut(impD$env, breaks = seq(0,1,length.out = 12)))
	
vp.plot(data=impD2, ind=6, textM=titleText, textL=legendText)
	
dev.off()
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model', 'vari.order_sjSDM_s1_m1_pa_no4_0.1.0.pdf'), width=8, height=5)
	
titleText=paste('variation partition, ',nrow(data), ' OTUs', sep='')
legendText = 'Order (occurrence>9)'
vp.plot(data=impD2, ind=5, textM=titleText, textL=legendText)
	
dev.off()
	
# anova graph
an = anova(model,cv = F)
#saveRDS(an, file = here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','sjSDM_anova_s1_m1_present_no4_0.1.0.RDS'))
an = readRDS(here('results','sjsdm-model-RDS', 'sjsdm-0.1.0.9000','sjSDM_anova_s1_m1_present_no4_0.1.0.RDS'))
print(an)
	
#pdf(here('results','sjsdm-graph','sjsdm-0.1.0.9000','sjsdm_s1_m1_pa_no4', 'lidar_model','anova_sjSDM_s1_m1_pa_no4_0.1.0.pdf'), width=5, height=5)
	
plot(an)
title(paste('analysis of variance, ',ncol(result$sigma), ' OTUs', sep=''))
	
dev.off()
	
```

