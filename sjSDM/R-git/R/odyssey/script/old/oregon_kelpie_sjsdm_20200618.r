---
---

```{r setup}
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv', 'tidyverse', 'scatterplot3d', 'gridBase', 'grid', 'ggcorrplot','here'), library, character.only=T)
	
lapply(c('reticulate','sjSDM'), library, character.only=T)
	
here::here()
getwd()
source(here("R", "source", "sjsdm_function.r"))
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
```

```{r convert data (go directly to 'load files'), eval=FALSE, include=FALSE}
# for backup, converted data are saved

# format OTU data, combine with lidar data 
#dy not lidar data. it's landsat data.  there will be lidar data later, so must name it correctly
# (remote sensing data)
mulspec.env = read.csv(here('..','..','HJA_scripts','10_eo_data','biodiversity_site_info_multispectral_2020-04-13.txt'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(mulspec.env)
	
# ('HJA_analyses_Kelpie/Kelpie_maps' folder) 
otu.env1.noS = read.csv(here('..','..','Kelpie_maps', 'outputs_minimap2_20200221_F2308_f0x2_q48_kelpie20200214_vsearch97','sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
otu.env1.spike = read.csv(here('..','..','Kelpie_maps', 'outputs_minimap2_20200221_F2308_f0x2_q48_kelpie20200214_vsearch97', 'sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
# 1173-26, more otus
	
names(otu.env1.noS)[1:26] == names(otu.env1.spike)[1:26]
	
names(otu.env1.noS)[1:26]=c('SiteName', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
names(otu.env1.spike)[1:26]=names(otu.env1.noS)[1:26]
	
str(otu.env1.spike[,1:27])
	

# ........ format mulspec.env .......
sort(mulspec.env$SiteName) == sort(unique(otu.env1.spike$SiteName))
	
#NDVI - normalized difference vegetation index (calculated using bands 4 and 5): NDVI = (NearIR-Red)/(NearIR+Red)
#       these values should range between -1 and 1. Values in these columns should be divided by 1000
#EVI - enhanced vegetation index (calculated using bands 4, 5, and 2):  2.5 * ((Band 5 – Band 4) / (Band 5 + 6 * Band 4 – 7.5 * Band 2 + 1))
#      the values in these columns should be divided by 1000
names(mulspec.env)   #sort()
# 4 NDVI, 4 EVI

data.frame(mulspec.env$NDVI_20180717, mulspec.env$EVI_20180818)
	
a = mulspec.env %>% select(starts_with('NDVI'), starts_with('EVI')) %>% rename(nor.NDVI_20180717=1, nor.EVI_20180717=5, nor.NDVI_20180726=2, nor.EVI_20180726=6, nor.NDVI_20180802=3, nor.EVI_20180802=7, nor.NDVI_20180818=4, nor.EVI_20180818=8)/1000
#13,14,27,28,41,42,55,56
	
a[,c(1,3,5,7)]
# mean of all 4 NDVI, EVI
mulspec.env= cbind(mulspec.env,a)
str(mulspec.env)
	
mulspec.env$mean.NDVI = base::rowMeans(select(mulspec.env, starts_with('nor.NDVI_')))
mulspec.env$mean.EVI = base::rowMeans(select(mulspec.env, starts_with('nor.EVI_')))
mulspec.env[,c(60,62,64,66,68)]
head(mulspec.env[,c(61,63,65,67,69)])
	
# ... calculate mean B, G, W values
# B -> Tasseled cap brightness; "G" -> Tasseled cap greenness; "W" -> Tasseled cap wetness
mulspec.env$mean.bright = base::rowMeans(select(mulspec.env, starts_with('B_')))
mulspec.env$mean.green = base::rowMeans(select(mulspec.env, starts_with('G_')))
mulspec.env$mean.wet = base::rowMeans(select(mulspec.env, starts_with('W_')))
	
# pdf(here('R','graph','describe_spec_indices.pdf'), height=5, width=5)
	
range(select(mulspec.env, starts_with('B_'), starts_with('G_'), starts_with('W_')))
plot(1:96, mulspec.env$B_20180818,ylim=c(-1519,3505),type='l', main='solid - brightness, dash - wetness, dotted - greenness', col='black')
lines(1:length(mulspec.env$SiteName), mulspec.env$B_20180717, col='black')
lines(1:length(mulspec.env$SiteName), mulspec.env$B_20180726, col='black')
lines(1:length(mulspec.env$SiteName), mulspec.env$B_20180802, col='black')
	
lines(1:length(mulspec.env$SiteName), mulspec.env$W_20180818, lty=2, col='black')
lines(1:length(mulspec.env$SiteName), mulspec.env$W_20180717, lty=2, col='red')
lines(1:length(mulspec.env$SiteName), mulspec.env$W_20180726, lty=2, col='blue')
lines(1:length(mulspec.env$SiteName), mulspec.env$W_20180802, lty=2, col='green')
	
lines(1:length(mulspec.env$SiteName), mulspec.env$G_20180818, lty=3, col='green')
lines(1:length(mulspec.env$SiteName), mulspec.env$G_20180717, lty=3, col='green')
lines(1:length(mulspec.env$SiteName), mulspec.env$G_20180726, lty=3, col='green')
lines(1:length(mulspec.env$SiteName), mulspec.env$G_20180802, lty=3, col='green')
	
# dev.off()
	

# ... explore OTU reads ...
names(otu.env1.spike)[1:27]
	
hist(otu.env1.spike[,30])
sort(unique(otu.env1.noS[,30]))
	
which(is.na(otu.env1.noS[,28]))#:dim(otu.env1.noS)[2]
which(is.na(otu.env1.noS[,29]))
# 181 237
	
table(is.na(otu.env1.noS[c(181,237),27:1173]))
	
table(is.na(otu.env1.spike[c(181,237),27:1173]))
	
otu.env1.noS$SiteName[c(181,237)]
	
# delete row 181, 237 -> "HOBO-040", "290415" 
dim(otu.env1.noS)
	
otu.env1.noS = otu.env1.noS[-c(181,237),]
	
dim(otu.env1.spike)
	
otu.env1.spike = otu.env1.spike[-c(181,237),]
	

# ....... scale variables .......
names(otu.env1.spike)[15:22]
	
a = select(otu.env1.spike,15:22)%>%rename(elevation.scale=1,canopy.ht.scale=2,min.T.scale=3, max.T.scale=4, precipitation.scale=5, metre.road.scale=6, metre.stream.scale=7, yrs.disturb.min.scale=8)%>% scale()
	
otu.env1.spike = cbind(otu.env1.spike[,1:26], data.frame(a), otu.env1.spike[,27:dim(otu.env1.spike)[2]])
dim(otu.env1.spike)
	
names(otu.env1.spike)[26:34]
	
otu.env1.noS = dplyr::left_join(otu.env1.noS, select(otu.env1.spike,'site_trap_period',ends_with('.scale')), by=c('site_trap_period', 'site_trap_period'), copy=F)
dim(otu.env1.noS)
names(otu.env1.noS)[1174:1194]
	
otu.env1.noS = otu.env1.noS[,c(1:26,1174:1181,27:1173)]
str(otu.env1.noS[,1:34])
	

# ..... combine multispectral-data .....
str(mulspec.env)
	
par(mfrow=c(1,2))
	
hist(mulspec.env$mean.NDVI)
hist(mulspec.env$mean.EVI)
	
par(mfrow=c(1,3))
	
hist(mulspec.env$mean.green)
hist(mulspec.env$mean.bright)
hist(mulspec.env$mean.wet)
	
names(mulspec.env)[68:69]
	
a = select(mulspec.env,"mean.NDVI", "mean.EVI",'mean.green','mean.bright','mean.wet') %>% rename(mean.NDVI.scale=1,mean.EVI.scale=2, mean.green.scale=3,mean.bright.scale=4,mean.wet.scale=5) %>% scale()
mulspec.env = cbind(mulspec.env,data.frame(a))
	
par(mfrow=c(1,2))
	
hist(mulspec.env$mean.NDVI.scale)
hist(mulspec.env$mean.EVI.scale)
	
par(mfrow=c(1,3))
	
hist(mulspec.env$mean.green.scale)
hist(mulspec.env$mean.bright.scale)
hist(mulspec.env$mean.wet.scale)
	
# row 181, 237 -> "HOBO-040", "290415" in OTU datasets
sort(mulspec.env$SiteName) == sort(unique(otu.env1.spike$SiteName))
	
names(mulspec.env)[c(1,70:71)]
	
otu.env1.spike = dplyr::left_join(otu.env1.spike, dplyr::select(mulspec.env,'SiteName',ends_with('.scale')), by=c('SiteName', 'SiteName'), copy=F) 
str(otu.env1.spike[,1:37])
names(otu.env1.spike[,1182:1186])
	
otu.env1.spike = otu.env1.spike[,c(1:34,1182:1186,35:1181)]
dim(otu.env1.spike)
	
otu.env1.noS = dplyr::left_join(otu.env1.noS, dplyr::select(mulspec.env,'SiteName',ends_with('.scale')), by=c('SiteName', 'SiteName'), copy=F) 
dim(otu.env1.noS)
	
str(otu.env1.noS[,1:37])
	
otu.env1.noS = otu.env1.noS[,c(1:34,1182:1186,35:1181)]
dim(otu.env1.noS)
	
# write.table(otu.env1.noS, here('kelpie','formatted_data','mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	
# write.table(otu.env1.spike, here('kelpie','formatted_data','mulspec2_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	
# write.table(mulspec.env, here('kelpie','rs_data','sjsdm_biodiversity2_site_info_multispectral_2020-04-13.csv'), row.names= F, sep=',')
	

```

```{r convert data2 (go directly to 'load files'), eval=FALSE, include=FALSE}
# for backup, data are already saved
# convert to present/rel.abun according to formatted data 
# ( continue with 'convert data')

# (remote sensing data, load data formatted for sjsdm)
mulspec.env = read.csv(here('kelpie','rs_data','sjsdm_biodiversity2_site_info_multispectral_2020-04-13.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(mulspec.env)
	
# ('formatted_data' folder, load data formatted for sjsdm) 
otu.env1.noS = read.csv(here('kelpie','formatted_data','mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
otu.env1.spike = read.csv(here('kelpie','formatted_data','mulspec2_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
	

# ..... convert to presence/absence .....
# . no spike
names(otu.env1.noS)[1:37]
	
a = select(otu.env1.noS, contains('__'))
dim(a)
dim(otu.env1.noS)
# 1147 + 36 = 1183
	
# !!! env.variables cannot have '__' in their names !!!
a = data.frame(lapply(data.frame(select(otu.env1.noS, contains('__'))>0), as.numeric))
otu.env1.noS.present = cbind(select(otu.env1.noS, -contains('__')), a)
	

#DY old code using base::scale()
# a = data.frame(sapply(otu.env1.noSlogLik(sjsdm.s1_m1_spike.relAbun)[,-c(1:36)],function(x) scale(x,center=F)))

#DY rescaling to quasi-probabilities, following Max's advice
# sjsdm is a probit model, X needs to be [0,1]
a = data.frame(sapply(select(otu.env1.noS, contains('__')),function(x) scales::rescale(log(x+0.001))))
	
# with log(), data points go towards either 0/1. is this better for a probit model ??? 
par(mfrow=c(1,2))
hist(a[,115])
hist(otu.env1.noS[,115+36])
	
otu.env1.noS.rel.abun = cbind(select(otu.env1.noS, -contains('__')), a)
	

# . spike
a = data.frame(lapply(data.frame(select(otu.env1.spike, contains('__'))>0), as.numeric))
table(a[,111])
otu.env1.spike.present = cbind(select(otu.env1.spike, -contains('__')), a)
names(otu.env1.spike.present)[1:40]
	
# hist(otu.env1.spike$R200_218__Insecta_Diptera_Syrphidae_Blera_Blera_scitula_BOLD_ABY7981_size.595)


#DY rescaling to quasi-probabilities, following Max's advice
a = data.frame(sapply(select(otu.env1.spike, contains('__')), function(x) scales::rescale(log(x+0.001))))
par(mfrow=c(1,2))
hist(a[,115])
hist(otu.env1.spike[,115+36])
	
otu.env1.spike.rel.abun = cbind(select(otu.env1.spike, -contains('__')), a)
range(otu.env1.spike.rel.abun[,-c(1:39)])
table(otu.env1.spike.rel.abun[,400])
	

#DY i usually comment out write.table() commands, to prevent running them incorrectly (e.g. when debugging)
# write.table(otu.env1.noS.rel.abun, here('kelpie','formatted_data','relAbun_mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
# write.table(otu.env1.noS.present, here('kelpie','formatted_data','present_mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	
# 
# write.table(otu.env1.spike.present, here('kelpie','formatted_data','present2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
# write.table(otu.env1.spike.rel.abun, here('kelpie','formatted_data','relAbun2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	

```


```{r load files, echo=FALSE}
# kelpie, remote sensing data 
mulspec.env = read.csv(here('kelpie','rs_data','sjsdm_biodiversity2_site_info_multispectral_2020-04-13.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(mulspec.env)
	
# ('formatted_data' folder, load data formatted for sjsdm, present) 
#otu.env1.noS.present = read.csv(here('kelpie','formatted_data','present_mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
#otu.env1.noS.rel.abun = read.csv(here('kelpie','formatted_data','relAbun_mulspec_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
#range(otu.env1.noS.rel.abun[,37:800])
	
# use only 'spike' data for analysis. for 0/1 it's the same as not spiked, but order of rows&cols is consistent with spike rel.abun 
otu.env1.present = read.csv(here('kelpie','formatted_data','present2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
otu.env1.rel.abun = read.csv(here('kelpie','formatted_data','relAbun2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
range(otu.env1.rel.abun[,40:800])
	
#names(otu.env1.spike.present)[40:47]==names(otu.env1.spike)[27:34]
#names(otu.env1.spike.rel.abun)[40:47]==names(otu.env1.spike)[27:34]
	
```

```{r explore remote sensing data }
# pdf(here('R','graph','describe_EVI_NDVI.pdf'), height=5, width=5)
	
range(mulspec.env[,60:67])
plot(1:96, mulspec.env[,61],ylim=c(0,3.3),type='l', main='solid - EVI, dash - NDVI')
lines(1:96, mulspec.env[,63], col='red')
lines(1:96, mulspec.env[,65], col='blue')
lines(1:96, mulspec.env[,67], col='green')
	
lines(1:96, mulspec.env[,60], lty=2, col='black')
lines(1:96, mulspec.env[,62], lty=2, col='red')
lines(1:96, mulspec.env[,64], lty=2, col='blue')
lines(1:96, mulspec.env[,66], lty=2, col='green')
	
# dev.off()
	
```

```{r subsets of data}
dataI.1.present = subset(otu.env1.present, session == 'S1' & trap == 'M1' )
	
dataII.1.present = subset(otu.env1.present, session == 'S2' & trap == 'M1' )
	
dataI.1.rel.abun = subset(otu.env1.rel.abun, session == 'S1' & trap == 'M1' )
	
dataII.1.rel.abun = subset(otu.env1.rel.abun, session == 'S2' & trap == 'M1' )
print(c(dim(dataI.1.rel.abun), dim(dataI.1.present)))
	
print(c(dim(dataII.1.rel.abun), dim(dataII.1.present)))
	
# .... if there's all zero OTUs ....
# .... S1
b = data.frame(otu=colnames(select(dataI.1.rel.abun, contains('__'))), zero=apply(select(dataI.1.rel.abun, contains('__')),2,sum)==0)
b$otu=as.character(b$otu)
dim(b)
	
dataI.1.rel.abun2 = dplyr::select(dataI.1.rel.abun, -contains('__'),b$otu[b$zero==F])
dim(dataI.1.rel.abun2)
dataI.1.rel.abun2[1:5,34:40]
	
b = data.frame(otu=colnames(select(dataI.1.present,contains('__'))), zero=apply(select(dataI.1.present,contains('__')),2,sum)==0)
b$otu=as.character(b$otu)
dim(b)
	
dataI.1.present2 = dplyr::select(dataI.1.present, -contains('__'),b$otu[b$zero==F])
dim(dataI.1.present2)
dataI.1.present2[1:5,34:40]
	
print(c(dim(dataI.1.rel.abun2), dim(dataI.1.present2)))
	

# ... session II
b = data.frame(otu=colnames(select(dataII.1.present,contains('__'))), zero=apply(select(dataII.1.present,contains('__')),2,sum)==0)
b$otu=as.character(b$otu)
table(b$zero)
	
dataII.1.present2 = dplyr::select(dataII.1.present, -contains('__'),b$otu[b$zero==F])
dim(dataII.1.present2)
dataII.1.present2[1:5,34:40]
	

b = data.frame(otu=colnames(select(dataII.1.rel.abun,contains('__'))), zero=apply(select(dataII.1.rel.abun,contains('__')),2,sum)==0)
b$otu=as.character(b$otu)
table(b$zero)
	
dataII.1.rel.abun2 = dplyr::select(dataII.1.rel.abun, -contains('__'),b$otu[b$zero==F])
dim(dataII.1.rel.abun2)
dataII.1.rel.abun2[1:5,34:40]
	
print(c(dim(dataII.1.rel.abun2), dim(dataII.1.present2)))
	

```

```{r sjsdm-model-spatial (rel.abun & 0/1)}
# ... session 1, Malaise I, present, no mulspec-data ...
# exploring spatial function
# . sjSDM-version = ‘0.0.2.9000’, but updated that spatial is available in 'sjSDM_cv'
scale.env1 = dataI.1.present2[ ,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale')]
str(scale.env1)
	
otu.data = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
par(mfrow=c(2,2))
	
hist(dataI.1.present2$UTM_E)
hist(dataI.1.present2$UTM_N)
	
plot(dataI.1.present2$UTM_E, dataI.1.present2$UTM_N)
	
XY = as.data.frame(select(dataI.1.present2, starts_with('UTM_')) %>% scale)
#table(XY$UTM_E == scale(dataI.1.present2$UTM_E))
hist(XY$UTM_E)
hist(XY$UTM_N)
	
hyper.s1_m1_present.space = sjSDM_cv(
	Y = as.matrix(otu.data),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, 
#	tune_steps = 3L, iter = 10L,
	 tune_steps = 20L, iter = 60L, 
	spatial = linear(as.matrix(XY), ~UTM_E:UTM_N)
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
	
 saveRDS(hyper.s1_m1_present.space, file = here('R','result','s-jSDM_tune_s1_m1_present.space_V0.0.2.9000_2.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_present.space, perf = "AUC")
# still has error:
#Error in mgcv::gam(stats::as.formula(form), data = x) : 
#  Model has more coefficients than data
	
a = as.data.frame(summary(hyper.s1_m1_present.space))
best = a[which(a$'AUC_test'==max(a$'AUC_test')),]
	
	
best = hyper.s1_m1_present.space
	
 saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_present.space_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s1_m1_present.space
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_present.space)
	
sjsdm.s1_m1_present.space = sjSDM(
	Y = as.matrix(otu.data),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale + canopy.ht.scale+ min.T.scale + max.T.scale + precipitation.scale + metre.road.scale + metre.stream.scale + yrs.disturb.min.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	spatial = linear(data=as.matrix(XY), ~UTM_E:UTM_N
	, 
	lambda = best[["lambda_spatial"]], alpha = best[["alpha_spatial"]]
	) , 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
logLik = logLik(sjsdm.s1_m1_present.space)
#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_present.space)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
#saveRDS(sjsdm.s1_m1_present.space, file = here('R','result','model_sjSDM_s1_m1_present.space_V0.0.2.9000.RDS'))
	
result = list(beta = coef(sjsdm.s1_m1_present.space), sigma = getCov(sjsdm.s1_m1_present.space), history = sjsdm.s1_m1_present.space$history, p = p, logLik=logLik)
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_present.space_V0.0.2.9000.RDS')) 
	

	
```



```{r sjsdm-model-only-mulspec (rel.abun & 0/1)}
# ... session 1, Malaise I, spike, 0/1 ...
# ... with multispectral environmental data
# . sjSDM-version = ‘0.0.2.9000’
scale.env1 = data.frame(select(dataI.1.present2, starts_with('mean')))
str(scale.env1)
	
otu.data = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
hyper.s1_m1_present.spec = sjSDM_cv(
	Y = as.matrix(otu.data),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
 saveRDS(hyper.s1_m1_present.spec, file = here('R','result','s-jSDM_tune_s1_m1_spike.present.spec_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_present.spec, perf = "AUC")
	
# saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.present.spec_V0.0.2.9000.RDS') )
	
# best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.present.spec_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s1_m1_present.spec
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_present.spec)
	
sjsdm.s1_m1_present.spec = sjSDM(
	Y = as.matrix(otu.data),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ mean.NDVI.scale +mean.EVI.scale +mean.green.scale + mean.bright.scale + mean.wet.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
	
saveRDS(sjsdm.s1_m1_present.spec, file = here('R','result','model_sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
	
sjsdm.s1_m1_present.spec = readRDS(here('R','result','model_sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
	
summary(result$p)$coefmat
	
a = 
sjSDM::plot(sjSDM::importance(sjsdm.s1_m1_present.spec))
	
an = anova(sjsdm.s1_m1_present.spec,cv = 2L)
plot(an)


logLik = logLik(sjsdm.s1_m1_present.spec)
#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_present.spec)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_present.spec), sigma = getCov(sjsdm.s1_m1_present.spec), history = sjsdm.s1_m1_present.spec$history, p = p, logLik=logLik)
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS')) 
	



# ... session 1, Malaise I, spike, rel.abun ...
# ... use multispectral environmental data
# . sjSDM-version = ‘0.0.2.9000’
names(dataI.1.rel.abun2[,c(1:3,24:39)])
# [1] "SiteName"              "UTM_E"                 "UTM_N"                
# [4] "trap"                  "session"               "site_trap_period"     
# [7] "elevation.scale"       "canopy.ht.scale"       "min.T.scale"          
#[10] "max.T.scale"           "precipitation.scale"   "metre.road.scale"     
#[13] "metre.stream.scale"    "yrs.disturb.min.scale" "mean.NDVI.scale"      
#[16] "mean.EVI.scale"
scale.env = dataI.1.rel.abun2[,c(1:3,24:36)]
str(scale.env)
	
scale.env1 = data.frame(select(dataI.1.rel.abun2, starts_with("mean")))
str(scale.env1)
	
otu.data = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
hyper.s1_m1_relAbun.spec = sjSDM_cv(
	Y = as.matrix(otu.data),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
saveRDS(hyper.s1_m1_relAbun.spec, file = here('R','result','s-jSDM_tune_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_relAbun.spec, perf = "AUC")
	
saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s1_m1_relAbun.spec
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_relAbun.spec)
	
sjsdm.s1_m1_relAbun.spec = sjSDM(
	Y = as.matrix(otu.data),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ mean.NDVI.scale +mean.EVI.scale +mean.green.scale + mean.bright.scale + mean.wet.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)

#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_relAbun.spec)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_relAbun.spec), sigma = getCov(sjsdm.s1_m1_relAbun.spec), history = sjsdm.s1_m1_relAbun.spec$history, p = p, logLik=logLik(sjsdm.s1_m1_relAbun.spec))
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS')) 
	
```


```{r sjsdm-analyse-correlation-only-mulspec (relAbun)}
# ... session 1, Malaise I, spike, rel.abun, only spec-data ...
# .. load 'sjSDM' result
best =  readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
	
#tune = readRDS(here('R','result','s-jSDM_tune_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
#plot(tune, perf = "AUC")
#rm(tune)
	
summary(result)
plot(result$history)
	
# ... correlation plot ...
dim(result$sigma)# cov <- result$sigma 
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
	
otu.table = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
str(otu.table[,1:5])
	
# extract min & max spp pairs of correlation
#source(here("R", "source", "sjsdm-analyse-functions.r"))
extract.minmax.cor(otu=otu.table, sigma=co.env.spp)
	
minmax = my.output[[2]]
	
#write.table(minmax, file=here::here('correlation_value_s1_m1_spike_relAbun_spec.csv'), row.names=F, sep=', ')
rm(minmax, my.output)
	

rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_oSpec','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5),sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	

```


```{r model-analyse-polygon-only-mulspec (relAbun)}
# ... session 1, Malaise I, spike, relAbun, only spec-data ...
# ..... Polygon Drawing .....
# . species association .
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
number=10
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
# variables needed for plotting. can skip. good for checking
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
summary(cuts)
cols = viridis::magma(10) 
# variables needed for plotting. can skip. good for checking
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_oSpec','species_covariance_circular.pdf'), height=8, width=8)
	
version.text = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep='')
otu.text = "sum of relAbun"
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl)
	
#dev.off()
	
```



```{r model-analyse-polygon-envir-only-spec (relAbun)}
# ... session 1, Malaise I, spike, relAbun, only spec-data ...
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ mean.NDVI.scale+mean.EVI.scale + mean.green.scale + mean.bright.scale + mean.wet.scale
# 10 variables, excluding intercept
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun.spec_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
beta = as.matrix(data.frame(result$beta)[,2:6])
effects = apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
str(effect_comb)
# variables index & value which has biggest coefficient
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
	
version.text = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep='')
otu.text = "relAbun. spp."
	
evnames =c('NDVI', 'EVI','green', 'brightness', 'wetness')
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_oSpec','max_environ_spp_cov.pdf'), height=8, width=8)
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle.env(sigma=sigma, version.text=version.text, evnames=evnames, otu.text=otu.text)
	
#dev.off()
	
```



```{r sjsdm-analyse-correlation-only-mulspec (0/1)}
# ... session 1, Malaise I, spike, present, only spec-data ...
# .. load 'sjSDM' result
result = readRDS(here('R','result','sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
	
summary(result)
plot(result$history)
	
# ... correlation plot ...
dim(result$sigma)# cov <- result$sigma 
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
	
otu.table = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
str(otu.table[,1:5])
	
# extract min & max spp pairs of correlation
source(here("R", "source", "sjsdm-analyse-functions.r"))
extract.minmax.cor(otu=otu.table, sigma=co.env.spp)
	
minmax = my.output[[2]]
	
#write.table(minmax, file=here::here('correlation_value_s1_m1_spike_present_spec.csv'), row.names=F, sep=', ')
rm(minmax, my.output)
	

rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_oSpec','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5),sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	

```

```{r model-analyse-polygon-only-mulspec (0/1)}
# ... session 1, Malaise I, spike, present, only spec-data ...
# ..... Polygon Drawing .....
# . species association .
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
number=10
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
# variables needed for plotting. can skip. good for checking
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
# variables needed for plotting. can skip. good for checking
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_oSpec','species_covariance_circular.pdf'), height=8, width=8)
	
version.text = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep='')
otu.text = "sum of presence"
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl)
	
#dev.off()
	
```


```{r model-analyse-polygon-envir-only-spec (0/1)}
# ... session 1, Malaise I, spike, present, only spec-data ...
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ mean.NDVI.scale+mean.EVI.scale + mean.green.scale + mean.bright.scale + mean.wet.scale
# 10 variables, excluding intercept
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.present.spec_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
beta = as.matrix(data.frame(result$beta)[,2:6])
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
str(effect_comb)
# variables index & value which has biggest coefficient
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
	
version.text = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep='')
otu.text = "presence. spp."
	
evnames =c('NDVI', 'EVI','green', 'brightness', 'wetness')
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_oSpec','max_environ_spp_cov.pdf'), height=8, width=8)
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle.env(sigma=sigma, version.text=version.text, evnames=evnames, otu.text=otu.text)
	
#dev.off()
	
```


```{r sjsdm-model-SII (rel.abun & 0/1) !!rerun new.data}
# models are already saved. go directly to 'model-analyse'

# according to nmds, precipitation, elevation, yrs.disturb.min, old.growth.str, T, canopy.height can be used for now

# ... session 2, Malaise I, spike, rel.abun ...
# . sjSDM-version = ‘0.0.2.9000’
names(dataII.1.rel.abun2[,c(1:3,24:36)])
# [1] "SiteName"              "UTM_E"                 "UTM_N"                
# [4] "trap"                  "session"               "site_trap_period"     
# [7] "elevation.scale"       "canopy.ht.scale"       "min.T.scale"          
#[10] "max.T.scale"           "precipitation.scale"   "metre.road.scale"     
#[13] "metre.stream.scale"    "yrs.disturb.min.scale" "mean.NDVI.scale"      
#[16] "mean.EVI.scale"
scale.env = dataII.1.rel.abun2[,c(1:3,24:36)]
str(scale.env)
	
scale.env1 = dataII.1.rel.abun2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale')]
str(scale.env1)
	
hyper.s2_m1_spike.relAbun = sjSDM_cv(
	Y = as.matrix(as.data.frame(dplyr::select(dataII.1.rel.abun2, contains('__')))),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
 saveRDS(hyper.s2_m1_spike.relAbun, file = here('R','result','s-jSDM_tune_s2_m1_spike.relAbun_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s2_m1_spike.relAbun, perf = "AUC")
	
saveRDS(best, file = here('R','result','s-jSDM_tune_best_s2_m1_spike.relAbun_V0.0.2.9000.RDS'))
	
## print overall results:
#hyper.s1_m1_spike.relAbun
## summary (mean values over CV for each tuning step)
#summary(hyper.s1_m1_spike.relAbun)
	
sjsdm.s2_m1_spike.relAbun = sjSDM(
	Y = as.matrix(as.data.frame(dplyr::select(dataII.1.rel.abun2, contains('__')))),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)

#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s2_m1_spike.relAbun)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s2_m1_spike.relAbun), sigma = getCov(sjsdm.s2_m1_spike.relAbun), history = sjsdm.s2_m1_spike.relAbun$history, p = p, logLik=logLik(sjsdm.s2_m1_spike.relAbun))
saveRDS(result, file = here('R','result','sjSDM_s2_m1_spike.relAbun_V0.0.2.9000.RDS')) 
	

# ... session 2, Malaise I, spike, 0/1 ...
# . sjSDM-version = ‘0.0.2.9000’
scale.env1 = dataII.1.present2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale')]
str(scale.env1)
	
hyper.s2_m1_spike.present = sjSDM_cv(
	Y = as.matrix(as.data.frame(dplyr::select(dataII.1.present2, contains('__')))),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
 saveRDS(hyper.s2_m1_spike.present, file = here('R','result','s-jSDM_tune_s2_m1_spike.present_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s2_m1_spike.present, perf = "AUC")
	
 saveRDS(best, file = here('R','result','s-jSDM_tune_best_s2_m1_spike.present_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s2_m1_spike.present
# summary (mean values over CV for each tuning step)
summary(hyper.s2_m1_spike.present)
	
sjsdm.s2_m1_spike.present = sjSDM(
	Y = as.matrix(as.data.frame(dplyr::select(dataII.1.present2, contains('__')))),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s2_m1_spike.present)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s2_m1_spike.present), sigma = getCov(sjsdm.s2_m1_spike.present), history = sjsdm.s2_m1_spike.present$history, p = p, logLik=logLik(sjsdm.s2_m1_spike.present))
 saveRDS(result, file = here('R','result','sjSDM_s2_m1_spike.present_V0.0.2.9000.RDS')) 
```	

```{r sjsdm-model-with-mulspec2 (rel.abun & 0/1)}
# ... session 1, Malaise I, 0/1 ...
# ... with both types of environmental data, more mulspec variables
# . sjSDM-version = ‘0.0.2.9000’
scale.env1 = dataI.1.spike.present2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale','mean.NDVI.scale','mean.EVI.scale', 'mean.green.scale', 'mean.bright.scale', 'mean.wet.scale')]
str(scale.env1)
	
otu.table = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
hyper.s1_m1_spike.present.W.spec2 = sjSDM_cv(
	Y = as.matrix(otu.table),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
 saveRDS(hyper.s1_m1_spike.present.W.spec2, file = here('R','result','s-jSDM_tune_s1_m1_spike.present.W.spec2_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_spike.present.W.spec2, perf = "AUC")
	
saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.present.W.spec2_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s1_m1_spike.present.W.spec2
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_spike.present.W.spec2)
	
sjsdm.s1_m1_spike.present.W.spec2 = sjSDM(
	Y = as.matrix(otu.table),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale +mean.NDVI.scale +mean.EVI.scale + mean.green.scale + mean.bright.scale + mean.wet.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
logLik = logLik(sjsdm.s1_m1_spike.present.W.spec2)
#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_spike.present.W.spec2)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_spike.present.W.spec2), sigma = getCov(sjsdm.s1_m1_spike.present.W.spec2), history = sjsdm.s1_m1_spike.present.W.spec2$history, p = p, logLik=logLik)
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.present.W.spec2_V0.0.2.9000.RDS')) 
	
```

```{r sjsdm-model-with-mulspec (rel.abun & 0/1)}
# ... session 1, Malaise I, spike, 0/1 ...
# ... with both types of environmental data
# . sjSDM-version = ‘0.0.2.9000’
scale.env1 = dataI.1.spike.present2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale','mean.NDVI.scale','mean.EVI.scale')]
str(scale.env1)
	
hyper.s1_m1_spike.present.W.spec = sjSDM_cv(
	Y = as.matrix(as.data.frame(dplyr::select(dataI.1.present2, contains('__')))),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
# saveRDS(hyper.s1_m1_spike.present.W.spec, file = here('R','result','s-jSDM_tune_s1_m1_spike.present.W.spec_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_spike.present.W.spec, perf = "AUC")
	
 #saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.present.W.spec_V0.0.2.9000.RDS') )
	
# print overall results:
hyper.s1_m1_spike.present.W.spec
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_spike.present.W.spec)
	
sjsdm.s1_m1_spike.present.W.spec = sjSDM(
	Y = as.matrix(as.data.frame(dplyr::select(dataI.1.present2, contains('__')))),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale +mean.NDVI.scale +mean.EVI.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
logLik = logLik(sjsdm.s1_m1_spike.present.W.spec)
#summary(model)
# calculate post-hoc p-values:
#p = getSe(sjsdm.s1_m1_spike.present.W.spec)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_spike.present.W.spec), sigma = getCov(sjsdm.s1_m1_spike.present.W.spec), history = sjsdm.s1_m1_spike.present.W.spec$history, p = p, logLik=logLik)
# saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.present.W.spec_V0.0.2.9000.RDS')) 
	

# ... session 1, Malaise I, spike, rel.abun ...
# ... use both types of environmental data
# . sjSDM-version = ‘0.0.2.9000’
names(dataI.1.rel.abun2[,c(1:3,24:39)])
# [1] "SiteName"              "UTM_E"                 "UTM_N"                
# [4] "trap"                  "session"               "site_trap_period"     
# [7] "elevation.scale"       "canopy.ht.scale"       "min.T.scale"          
#[10] "max.T.scale"           "precipitation.scale"   "metre.road.scale"     
#[13] "metre.stream.scale"    "yrs.disturb.min.scale" "mean.NDVI.scale"      
#[16] "mean.EVI.scale"
scale.env = dataI.1.rel.abun2[,c(1:3,24:36)]
str(scale.env)
	
scale.env1 = dataI.1.rel.abun2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale', "mean.NDVI.scale", "mean.EVI.scale")]
	
otu.data = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
hyper.s1_m1_spike.relAbun.W.spec = sjSDM_cv(
	Y = as.matrix(otu.data),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
saveRDS(hyper.s1_m1_spike.relAbun.W.spec, file = here('R','result','s-jSDM_tune_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_spike.relAbun.W.spec, perf = "AUC")
	
saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS') )
	
# print overall results:
#hyper.s1_m1_spike.relAbun.W.spec
## summary (mean values over CV for each tuning step)
#summary(hyper.s1_m1_spike.relAbun.W.spec)
	
sjsdm.s1_m1_spike.relAbun.W.spec = sjSDM(
	Y = as.matrix(otu.data),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale+mean.NDVI.scale +mean.EVI.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)

#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_spike.relAbun.W.spec)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_spike.relAbun.W.spec), sigma = getCov(sjsdm.s1_m1_spike.relAbun.W.spec), history = sjsdm.s1_m1_spike.relAbun.W.spec$history, p = p, logLik=logLik(sjsdm.s1_m1_spike.relAbun.W.spec))
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS')) 
	
```

```{r sjsdm-analyse-correlation-with-mulspec (rel.abun)}
# ... session 1, Malaise I, spike, rel.abun, with spec-data ...
# .. load 'sjSDM_cv' & 'sjSDM' result
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS') )
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS'))
	
#tune = readRDS(here('R','result','s-jSDM_tune_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS') )
#plot(tune, perf='AUC')
#rm(tune)
	
summary(result)
plot(result$history)
	
# . cov <- result$sigma 
dim(result$sigma)
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
otu.table = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
str(otu.table[,1:5])
	
# extract min & max spp pairs of correlation
#source(here("R", "source", "sjsdm-analyse-functions.r"))
extract.minmax.cor(otu=otu.table, sigma=co.env.spp)
	
minmax = my.output[[2]]
head(minmax)
	
#write.table(minmax, file=here::here('correlation_value_s1_m1_spike_relAbun_W.spec.csv'), row.names=F, sep=', ')
rm(minmax, my.output)
	

rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_wSpec','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5),sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	

```


```{r model-analyse-polygon-with-spec (rel.abun)}
# ... session 1, Malaise I, spike, rel.abun, with spec-data ...
# ..... Polygon Drawing .....
# . species association 
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun.W.spec_V0.0.2.9000.RDS'))
	
number=10
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
# variables needed for plotting. can skip. good for checking
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
summary(cuts)
cols = viridis::magma(10) 
# variables needed for plotting. can skip. good for checking
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_wSpec','species_covariance_circular.pdf'), height=8, width=8)
	
version.text = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep='')
otu.text = "sum of relAbun"
	
source(here("R", "source", "sjsdm-analyse-functions.r"))
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl)
	
#dev.off()
	
```

```{r model-analyse-polygon-envir-with-spec (rel.abun)}
# ... session 1, Malaise I, spike, rel.abun, with spec-data ...
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale+mean.NDVI.scale+mean.EVI.scale,  
# 10 variables, excluding intercept

beta = as.matrix(data.frame(result$beta)[,2:11])
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
# variables index & value which has biggest coefficient
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
	
version.text = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep='')
otu.text = "relAbun. spp."
	
evnames=c("ele","canopy","min.T","max.T","preci","road","stream","disturb", 'NDVI', 'EVI')
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun_wSpec','max_environ_spp_cov.pdf'), height=8, width=8)
	
cov.circle.env(sigma=sigma, version.text=version.text, evnames=evnames, otu.text=otu.text)
	
#dev.off()
	
```

```{r sjsdm-analyse-correlation-with-mulspec (0/1)}
# ... session 1, Malaise I, spike, present, with spec-data ...
# .. load 'sjSDM_cv' & 'sjSDM' result
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.present.W.spec.RDS') )
result = readRDS(here('R','result','sjSDM_s1_m1_spike.present.W.spec_V0.0.2.9000.RDS'))
	
summary(result)
plot(result$history)
	
# . cov <- result$sigma 
dim(result$sigma)
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
# extract min & max spp pairs of correlation
spp.names<-colnames(dataI.1.present2[,-c(1:36)])
	
corvalue=data.frame(cor=numeric(), rowname=character(), colname=character())
for (j in 1:dim(sigma)[1]) {
	print(j)
	a = data.frame(cor=sigma[,j], rowname=as.character(seq(1, dim(sigma)[1])), colname=as.character(rep(j,dim(sigma)[1])))
	
	corvalue = rbind(corvalue, a)
	rm(a)
}
	
#write.table(corvalue, file=here::here('correlation_value_s1_m1_spike_present.csv'), row.names=F, sep=', ')
	
dim(corvalue)
corvalue = subset(corvalue, rowname!=colname)
corvalue = corvalue[order(corvalue$cor),]
corvalue.minmax = corvalue[c(1:20,(dim(corvalue)[1]-19):dim(corvalue)[1]),]
corvalue.minmax = corvalue.minmax[seq(1,40,2), ]
	
corvalue.minmax$spp1 = spp.names[corvalue.minmax$rowname]
corvalue.minmax$spp2 = spp.names[corvalue.minmax$colname]
	
#write.table(corvalue.minmax, file=here::here('correlation_value_s1_m1_spike_present_W.spec.csv'), row.names=F, sep=', ')
	

rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_wSpec','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5),sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	

```


```{r model-analyse-polygon-with-spec (0/1)}
# ... session 1, Malaise I, spike, present, with spec-data ...
# ..... Polygon Drawing .....
# . species association .
number=10
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_wSpec','species_covariance_circular.pdf'), height=8, width=8)
	
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -2, y = 6, pos = 3, xpd = NA, labels = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep=''))
#text(x = -6, y = 5.7, pos = 3, xpd = NA, labels = "A", font = 2, cex = 1.5)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,355,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
counter = 1
coords = cbind(xx, yy, angles)
for(i in 2:n) {
	for(j in 1:(i-1)){
      if(to_plot[counter]) {
		print (c(i,j,counter))
		add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, lineSeq = lineSeq)
	}
	counter = counter + 1
  }
}
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
# .. add abundance legends
lineSeq = 5.0
for(i in 1:length(OTU_log)){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = as.character(cuts[i]), lend = 1)
}
	
add_legend(viridis::viridis(11), angles = c(140,110),radius = 5.4)
text(cos(deg2rad(123))*(lineSeq+1), sin(deg2rad(123))*(lineSeq+1.2), labels = "covariance", pos = 2, xpd = NA)
	
add_legend(cols = cols, range = c(2, 850), angles = c(70,40),radius = 5.4)
text(cos(deg2rad(53))*(lineSeq+1), sin(deg2rad(55))*(lineSeq+1.1), labels = "low to high", pos = 4, xpd = NA) 
text(cos(deg2rad(64))*(lineSeq+1.3), sin(deg2rad(62))*(lineSeq+1.1), labels = "sum of presence", pos = 4, xpd = NA) 
	
### arrows
segments(x0 = cos(deg2rad(-1))*(lineSeq-0.2), x1 = cos(deg2rad(-1))*(lineSeq+0.9), y0 = sin(deg2rad(-1))*(lineSeq-0.2), y1 = sin(deg2rad(-1))*(lineSeq+0.9), xpd = NA)
segments(x0 = cos(deg2rad(356))*(lineSeq-0.2), x1 = cos(deg2rad(356))*(lineSeq+0.9), 
         y0 = sin(deg2rad(356))*(lineSeq-0.2), y1 = sin(deg2rad(356))*(lineSeq+0.9), xpd = NA)
	
# first
angles = seq(150,195,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
end = curve_text(195+3, "Species",lineSeq = lineSeq+0.6,reverse = TRUE)
	
# second
angles = seq(rad2deg(end)+3,rad2deg(end)+45+8,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
arrow_angle = max(angles)-2.8
polygon(x = c(cos(deg2rad(arrow_angle))*(lineSeq+0.5), cos(deg2rad(arrow_angle))*(lineSeq+0.7), cos(deg2rad(max(angles)))*(lineSeq+0.6), cos(deg2rad(arrow_angle))*(lineSeq+0.5)),
        y = c(sin(deg2rad(arrow_angle))*(lineSeq+0.5), sin(deg2rad(arrow_angle))*(lineSeq+0.7), sin(deg2rad(max(angles)))*(lineSeq+0.6), sin(deg2rad(arrow_angle))*(lineSeq+0.5)),col = "black", xpd = NA)
	
#dev.off()
	
```


```{r model-analyse-polygon-envir-with-spec (0/1)}
# ... session 1, Malaise I, spike, present, with spec-data ...
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale+mean.NDVI.scale+mean.EVI.scale,  
# 10 variables, excluding intercept

beta = as.matrix(data.frame(result$beta)[,2:11])
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
# variables index & value which has biggest coefficient
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
sigmas = sigma[upper.tri(sigma)]
number=10
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 3.5
nseg = 100
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_present_wSpec','max_environ_spp_cov.pdf'), height=8, width=8)
	
#Drawing figure parameter
par( mar = c(1,2,2.1,2)+0.1)
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -1.7, y = 5.5, pos = 3, xpd = NA, labels = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep=''), font = 2, cex = 1)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,360,length.out = n+1)[1:(n)] # for all otus
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
##inside circle
counter = 1
coords = cbind(xx, yy, angles)
	
for(i in 2:n) {
  for(j in 1:(i-1)){
    
      if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, species = T, lineSeq = 3.5, lwd = 1.3)
      counter = counter + 1
  }
}

lineSeq = 4.0
for(i in 1:n){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = effect_comb$abun[i], lend = 1)
}
lineSeq = 3.5
	
##outside circle 
#formula = ~elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale+mean.NDVI.scale+mean.EVI.scale,  
evnames=c("ele","canopy","min.T","max.T","preci","road","stream","disturb", 'NDVI', 'EVI')
colourCount = length(unique(evnames))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))
cols=getPalette(colourCount)
	
coords = data.frame(cbind(xx, yy, angles))
effect_comb=effect_comb[,-c(3,4)]
effect_comb2 = effect_comb
effect_comb2[,2] = ff(effect_comb[,2])
effect_comb2 = cbind(effect_comb2, effect_comb[,2])
effect_comb2 = data.frame(effect_comb2)
	
for(i in sort(unique(max_effects))) {
  sub<- coords %>% filter(effect_comb2$max_effects==i)
  sub_eff <- effect_comb2 %>% filter(max_effects==i)
  from <- sub[1,3]
  to <- sub[nrow(sub),3]

  x = c((3.6+1.5*(sub_eff[,2]))*cos(deg2rad(sub[,3]) ), 
        rev((3.6+1.5/2)*cos(deg2rad(sub[,3]))))
  
  y = c((3.6+1.5*(sub_eff[,2]))*sin(deg2rad(sub[,3])),
        rev((3.6+1.5/2)*sin(deg2rad(sub[,3]))))
  
  angleName = (from+to)/2
  if(angleName > 180) {reverse = TRUE} else {reverse = FALSE}
  ###environment variable text
  curve_text(angleName, label = evnames[i],reverse = reverse,lineSeq = 5.5, middle = TRUE, extend = 1.1, col = cols[i])
  ###environment variable bar
 if(i == 8) polygon(x-0.1, y, xpd = NA,col = cols[i])
#  else if(i == 5) polygon(x, y+0.55, xpd = NA,col = cols[i])
#  else if(i == 7) polygon(x-0.6, y-0.2, xpd = NA,col = cols[i])
#  else if(i == 9) polygon(x, y-0.55, xpd = NA,col = cols[i])
  #else if(i == 12) polygon(x-0.2, y-0.5, xpd = NA,col = cols[i])
 else 
polygon(x, y, xpd = NA,col = cols[i])
  
  ###environment variable range number
  text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[1,3]+4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[1,3]+4)), 
         xpd = NA, labels = round(min(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
  
   text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[nrow(sub),3]-4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[nrow(sub),3]-4)), 
         xpd = NA, labels = round(max(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
}
	
###legend of bar
rec_cols = viridis::viridis(11)
x = seq(3,5, length.out = 12)
for(i in 1:length(rec_cols)){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = rec_cols[i], xpd = NA, border = NA)
}
text(x[1],-5.2, labels = -1)
text(x[11],-5.2, labels = +1)
	
#abun=as.character(abun)
abun1 = as.character(viridis::magma(10))
x = seq(-5.5,-3, length.out = 11)
for(i in 1:unique(length(abun))){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = abun1[i], xpd = NA, border = NA)
  text(x= x[1]-0.2, y=-5.2, labels = "2", pos = 4, xpd = NA)
  text(x= x[10]-0.2, y=-5.2, labels = '850', pos = 4, xpd = NA)
}
text(x=-5.3, y=-5.3, labels = "presence. spp.", pos = 4, xpd = NA)
	
#dev.off()
	
```

```{r sjsdm-model (rel.abun & 0/1)}
# models are already saved. go directly to 'model-analyse'

# according to nmds, precipitation, elevation, yrs.disturb.min, old.growth.str, T, canopy.height can be used for now

# ... session 1, Malaise I, spike, rel.abun ...
# . sjSDM-version = ‘0.0.2.9000’
names(dataI.1.rel.abun2[,c(1:3,24:39)])
# [1] "SiteName"              "UTM_E"                 "UTM_N"                
# [4] "trap"                  "session"               "site_trap_period"     
# [7] "elevation.scale"       "canopy.ht.scale"       "min.T.scale"          
#[10] "max.T.scale"           "precipitation.scale"   "metre.road.scale"     
#[13] "metre.stream.scale"    "yrs.disturb.min.scale" "mean.NDVI.scale"      
#[16] "mean.EVI.scale"
scale.env = dataI.1.rel.abun2[,c(1:3,24:36)]
str(scale.env)
	
scale.env1 = dataI.1.rel.abun2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale')]
str(scale.env1)
	
otu.data = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
hyper.s1_m1_spike.relAbun = sjSDM_cv(
	Y = as.matrix(otu.data),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 20L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
 saveRDS(hyper.s1_m1_spike.relAbun, file = here('R','result','s-jSDM_tune_s1_m1_spike.relAbun_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_spike.relAbun, perf = "AUC")
	
 saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun_V0.0.2.9000.RDS') )
	
# print overall results:
#hyper.s1_m1_spike.relAbun
## summary (mean values over CV for each tuning step)
#summary(hyper.s1_m1_spike.relAbun)
	
sjsdm.s1_m1_spike.relAbun = sjSDM(
	Y = as.matrix(otu.data),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)

#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_spike.relAbun)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_spike.relAbun), sigma = getCov(sjsdm.s1_m1_spike.relAbun), history = sjsdm.s1_m1_spike.relAbun$history, p = p, logLik=logLik(sjsdm.s1_m1_spike.relAbun))
 saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.relAbun_V0.0.2.9000.RDS')) 
	

# ... session 1, Malaise I, spike, 0/1 ...
# . sjSDM-version = ‘0.0.2.9000’
scale.env1 = dataI.1.present2[,c('elevation.scale','canopy.ht.scale','min.T.scale','max.T.scale','precipitation.scale','metre.road.scale','metre.stream.scale','yrs.disturb.min.scale')]
str(scale.env1)
	
hyper.s1_m1_spike.present = sjSDM_cv(
	Y = as.matrix(as.data.frame(dplyr::select(dataI.1.present2, contains('__')))),
	env = as.matrix(scale.env1),
	biotic = bioticStruct(on_diag=FALSE), tune='random', n_cores = NULL,
	CV=3L, tune_steps = 25L, iter = 60L 
	# iter -> optimization step; lambda -> regularization strength (multiplied constant)
	)
# saveRDS(hyper.s1_m1_spike.present, file = here('R','result','s-jSDM_tune_s1_m1_spike.present_V0.0.2.9000.RDS') )
	
# visualize tuning and best points:
best = plot(hyper.s1_m1_spike.present, perf = "AUC")
	
# print overall results:
hyper.s1_m1_spike.present
# summary (mean values over CV for each tuning step)
summary(hyper.s1_m1_spike.present)
	
sjsdm.s1_m1_spike.present = sjSDM(
	Y = as.matrix(as.data.frame(dplyr::select(dataI.1.present2, contains('__')))),
	iter = 60L, step_size = 27L, link = "probit",
	env = linear(data = as.matrix(scale.env1), formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,
	lambda = best[["lambda_coef"]], alpha = best[["alpha_coef"]]
	), 
	biotic = bioticStruct(lambda = best[["lambda_cov"]], alpha = best[["alpha_cov"]], on_diag = FALSE
	)
)
logLik = logLik(sjsdm.s1_m1_spike.present)
#summary(model)
# calculate post-hoc p-values:
p = getSe(sjsdm.s1_m1_spike.present)
#summary(p)
#summary.p=summary(p)
#plot(sjsdm.s1_m1_spike.relAbun$history)
#save result
result = list(beta = coef(sjsdm.s1_m1_spike.present), sigma = getCov(sjsdm.s1_m1_spike.present), history = sjsdm.s1_m1_spike.present$history, p = p, logLik=logLik)
# saveRDS(result, file = here('R','result','sjSDM_s1_m1_spike.present_V0.0.2.9000.RDS')) 
	

```


```{r sjsdm-model (old version), eval=FALSE, include=FALSE}
#DY use eval=FALSE, include=FALSE
# old sjsdm version, even syntax of model changed!
# just for backup! don't bother to run!

# ... session 1, Malaise I, spike, rel.abun ...
lrs = seq(-18, -1, length.out = 7);f = function(x) 2^x;lrs = f(lrs)
result = vector("list", 7)
	
scale.env = dataI.1.rel.abun2[,c(1:3,24:36)]
str(scale.env)
	
for(i in 1:7) { 
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.rel.abun2[,-c(1:36)]),
               formula = ~ elevation.scale+min.T.scale+max.T.scale+precipitation.scale+mean.NDVI.scale+mean.EVI.scale,  
               learning_rate = 0.01, 
               iter = 100L,
               step_size = 27L, l1_coefs = 0.07*lrs[i], l2_coefs = 0.03*lrs[i], # should be normally tuned...
               sampling = 100L,l1_cov = lrs[i], l2_cov = lrs[i])
  loss=unlist(model$logLik)
  history=model$history
  weights = list(beta = coef(model), sigma = getCov(model),loss=loss,history=history)
  rm(model)
  result[[i]] = weights
 }
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
  saveRDS(result, file = here('R','result','s-jSDM_result_s1_m1_spike.relAbun_lidar.RDS') )
	
```

```{r sjsdm-analyse-correlation (0/1)} 
# ... session 1, Malaise I, spike, present ...
# .. load 'sjSDM_cv' result
#hyper = readRDS(here('R','result','s-jSDM_tune_s1_m1_spike.present.RDS'))
	
#summary(hyper)  # mean values over CV for each tuning step
	
## visualize tuning and best points:
#best = plot(hyper, perf = "AUC")
	
#rm(hyper)
	
#saveRDS(best, file = here('R','result','s-jSDM_tune_best_s1_m1_spike.present.RDS'))
	
# sjSDM_cv result is too big to upload to git, saved only alpha & lambda coef,cov !!!
	
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.present.RDS'))
	
# . load 'sjSDM' result
result = readRDS(here('R','result','sjSDM_s1_m1_spike.present_V0.0.2.9000.RDS'))
	
summary(result)
plot(result$history)
	
# . cov <- result$sigma 
dim(result$sigma)
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
# extract min & max spp pairs of correlation
spp.names<-colnames(dataI.1.present2[,-c(1:36)])
	
corvalue=data.frame(cor=numeric(), rowname=character(), colname=character())
for (j in 1:dim(sigma)[1]) {
	print(j)
	a = data.frame(cor=sigma[,j], rowname=as.character(seq(1, dim(sigma)[1])), colname=as.character(rep(j,dim(sigma)[1])))
	
	corvalue = rbind(corvalue, a)
	rm(a)
}
	
#write.table(corvalue, file=here::here('correlation_value_s1_m1_spike_present.csv'), row.names=F, sep=', ')
	
dim(corvalue)
corvalue = subset(corvalue, rowname!=colname)
corvalue = corvalue[order(corvalue$cor),]
corvalue.minmax = corvalue[c(1:20,(dim(corvalue)[1]-19):dim(corvalue)[1]),]
corvalue.minmax = corvalue.minmax[seq(1,40,2), ]
	
corvalue.minmax$spp1 = spp.names[corvalue.minmax$rowname]
corvalue.minmax$spp2 = spp.names[corvalue.minmax$colname]
	
#write.table(corvalue.minmax, file=here::here('correlation_value_s1_m1_spike_present.csv'), row.names=F, sep=', ')
	


	
rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_present','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5),sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	
```

```{r model-analyse-polygon (0/1)}
# ... session 1, Malaise I, spike, present ...
# ..... Polygon Drawing .....
# . species association .
number=10
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.present2, contains('__')))
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_present','species_covariance_circular.pdf'), height=8, width=8)
	
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -2, y = 6, pos = 3, xpd = NA, labels = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep=''))
#text(x = -6, y = 5.7, pos = 3, xpd = NA, labels = "A", font = 2, cex = 1.5)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,355,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
counter = 1
coords = cbind(xx, yy, angles)
for(i in 2:n) {
	for(j in 1:(i-1)){
      if(to_plot[counter]) {
		print (c(i,j,counter))
		add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, lineSeq = lineSeq)
	}
	counter = counter + 1
  }
}
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
# .. add abundance legends
lineSeq = 5.0
for(i in 1:length(OTU_log)){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = as.character(cuts[i]), lend = 1)
}
	
add_legend(viridis::viridis(11), angles = c(140,110),radius = 5.4)
text(cos(deg2rad(123))*(lineSeq+1), sin(deg2rad(123))*(lineSeq+1.2), labels = "covariance", pos = 2, xpd = NA)
	
add_legend(cols = cols, range = c(2, 850), angles = c(70,40),radius = 5.4)
text(cos(deg2rad(53))*(lineSeq+1), sin(deg2rad(55))*(lineSeq+1.1), labels = "low to high", pos = 4, xpd = NA) 
text(cos(deg2rad(64))*(lineSeq+1.3), sin(deg2rad(62))*(lineSeq+1.1), labels = "sum of presence", pos = 4, xpd = NA) 
	
### arrows
segments(x0 = cos(deg2rad(-1))*(lineSeq-0.2), x1 = cos(deg2rad(-1))*(lineSeq+0.9), y0 = sin(deg2rad(-1))*(lineSeq-0.2), y1 = sin(deg2rad(-1))*(lineSeq+0.9), xpd = NA)
segments(x0 = cos(deg2rad(356))*(lineSeq-0.2), x1 = cos(deg2rad(356))*(lineSeq+0.9), 
         y0 = sin(deg2rad(356))*(lineSeq-0.2), y1 = sin(deg2rad(356))*(lineSeq+0.9), xpd = NA)
	
# first
angles = seq(150,195,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
end = curve_text(195+3, "Species",lineSeq = lineSeq+0.6,reverse = TRUE)
	
# second
angles = seq(rad2deg(end)+3,rad2deg(end)+45+8,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
arrow_angle = max(angles)-2.8
polygon(x = c(cos(deg2rad(arrow_angle))*(lineSeq+0.5), cos(deg2rad(arrow_angle))*(lineSeq+0.7), cos(deg2rad(max(angles)))*(lineSeq+0.6), cos(deg2rad(arrow_angle))*(lineSeq+0.5)),
        y = c(sin(deg2rad(arrow_angle))*(lineSeq+0.5), sin(deg2rad(arrow_angle))*(lineSeq+0.7), sin(deg2rad(max(angles)))*(lineSeq+0.6), sin(deg2rad(arrow_angle))*(lineSeq+0.5)),col = "black", xpd = NA)
	
#dev.off()
	
```


```{r model-analyse-polygon-envir (0/1)}
# ... session 1, Malaise I, spike, present ...
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
# 8 variables, excluding intercept

beta = as.matrix(data.frame(result$beta)[,2:9])
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
# variables index & value which has biggest coefficient
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
sigmas = sigma[upper.tri(sigma)]
number=10
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 3.5
nseg = 100
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_present','max_environ_spp_cov.pdf'), height=8, width=8)
	
#Drawing figure parameter
par( mar = c(1,2,2.1,2)+0.1)
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -1.7, y = 5.5, pos = 3, xpd = NA, labels = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep=''), font = 2, cex = 1)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,360,length.out = n+1)[1:(n)] # for all otus
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
##inside circle
counter = 1
coords = cbind(xx, yy, angles)
	
for(i in 2:n) {
  for(j in 1:(i-1)){
    
      if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, species = T, lineSeq = 3.5, lwd = 1.3)
      counter = counter + 1
  }
}

lineSeq = 4.0
for(i in 1:n){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = effect_comb$abun[i], lend = 1)
}
lineSeq = 3.5
	
##outside circle 
#formula = ~elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
evnames=c("ele","canopy","min.T","max.T","preci","road","stream","disturb")
colourCount = length(unique(evnames))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
cols=getPalette(colourCount)
	
coords = data.frame(cbind(xx, yy, angles))
effect_comb=effect_comb[,-c(3,4)]
effect_comb2 = effect_comb
effect_comb2[,2] = ff(effect_comb[,2])
effect_comb2 = cbind(effect_comb2, effect_comb[,2])
effect_comb2 = data.frame(effect_comb2)
	
for(i in sort(unique(max_effects))) {
  sub<- coords %>% filter(effect_comb2$max_effects==i)
  sub_eff <- effect_comb2 %>% filter(max_effects==i)
  from <- sub[1,3]
  to <- sub[nrow(sub),3]

  x = c((3.6+1.5*(sub_eff[,2]))*cos(deg2rad(sub[,3]) ), 
        rev((3.6+1.5/2)*cos(deg2rad(sub[,3]))))
  
  y = c((3.6+1.5*(sub_eff[,2]))*sin(deg2rad(sub[,3])),
        rev((3.6+1.5/2)*sin(deg2rad(sub[,3]))))
  
  angleName = (from+to)/2
  if(angleName > 180) {reverse = TRUE} else {reverse = FALSE}
  ###environment variable text
  curve_text(angleName, label = evnames[i],reverse = reverse,lineSeq = 5.5, middle = TRUE, extend = 1.1, col = cols[i])
  ###environment variable bar
 if(i == 8) polygon(x-0.1, y, xpd = NA,col = cols[i])
#  else if(i == 5) polygon(x, y+0.55, xpd = NA,col = cols[i])
#  else if(i == 7) polygon(x-0.6, y-0.2, xpd = NA,col = cols[i])
#  else if(i == 9) polygon(x, y-0.55, xpd = NA,col = cols[i])
  #else if(i == 12) polygon(x-0.2, y-0.5, xpd = NA,col = cols[i])
 else 
polygon(x, y, xpd = NA,col = cols[i])
  
  ###environment variable range number
  text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[1,3]+4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[1,3]+4)), 
         xpd = NA, labels = round(min(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
  
   text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[nrow(sub),3]-4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[nrow(sub),3]-4)), 
         xpd = NA, labels = round(max(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
}
	
###legend of bar
rec_cols = viridis::viridis(11)
x = seq(3,5, length.out = 12)
for(i in 1:length(rec_cols)){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = rec_cols[i], xpd = NA, border = NA)
}
text(x[1],-5.2, labels = -1)
text(x[11],-5.2, labels = +1)
	
#abun=as.character(abun)
abun1 = as.character(viridis::magma(10))
x = seq(-5.5,-3, length.out = 11)
for(i in 1:unique(length(abun))){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = abun1[i], xpd = NA, border = NA)
  text(x= x[1]-0.2, y=-5.2, labels = "2", pos = 4, xpd = NA)
  text(x= x[10]-0.2, y=-5.2, labels = '850', pos = 4, xpd = NA)
}
text(x=-5.3, y=-5.3, labels = "presence. spp.", pos = 4, xpd = NA)
	
#dev.off()
	
```


```{r sjsdm-analyse-correlation (rel.abun)} 
# ... session 1, Malaise I, spike, rel.abun ...
# .. load 'sjSDM_cv', 'sjSDM' result
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun_V0.0.2.9000.RDS') )
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun_V0.0.2.9000.RDS')) 
	
summary(result)
plot(result$history)
	
#tune = readRDS(here('R','result','s-jSDM_tune_s1_m1_spike.relAbun_V0.0.2.9000.RDS') )
#plot(tune, perf='AUC')
#rm(tune)
	
# . cov <- result$sigma 
dim(result$sigma)
# [1] 850 850
co.env.spp <- cov2cor(result$sigma)
	
otu.table = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
str(otu.table[,1:5])
	
# extract min & max spp pairs of correlation
extract.minmax.cor(otu=otu.table, sigma=co.env.spp)
	
minmax = my.output[[2]]
table(minmax$spp1)
table(minmax$spp2)
	
#write.table(minmax, file=here::here('correlation_value_s1_m1_spike_relAbun.csv'), row.names=F, sep=', ')
rm(minmax, my.output)
	
rownames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
colnames(co.env.spp) <- 1:dim(result$sigma)[1]   #spp.names
	
#ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
cut.t = cut(co.env.spp, breaks = seq(-1,1,length.out = 12))
summary(cut.t)
	
rm(cut.t)
	
# species correlation plot
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','species_correlation.pdf'), height=20, width=20)
	
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T, title=paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', best[['lambda_coef']],sep='')) 
# > packageVersion('sjSDM')
#[1] ‘0.0.2.9000’
	
#dev.off()
	

```

```{r model-analyse-polygon (rel.abun)}
# ... session 1, Malaise I, spike, rel.abun, no Spec-data ...
# ..... Polygon Drawing .....
# . species association .
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
number=10
	
sigma = re_scale(result$sigma)[order(apply(otu.tbl, 2, sum)), order(apply(otu.tbl, 2, sum))]
	
# variables needed for plotting. can skip. good for checking
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result$sigma)
lineSeq = 4.7
nseg = 100
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
summary(cuts)
cols = viridis::magma(10) 
# variables needed for plotting. can skip. good for checking
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','species_covariance_circular.pdf'), height=8, width=8)
	
version.text = paste('sjSDM version: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']], 5),sep='')
otu.text = "sum of relAbun"
	
cov.circle(version.text=version.text, otu.text=otu.text, sigma=sigma, otu.tbl=otu.tbl)
	
#dev.off()
	
```

```{r model-analyse-polygon-envir (rel.abun)}
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
# 8 variables, excluding intercept
best = readRDS(here('R','result','s-jSDM_tune_best_s1_m1_spike.relAbun_V0.0.2.9000.RDS'))
result = readRDS(here('R','result','sjSDM_s1_m1_spike.relAbun_V0.0.2.9000.RDS'))
	
otu.tbl = as.data.frame(dplyr::select(dataI.1.rel.abun2, contains('__')))
	
beta = as.matrix(data.frame(result$beta)[,2:9])
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result$sigma)# number of otus
max_effects= apply(beta ,1, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[i,max_effects[i]] )))
# variables index & value which has biggest coefficient
	
OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
summary(cuts)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(otu.tbl)[2]), sum = apply(otu.tbl, 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
	
version.text = paste('sjSDM: ',packageVersion('sjSDM'),', kelpie20200214, alpha: ', best[["alpha_coef"]],', lambda: ', round(best[['lambda_coef']],5), sep='')
otu.text = "relAbun. spp."
	
evnames=c("ele","canopy","min.T","max.T","preci","road","stream","disturb")
	
#pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','max_environ_spp_cov.pdf'), height=8, width=8)
	
cov.circle.env(sigma=sigma, version.text=version.text, evnames=evnames, otu.text=otu.text)
	
#dev.off()
	

```


```{r sjsdm-analyse-correlation (old version), eval=FALSE, include=FALSE}
# model result from old sjsdm version
# don't bother to run! just for backup!
# s1, m1, spiked data 
# formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
result<-readRDS(here('R','result','s-jSDM_result_s1_m1_spike.relAbun.RDS'))
str(result)
	
lrs = seq(-18, -1, length.out = 7);f = function(x) 2^x;lrs = f(lrs)

  
# ... look at the loss
result[[1]]$loss
	
loss.all<-data.frame(LogLik = sapply(result, function(r) r$loss[1]))
loss.all$lrs<-lrsdataI.1.rel.abun
str(loss.all)
	
loss.2 = data.frame(sapply(result, function(r) r$loss))
loss.2 = unlist(loss.2)
loss.2<-data.frame(LogLik = loss.2[seq(1,length(loss.2),by=2)], regulation=loss.2[seq(2,length(loss.2),by=2)])
data.frame(loss.2,loss.all$lrs)
	
pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','LogLikelihood.pdf'), height=5, width=5)
	
plot(x=loss.all$lrs,y=loss.all$LogLik,ylab = "LogLik",xlab = "lrs")
	
dev.off()
	

# . cov <- result$sigma 
dim(result[[5]]$sigma)
# [1] 850 850
co.env.spp <- cov2cor(result[[4]]$sigma)
spp.names<-colnames(dataI.1.rel.abun2[,-c(1:36)])
rownames(co.env.spp)<-spp.names
colnames(co.env.spp)<-spp.names
	
rownames(co.env.spp) <- 1:850
colnames(co.env.spp) <- 1:850
	
ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
ggcorrplot(co.env.spp, hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)
	
overall = apply(abind::abind(lapply(result, function(r) cov2cor(r$sigma)), along = 0L), 2:3, sum)
	
# sum of all models with diff lrs
pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','7lrs_spp_corr.pdf'), height=20, width=20)
	
ggcorrplot(cov2cor(overall), hc.order = TRUE, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,title=paste('lrs ',formatC(lrs[1],format='e',3),' to ',lrs[7],', sum of 7 models',sep=''))  
	
dev.off()
	
# .... heatmap
co.env.spp4<-cov2cor(result[[4]]$sigma)
range(co.env.spp4)
co.env.spp1<-cov2cor(result[[1]]$sigma)
co.env.spp7<-cov2cor(result[[7]]$sigma)
	
pdf(here('R','graph','sjsdm_s1_m1_spike_relAbun','spp_corr_heatmap_3lrs.pdf'), height=5, width=15)
	
par(mfrow=c(1,3))
cols = (colorRampPalette(c("blue", "white", "red")))(10)
graphics::image(co.env.spp1[indices$rowInd, indices$rowInd], col = cols, main=paste('lrs ',formatC(lrs[1],format='e',3),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp1),4),' - ', max(co.env.spp1),sep=''), ylab='blue (min value) to red (max value)')
image(co.env.spp4[indices$rowInd, indices$rowInd], col = cols, main=paste('species correlation (850 spp), lrs ',formatC(lrs[4],format='e',3),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp4),4),' - ', max(co.env.spp4),sep=''))
image(co.env.spp7[indices$rowInd, indices$rowInd], col = cols, main=paste('lrs ',round(lrs[7],8),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp7),4),' - ', max(co.env.spp7),sep=''))
	
dev.off()
	

```

```{r model-analyse-polygon (old version), eval=FALSE, include=FALSE}
# model result from old sjsdm version
# don't bother to run! just for backup!
# ..... Polygon Drawing .....
# . species association .
# result[1][4][7], lrs 
number=10
lr_step=7 # change here for different result list
sigma = re_scale(result[[lr_step]]$sigma)[order(apply(dataI.1.rel.abun2[,-c(1:36)], 2, sum)), order(apply(dataI.1.rel.abun2[,-c(1:36)], 2, sum))]

sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result[[lr_step]]$sigma)
lineSeq = 4.7
nseg = 100
	
pdf(paste(here('R','graph','sjsdm_s1_m1_spike_relAbun','spp_cov_lrs'),formatC(lrs[lr_step],format=NULL,3),'.pdf',sep=''), height=8, width=8)
	
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = 0, y = 5.7, pos = 3, xpd = NA, labels = paste("Penalty: ",formatC(lrs[lr_step],format=NULL,4),sep=''))
text(x = -6, y = 5.7, pos = 3, xpd = NA, labels = "A", font = 2, cex = 1.5)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,355,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
counter = 1
coords = cbind(xx, yy, angles)
for(i in 2:n) {
	for(j in 1:(i-1)){
      if(to_plot[counter]) {
		print (c(i,j,counter))
		add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, lineSeq = lineSeq)
	}
	counter = counter + 1
  }
}
# ??? correlation plotted, why legend says 'covariance' 
	
OTU_log = log(sort(apply(dataI.1.rel.abun2[, -c(1:36)], 2, sum)))
range(OTU_log)
OTU_log[1]=0 
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) 
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(dataI.1.rel.abun2[,-c(1:36)])[2]), sum = apply(dataI.1.rel.abun2[,-c(1:36)], 2, sum))
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
OTU_sort_abun$abun<-abun
OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
# .. add abundance legends
lineSeq = 5.0
for(i in 1:length(OTU_log)){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = as.character(cuts[i]), lend = 1)
}
	
add_legend(viridis::viridis(11), angles = c(140,110),radius = 5.4)
text(cos(deg2rad(123))*(lineSeq+1), sin(deg2rad(123))*(lineSeq+1.2), labels = "covariance", pos = 2, xpd = NA)
	
add_legend(cols = cols, range = c(2, 850), angles = c(70,40),radius = 5.4)
text(cos(deg2rad(53))*(lineSeq+1), sin(deg2rad(55))*(lineSeq+1.1), labels = "low to high", pos = 4, xpd = NA) 
text(cos(deg2rad(64))*(lineSeq+1.3), sin(deg2rad(62))*(lineSeq+1.1), labels = "sum of rel. abun.", pos = 4, xpd = NA) 
	
### arrows
segments(x0 = cos(deg2rad(-1))*(lineSeq-0.2), x1 = cos(deg2rad(-1))*(lineSeq+0.9), y0 = sin(deg2rad(-1))*(lineSeq-0.2), y1 = sin(deg2rad(-1))*(lineSeq+0.9), xpd = NA)
segments(x0 = cos(deg2rad(356))*(lineSeq-0.2), x1 = cos(deg2rad(356))*(lineSeq+0.9), 
         y0 = sin(deg2rad(356))*(lineSeq-0.2), y1 = sin(deg2rad(356))*(lineSeq+0.9), xpd = NA)
	
# first
angles = seq(150,195,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
end = curve_text(195+3, "Species",lineSeq = lineSeq+0.6,reverse = TRUE)
	
# second
angles = seq(rad2deg(end)+3,rad2deg(end)+45+8,length.out = n+1)[1:(n)]
xx = cos(deg2rad(angles))*(lineSeq+0.6)
yy = sin(deg2rad(angles))*(lineSeq+0.6)
lines(xx, yy, xpd = NA)
arrow_angle = max(angles)-2.8
polygon(x = c(cos(deg2rad(arrow_angle))*(lineSeq+0.5), cos(deg2rad(arrow_angle))*(lineSeq+0.7), cos(deg2rad(max(angles)))*(lineSeq+0.6), cos(deg2rad(arrow_angle))*(lineSeq+0.5)),
        y = c(sin(deg2rad(arrow_angle))*(lineSeq+0.5), sin(deg2rad(arrow_angle))*(lineSeq+0.7), sin(deg2rad(max(angles)))*(lineSeq+0.6), sin(deg2rad(arrow_angle))*(lineSeq+0.5)),col = "black", xpd = NA)
	
dev.off()
	
```

```{r model-analyse-polygon-envir (old version), eval=FALSE, include=FALSE}
# model result from old sjsdm version
# don't bother to run! just for backup!
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
# result[1][4][7]
lr_step = 7 # change here for different result list
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
# 8 variables, excluding intercept
beta = as.matrix(result[[lr_step]]$beta[[1]])[2:9,]
effects= apply(beta, 1, function(o) sum(abs(o)))
	
turn_over = 1
n = ncol(result[[lr_step]]$sigma)# number of otus
max_effects= apply(beta ,2, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) beta[max_effects[i],i] )))
# variables index & value which has biggest coefficient
	
sppname3=seq(1:20); sppname3=as.character(sppname3)
effect_comb$name <- rep(sppname3, len=850)
effect_comb$abun <- OTU_sort_abun$abun
effect_comb$abun<-as.character(effect_comb$abun)
effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
effect_comb = effect_comb[effect_comb_ind,]
	
sigma = re_scale(result[[lr_step]]$sigma)[effect_comb_ind, effect_comb_ind]
sigmas = sigma[upper.tri(sigma)]
number=10
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result[[lr_step]]$sigma)
lineSeq = 3.5
nseg = 100
	
pdf(paste(here('R','graph','sjsdm_s1_m1_spike_relAbun','max_environ_lrs'),formatC(lrs[lr_step],format=NULL,3),'.pdf',sep=''), height=8, width=8)
	
#Drawing figure parameter
par( mar = c(1,2,2.1,2)+0.1)
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -3.5, y = 5.3, pos = 3, xpd = NA, labels = paste("Penalty: ",formatC(lrs[lr_step],format=NULL,3), ', loss: ', formatC(loss.all[lr_step,1],7),sep=''), font = 2, cex = 1)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,360,length.out = n+1)[1:(n)] # for all otus
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
##inside circle
counter = 1
coords = cbind(xx, yy, angles)
	
for(i in 2:n) {
  for(j in 1:(i-1)){
    
      if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, species = T, lineSeq = 3.5, lwd = 1.3)
      counter = counter + 1
  }
}

lineSeq = 4.0
for(i in 1:n){
  p1 = coords[i,]
  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = effect_comb$abun[i], lend = 1)
}
lineSeq = 3.5
	
##outside circle 
#formula = ~elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
evnames=c("ele","canopy","min.T","max.T","preci","road","stream","disturb")
colourCount = length(unique(evnames))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
cols=getPalette(colourCount)

coords = data.frame(cbind(xx, yy, angles))
effect_comb=effect_comb[,-c(3,4)]
effect_comb2 = effect_comb
effect_comb2[,2] = ff(effect_comb[,2])
effect_comb2 = cbind(effect_comb2, effect_comb[,2])
effect_comb2 = data.frame(effect_comb2)
	
for(i in sort(unique(max_effects))) {
  sub<- coords %>% filter(effect_comb2$max_effects==i)
  sub_eff <- effect_comb2 %>% filter(max_effects==i)
  from <- sub[1,3]
  to <- sub[nrow(sub),3]

  x = c((3.6+1.5*(sub_eff[,2]))*cos(deg2rad(sub[,3]) ), 
        rev((3.6+1.5/2)*cos(deg2rad(sub[,3]))))
  
  y = c((3.6+1.5*(sub_eff[,2]))*sin(deg2rad(sub[,3])),
        rev((3.6+1.5/2)*sin(deg2rad(sub[,3]))))
  
  angleName = (from+to)/2
  if(angleName > 180) {reverse = TRUE} else {reverse = FALSE}
  ###environment variable text
  curve_text(angleName, label = evnames[i],reverse = reverse,lineSeq = 5.5, middle = TRUE, extend = 1.1, col = cols[i])
  ###environment variable bar
 if(i == 8) polygon(x-0.1, y, xpd = NA,col = cols[i])
#  else if(i == 5) polygon(x, y+0.55, xpd = NA,col = cols[i])
#  else if(i == 7) polygon(x-0.6, y-0.2, xpd = NA,col = cols[i])
#  else if(i == 9) polygon(x, y-0.55, xpd = NA,col = cols[i])
  #else if(i == 12) polygon(x-0.2, y-0.5, xpd = NA,col = cols[i])
 else 
polygon(x, y, xpd = NA,col = cols[i])
  
  
  ###environment variable range number
  text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[1,3]+4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[1,3]+4)), 
         xpd = NA, labels = round(min(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
  
   text(srt = 0, 
         x = (3.6+1.5)*cos(deg2rad(sub[nrow(sub),3]-4)), 
         y =  (3.6+1.5)*sin(deg2rad(sub[nrow(sub),3]-4)), 
         xpd = NA, labels = round(max(sub_eff[,3]), 2), col = cols[i], cex = 0.8)
}
	
###legend of bar
rec_cols = viridis::viridis(11)
x = seq(3,5, length.out = 12)
for(i in 1:length(rec_cols)){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = rec_cols[i], xpd = NA, border = NA)
}
text(x[1],-5.2, labels = -1)
text(x[11],-5.2, labels = +1)
	
#abun=as.character(abun)
abun1 = as.character(viridis::magma(10))
x = seq(-5.5,-3, length.out = 11)
for(i in 1:unique(length(abun))){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = abun1[i], xpd = NA, border = NA)
  text(x= x[1]-0.2, y=-5.2, labels = "2", pos = 4, xpd = NA)
  text(x= x[10]-0.2, y=-5.2, labels = '850', pos = 4, xpd = NA)
}
text(x=-5.3, y=-5.3, labels = "rel. abun. spp.", pos = 4, xpd = NA)
	
dev.off()
	

```



