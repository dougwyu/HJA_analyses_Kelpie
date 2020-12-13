# Yuanheng
# create: Mar 23, 2020
# last modified:  Apr , 2020
# apply sj-sdm & hmsc 
# with data (in'data_Feb_25/', 'kelpie20200214')


# revised working directory code
library(here)
here() # this reads working directory as the pathname of the 9_ecological_analyses.Rproj file. It therefore varies transparently from user to user. You can then use 'here' as a drop-in for file.path(), as i show below
# [1] "/Users/Negorashi2011/Dropbox/Working_docs/Luo_Mingjie_Oregon/HJA_analyses_Kelpie/HJA_scripts/9_ecological_analyses"

# .................................................................
# rm(list=ls())
# setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon")
# getwd()
	
source('R/source/sjsdm_function.r')
# missing from repo
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse','scatterplot3d', 'gridBase','grid', 'ggcorrplot', 'here'),library,character.only=T)
	
lapply(c('Hmsc','sjSDM', 'reticulate'),library, character.only=T)
	
# ...................... kelpie, lidar data ..............................
# (lidar data)
# lidar.env = read.csv('kelpie/lidar_data/biodiversity_site_info_multispectral_2020-04-13.txt', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
lidar.env = read.csv(here("..", "10_eo_data", "biodiversity_site_info_multispectral_2020-04-13.txt"), header=T, sep=',', stringsAsFactors = F, na.strings='NA') # this formulation works on any computer with RStudio

str(lidar.env)
	
# ('data_Feb_25' folder) 
# otu.env1.noS = read.csv('kelpie/data_Feb_25/sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
otu.env1.noS = read.csv(here("..", "..", "Kelpie_maps", "outputs_minimap2_20200221_F2308_f0x2_q48_kelpie20200214_vsearch97", "sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv"), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
# otu.env1.spike = read.csv('kelpie/data_Feb_25/sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
otu.env1.spike = read.csv(here("..", "..", "Kelpie_maps", "outputs_minimap2_20200221_F2308_f0x2_q48_kelpie20200214_vsearch97", "sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv"), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
# 1173-26, more otus
	
names(otu.env1.noS)[1:26] == names(otu.env1.spike)[1:26]
	
names(otu.env1.noS)[1:26]=c('SiteName', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
names(otu.env1.spike)[1:26]=names(otu.env1.noS)[1:26]
	
str(otu.env1.spike[,1:27])
	
# ........ formate lidar.env .......
sort(lidar.env$SiteName) == sort(unique(otu.env1.spike$SiteName))
	
#NDVI - normalized difference vegetation index (calculated using bands 4 and 5): NDVI = (NearIR-Red)/(NearIR+Red)
#       these values should range between -1 and 1. Values in these columns should be divided by 1000
#EVI - enhanced vegetation index (calculated using bands 4, 5, and 2):  2.5 * ((Band 5 – Band 4) / (Band 5 + 6 * Band 4 – 7.5 * Band 2 + 1))
#      the values in these columns should be divided by 1000
names(lidar.env)   #sort()
# 4 NDVI, 4 EVI
data.frame(lidar.env$nor.NDVI_20180717, lidar.env$EVI_20180818)
	
a = lidar.env %>% select(13,14,27,28,41,42,55,56) %>% rename(nor.NDVI_20180717=1, nor.EVI_20180717=2, nor.NDVI_20180726=3, nor.EVI_20180726=4, nor.NDVI_20180802=5, nor.EVI_20180802=6, nor.NDVI_20180818=7, nor.EVI_20180818=8)/1000
	
a[,c(1,3,5,7)]
# mean of all 4 NDVI, EVI
lidar.env= cbind(lidar.env, a)
str(lidar.env)
	
# pdf('kelpie/R/graph/lidar_describe_EVI_NDVI.pdf', height=5, width=5)
pdf(here("..", "..", "sjSDM", "lidar_describe_EVI_NDVI.pdf"), height=5, width=5)
	
range(lidar.env[,60:67])
plot(1:96, lidar.env[,61], ylim=c(0,3.3),type='l', main='solid - EVI, dash - NDVI')
lines(1:96, lidar.env[,63], col='red')
lines(1:96, lidar.env[,65], col='blue')
lines(1:96, lidar.env[,67], col='green')
	
lines(1:96, lidar.env[,60], lty=2, col='black')
lines(1:96, lidar.env[,62], lty=2, col='red')
lines(1:96, lidar.env[,64], lty=2, col='blue')
lines(1:96, lidar.env[,66], lty=2, col='green')
	
dev.off()
	
lidar.env$mean.NDVI = base::rowMeans(lidar.env[,c(60,62,64,66)])
lidar.env$mean.EVI = base::rowMeans(lidar.env[,c(61,63,65,67)])
lidar.env[,c(60,62,64,66,68)]
lidar.env[,c(61,63,65,67,69)]
	
# ... explore OTU reads ...
names(otu.env1.spike)[1:27]
	
hist(otu.env1.spike[,30])
sort(unique(otu.env1.noS[,30]))
	
which(is.na(otu.env1.noS[,28]))#:dim(otu.env1.noS)[2]
which(is.na(otu.env1.noS[,29]))
# 181 237+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "s-jSDM_result_s1_m1_spike2.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	
	
table(is.na(otu.env1.noS[c(181,23	
for(i in 1:7) { 
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.spike2[,-c(1:36)]),
               formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "R/result/s-jSDM_result_s1_m1_spike2.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	7),27:1173]))
	
table(is.na(otu.env1.spike[c(181,237),27:1173]))
	
otu.env1.noS$SiteName[c(181,237)]
	
# delete row 181, 237 -> "HOBO-040", "290415" 
dim(otu.env1.noS)
	
otu.env1.noS = otu.env1.noS[-c(181,237),]
	
dim(otu.env1.spike)
	
otu.env1.spike = otu.env1.spike[-c(181,237),]
	
# ....... scale variables .......
a = select(otu.env1.spike,15:22)%>%rename(elevation.scale=1,canopy.ht.scale=2,min.T.scale=3, max.T.scale=4, precipitation.scale=5, metre.road.scale=6, metre.stream.scale=7, yrs.disturb.min.scale=8)%>% scale()
		
for(i in 1:7) { 
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.spike2[,-c(1:36)]),
               formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "R/result/s-jSDM_result_s1_m1_spike2.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	
otu.env1.spike =+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "s-jSDM_result_s1_m1_spike2.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	 cbind(otu.env1.spike[,1:26], data.frame(a), otu.env1.spike[,27:dim(otu.env1.spike)[2]])
dim(otu.env1.spike)
	
otu.env1.noS = dplyr::left_join(otu.env1.noS, otu.env1.spike[,26:34], by=c('site_trap_period', 'site_trap_period'), copy=F)
otu.env1.noS = otu.env1.noS[,c(1:26,1174:1181,27:1173)]
str(otu.env1.noS[,1:34])
	
# ..... combine lidar .....	
for(i in 1:7) { 
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.spike2[,-c(1:36)]),
               formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "R/result/s-jSDM_result_s1_m1_spike2.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	
str(lidar.env)
	
hist(lidar.env$mean.NDVI)
hist(lidar.env$mean.EVI)
	
a = select(lidar.env,68:69) %>% rename(mean.NDVI.scale=1,mean.EVI.scale=2) %>% scale()
lidar.env = cbind(lidar.env,data.frame(a))
	
hist(lidar.env$mean.NDVI.scale)
hist(lidar.env$mean.EVI.scale)
	
# row 181, 237 -> "HOBO-040", "290415" in OTU datasets
lidar.env$SiteName == sort(unique(otu.env1.noS$SiteName))
	
otu.env1.spike = dplyr::left_join(otu.env1.spike, lidar.env[,c(1,70:71)], by=c('SiteName', 'SiteName'), copy=F) 
str(otu.env1.spike[,1:37])
	
otu.env1.spike = otu.env1.spike[,c(1:34,1182:1183,35:1181)]
dim(otu.env1.spike)
	
otu.env1.noS = dplyr::left_join(otu.env1.noS, lidar.env[,c(1,70:71)], by=c('SiteName', 'SiteName'), copy=F) 
str(otu.env1.noS[,1:37])
	
otu.env1.noS = otu.env1.noS[,c(1:34,1182:1183,35:1181)]
dim(otu.env1.noS)
	
write.table(otu.env1.noS, 'kelpie/data_Feb_25/lidar_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv', row.names=F, sep=',')
	
write.table(otu.env1.spike, 'kelpie/data_Feb_25/lidar_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', row.names=F, sep=',')
	
write.table(lidar.env, 'kelpie/lidar_data/sjsdm_biodiversity_site_info_multispectral_2020-04-13.csv', row.names= F, sep=',')
	
# ......................... subsets of data ..................................
dataI.1.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M1' )
dataI.2.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M2' )
print(c(dim(dataI.1.spike), dim(dataI.2.spike)))
	
dataII.1.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M1' )
dataII.2.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M2' )
print(c(dim(dataII.1.spike), dim(dataII.2.spike)))
	
# .... if there's all zero OTUs ....
a = data.frame(c.name=names(dataI.1.spike), zero = c(rep('env',length=36), rep('F', length= dim(dataI.1.spike)[2]-36)))
a$zero=as.character(a$zero)
for (i in 37:dim(dataI.1.spike)[2]) {
	print(i)
	if (sum(dataI.1.spike[,i])==0) {
		a$zero[i]='T'
		}
}
	
dataI.1.spike2 = dataI.1.spike[,c(1:36,which(a$zero=='F'))]
dataI.1.spike2[1:5,34:40]
	

# according to nmds, precipitation, elevation, yrs.disturb.min, old.growth.str, T, canopy.height can be used for now

# ........................ sjSDM ..............................
# ... session 1, Malaise I, spike ...
lrs = seq(-18, -1, length.out = 7);f = function(x) 2^x;lrs = f(lrs)
result = vector("list", 7)
	
scale.env = dataI.1.spike2[,c(1:3,24:36)]
str(scale.env)
	
for(i in 1:7) { 
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.spike2[,-c(1:36)]),
               formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "R/result/s-jSDM_result_s1_m1_spike2.RDS") 
# "R/result/s-jSDM_result_s1_m1_spike.RDS" -> with all OTUs (1147)
	
# ... session 1, Malaise I, spike, with lidar ...
lrs = seq(-18, -1, length.out = 7);f = function(x) 2^x;lrs = f(lrs)
result = vector("list", 7)
	
scale.env = dataI.1.spike2[,c(1:3,24:36)]
str(scale.env)
	
for(i in 1:7) { #1:30
model = sjSDM(X = scale.env, 
               Y = as.matrix(dataI.1.spike2[,-c(1:36)]),
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
  saveRDS(result, file = "s-jSDM_result_s1_m1_spike_lidar2.RDS") 
# "s-jSDM_result_s1_m1_spike_lidar.RDS" -> with all OTUs
	





	
# ... all data, spiked ...
lrs = seq(-18, -1, length.out = 30);f = function(x) 2^x;lrs = f(lrs)
result = vector("list", 30)
	

 for(i in c(1,5,10,15,20,25,30)) { #1:30
model = sjSDM(X = scale.env, 
               Y = as.matrix(otu.env1.spike[,-c(1:26)]),
               formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
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
  saveRDS(result, file = "s-jSDM_result_otu_env1_spike.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	
# . with Lidar data
str(lidar.env)
	
# deleted row 181, 237 -> "HOBO-040", "290415" in OTU datasets
lidar.env$SiteName == sort(unique(scale.env$SiteName))
	
table(scale.env$SiteName)
	
scale.env = dplyr::left_join(scale.env, lidar.env[,c(1,68,69)], by=c('SiteName', 'SiteName'), copy=F) 
str(scale.env)
	
hist(scale.env$mean.NDVI)
hist(scale.env$mean.EVI)
	
a = select(scale.env,19:20) %>% rename(mean.NDVI.scale=1,mean.EVI.scale=2) %>% scale()
scale.env = cbind(scale.env,data.frame(a))
	
hist(scale.env$mean.NDVI.scale)
hist(scale.env$mean.EVI.scale)
	
result2 = vector("list", 30)
	
 for(i in c(1,5,10,15,20,25,30)) { #1:30
model = sjSDM(X = scale.env, 
               Y = as.matrix(otu.env1.spike[,-c(1:26)]),
               formula = ~ elevation.scale+min.T.scale+max.T.scale+precipitation.scale+mean.NDVI.scale+mean.EVI.scale,  
               learning_rate = 0.01, 
               iter = 100L,
               step_size = 27L, l1_coefs = 0.07*lrs[i], l2_coefs = 0.03*lrs[i], # should be normally tuned...
               sampling = 100L,l1_cov = lrs[i], l2_cov = lrs[i])
  loss=unlist(model$logLik)
  history=model$history
  weights = list(beta = coef(model), sigma = getCov(model),loss=loss,history=history)
  rm(model)
  result2[[i]] = weights
 }
  saveRDS(result2, file = "s-jSDM_result_otu_env1_spike_lidar1.RDS") 
  #iter = 1000L, step_size = 6L, sampling = 100L, learning_rate = 0.0003
	







# ........................ hmsc ..............................
Y.I.1.noS = as.matrix(dataI.1.noS[,-(1:26)])
df.I.1.noS = data.frame(dataI.1.noS[,c(1:6,15:26)]) 
	
md.I.1.noS = Hmsc(Y=Y.I.1.noS, XData = df.I.1.noS, XFormula = ~yrs.disturb.min+precipitation+max.T, distr='lognormal poisson')
	
# fit HMSC model with Bayesian inference -> sampleMcmc
m.I.1.noS = sampleMcmc(md.I.1.noS,
					thin = 5, samples = 1000,
							  # num of samples to obtain per chain
					transient = 500*5, nChains = 2, 
					verbose = 500*5
					)
	
save(m.I.1.noS, file= "m.I.1.noS")
	
m.I.1.noS.2 = sampleMcmc(md.I.1.noS,
					thin = 5, samples = 5000,
							  # num of samples to obtain per chain
					transient = 500*5, nChains = 3, 
					verbose = 500*5
					)
	
save(m.I.1.noS.2, file= "m.I.1.noS.2")
	
m.I.1.noS.3 = sampleMcmc(md.I.1.noS,
					thin = 5, samples = 3000,
							  # num of samples to obtain per chain
					transient = 500*5, nChains = 3, 
					verbose = 500*5
					)
	
save(m.I.1.noS.3, file= "m.I.1.noS.3")
	
mpost = convertToCodaObject(m.normal)
	
load('m.I.1.noS.3')
	
# .... hmsc with coordinates ....
studyDesign = data.frame(sample=as.factor(dataI.1.noS$site))
xycoords = as.matrix(dataI.1.noS[,2:3], ncol=2)
rownames(xycoords) = studyDesign$sample
rL = HmscRandomLevel(sData=xycoords)
	
md.I.1.noS.II.poisson = Hmsc(Y=Y.I.1.noS, XData = df.I.1.noS, XFormula=~yrs.disturb.min+precipitation+max.T, studyDesign=studyDesign, ranLevels=list('sample'=rL), distr='lognormal poisson')
	
md.I.1.noS.II.normal = Hmsc(Y=Y.I.1.noS, XData = df.I.1.noS, XFormula=~yrs.disturb.min+precipitation+max.T, studyDesign=studyDesign, ranLevels=list('sample'=rL))
	
m.I.1.noS.II.n.1 = sampleMcmc(md.I.1.noS.II.normal,
					thin = 5, samples = 1000,
							  # num of samples to obtain per chain
					transient = 500*5, nChains = 3, 
					verbose = 500*5
					)
	
save(m.I.1.noS.II.n.1, file= "m.I.1.noS.II.n.1")
	
m.I.1.noS.II.p.1 = sampleMcmc(md.I.1.noS.II.poisson,
					thin = 5, samples = 1000,
							  # num of samples to obtain per chain
					transient = 500*5, nChains = 3, 
					verbose = 500*5
					)
	
save(m.I.1.noS.II.p.1, file= "m.I.1.noS.II.p.1")
	
# .......... analyze hmsc .............
# normal distribution
load('m.I.1.noS')
	
preds.m.I.1.noS = computePredictedValues(m.I.1.noS)
MF.m.I.1.noS = evaluateModelFit(hM=m.I.1.noS, predY=preds.m.I.1.noS)
	
str(MF.m.I.1.noS)
# RMSE: root-mean-square error
plot(MF.m.I.1.noS$O.AUC, MF.m.I.1.noS$O.TjurR2)
# O.AUC = 0.5, O.TjurR2 = 0, failed i think
	
	
load('m.I.1.noS.3')
	
preds.m.I.1.noS.3 = computePredictedValues(m.I.1.noS.3)
load('m.I.1.noS.2')
	
preds.m.I.1.noS.2 = computePredictedValues(m.I.1.noS.2)
# Error: cannot allocate vector of size 11.3 Gb
	
load('m.I.1.noS.3')
	
preds.m.I.1.noS.3 = computePredictedValues(m.I.1.noS.3)
# Error: cannot allocate vector of size 6.8 Gb
# (probably cuz running hmsc at the same time. memory in my laptop is used up)
MF.m.I.1.noS.3 = evaluateModelFit(hM=m.I.1.noS.3, predY=preds.m.I.1.noS.3)
	
mpost = convertToCodaObject(m.noS.inf)
eff.beta=effectiveSize(mpost$Beta)					# effective sample size
hist(eff.beta)
	
psrf.beta=gelman.diag(mpost$Beta, multivariate=T)$psrf		# potential scale reduction factor




