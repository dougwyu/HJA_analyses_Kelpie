# Yuanheng
# create: Jul 20, 2020
# data from 'explore_nmds_lidar_envir.rmd'
# scale lidar variables


```{r set-up}
setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon/R-git")
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse','scatterplot3d', 'gridBase','grid', 'here'),library,character.only=T)
  
here::here()
# .../R-git !!!
getwd()
	
```

```{r read-in data}
otu.env1.spike.present2 = read.table(here('kelpie','formatted_data','present_lidar_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), sep=',', header=T, na.strings='NA')
	
otu.env1.spike.relAbun2 = read.table(here('kelpie','formatted_data','relAbun_lidar_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), sep=',', header=T, na.strings='NA')
	
# can skip if directly goes to 'subsets of data'
# 0/1 OTU
otu.env1.spike.present = read.csv(here('kelpie','formatted_data','present2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(otu.env1.spike.present[,1:26])
	
# relAbun OTU
otu.env1.spike.relAbun = read.csv(here('kelpie','formatted_data','relAbun2_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(otu.env1.spike.relAbun[,1:26])
	
```


```{r scale-lidar-variables}
dim(otu.env1.spike.present2)
dim(otu.env1.spike.relAbun2)
	
names(otu.env1.spike.present2)[1:50] == names(otu.env1.spike.present2)[1:50]
table(names(otu.env1.spike.present2)[51:1197] == names(otu.env1.spike.present2)[51:1197])
	
a = data.frame(select(otu.env1.spike.present2,40:50)%>%rename(l_Cover_2m_4m.scale=1, l_Cover_2m_4m_all.scale=2, l_Cover_2m_max.scale=3, l_Cover_2m_max_all.scale=4, l_Cover_4m_16m.scale=5, l_Cover_4m_16m_all.scale=6, l_p25.scale=7, l_p25_all.scale=8, l_p95.scale=9, l_p95_all.scale=10, l_rumple.scale=11)%>% scale())
str(a)
	
par(mfrow=c(1,2))
hist(a$l_Cover_2m_4m.scale)
hist(otu.env1.spike.present2$l_Cover_2m_4m)
	
otu.env1.spike.present2.scale = cbind(otu.env1.spike.present2[,1:39], a, otu.env1.spike.present2[,51:dim(otu.env1.spike.present2)[2]])
dim(otu.env1.spike.present2.scale)
dim(otu.env1.spike.present2)
	
names(otu.env1.spike.present2.scale)[26:51]
	
otu.env1.spike.relAbun2.scale = dplyr::left_join(otu.env1.spike.relAbun2, select(otu.env1.spike.present2.scale, 'site_trap_period',starts_with('l_')), by=c('site_trap_period', 'site_trap_period'), copy=F)
dim(otu.env1.spike.relAbun2.scale)
names(otu.env1.spike.relAbun2.scale)[1197:1208]
	
otu.env1.spike.relAbun2.scale = otu.env1.spike.relAbun2.scale[,c(1:39,1198:1208,51:1197)]
str(otu.env1.spike.relAbun2.scale[,1:50])
	
# write.table(otu.env1.spike.present2.scale, here('kelpie','formatted_data','present2_lidar_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	
# write.table(otu.env1.spike.relAbun2.scale, here('kelpie','formatted_data','relAbun2_lidar_mulspec_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv'), row.names=F, sep=',')
	

