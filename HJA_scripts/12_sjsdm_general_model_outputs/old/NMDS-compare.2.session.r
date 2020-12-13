# Yuanheng
# create: Mar 09, 2020
# last modified:  , 2020
# apply NMDS for 2 sessions with sites that appear in both sessions
# with data (in'data_Feb_25/', 'kelpie20200214')


# .................................................................
rm(list=ls())
setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon/kelpie")
getwd()
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse','scatterplot3d', 'gridBase','grid'),library,character.only=T)
	
# ...................... read-in data ..........................
# ('data_Feb_25' folder) 
otu.env1.noS = read.csv('data_Feb_25/sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
otu.env1.spike = read.csv('data_Feb_25/sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
# 1173-26, more otus
	
names(otu.env1.noS)[1:26] == names(otu.env1.spike)[1:26]
	
names(otu.env1.noS)[1:26]=c('site', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
names(otu.env1.spike)[1:26]=c('site', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
	
str(otu.env1.noS[,1:26])
	
otu.env1.noS$hja = as.character(otu.env1.noS$hja)
otu.env1.noS$yrs.disturb.level = as.character(otu.env1.noS$yrs.disturb.level)
otu.env1.noS$poly.ID = as.character(otu.env1.noS$poly.ID) 
otu.env1.noS$point.ID = as.character(otu.env1.noS$point.ID)
	
otu.env1.spike$hja = as.character(otu.env1.spike$hja)
otu.env1.spike$yrs.disturb.level = as.character(otu.env1.spike$yrs.disturb.level)
otu.env1.spike$poly.ID = as.character(otu.env1.spike$poly.ID) 
otu.env1.spike$point.ID = as.character(otu.env1.spike$point.ID)
	
# .. subsets of data (spike or not) ..
# no spike
table(otu.env1.noS$session)
	
dataI.1.1.noS = subset(otu.env1.noS, session == 'S1' & trap == 'M1' )
dataI.1.2.noS = subset(otu.env1.noS, session == 'S1' & trap == 'M2' )
	
dataI.2.1.noS = subset(otu.env1.noS, session == 'S2' & trap == 'M1' )
dataI.2.2.noS = subset(otu.env1.noS, session == 'S2' & trap == 'M2' )
	
print(c(dim(dataI.1.1.noS),dim(dataI.1.2.noS), dim(dataI.2.1.noS), dim(dataI.2.2.noS)))
# 88+33+91+25=237
	
# with spike
table(otu.env1.spike$session)
	
dataI.1.1.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M1' )
dataI.1.2.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M2' )
	
dataI.2.1.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M1' )
dataI.2.2.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M2' )
	
print(c(dim(dataI.1.1.spike),dim(dataI.1.2.spike), dim(dataI.2.1.spike), dim(dataI.2.2.spike)))
# 88+33+91+25=237
	
# ..... find sites appeared in both session .....
# ... (Malaise 1) ...
names(otu.env1.noS)[1:26]
sort(unique(otu.env1.noS$site)) == sort(unique(otu.env1.spike$site))
	
par(mfrow=c(1,2))
plot(table(dataI.1.1.noS$site))
plot(table(dataI.2.1.noS$site))
	
match.s12 = data.frame(index.s1 = rep(0,length=dim(dataI.2.1.spike)[1]), index.s2 = rep(0,length=dim(dataI.2.1.spike)[1]), site = rep('NULL',length=dim(dataI.2.1.spike)[1]))
match.s12$site = as.character(match.s12$site)
match.s12
	
for (i in 1:dim(dataI.1.1.noS)[1]) {
	a = match(dataI.1.1.noS$site[i], dataI.2.1.noS$site, nomatch=999)
	if (a==999) {
		print(i)
	}
	match.s12$site[i] = dataI.1.1.noS$site[i]
	match.s12$index.s1[i] = i
	match.s12$index.s2[i] = a
}
	
# ... data with sites in both session ...
# . (malaise 1) .
print(c(dim(dataI.1.1.noS), dim(dataI.2.1.noS)))
	
data.m1.s1.noS= dataI.1.1.noS[match.s12$index.s1[match.s12$index.s2!=999], ]
	
data.m1.s2.noS = dataI.2.1.noS[match.s12$index.s2[match.s12$index.s2!=999], ]
	
print(c(dim(data.m1.s1.noS), dim(data.m1.s2.noS)))
	
# ...............................................................
# ................. NMDS with env.variable no spike......................
# . (88 sites) no spike .
dist.m1.s1.noS = vegdist(wisconsin(data.m1.s1.noS[,-(1:26)]))
str(dist.m1.s1.noS)
#86*85/2=3655
	
set.seed(756)
count.nmds.m1.s1.noS = metaMDS(dist.m1.s1.noS, distance='bray', try=30, trymax=2000, k=3, trace=T,plot=T, binary=F)
#  Procrustes: rmse 0.001713051  max resid 0.009046282
	
dist.m1.s2.noS = vegdist(wisconsin(data.m1.s2.noS[,-(1:26)]))
str(dist.m1.s2.noS)
#86*85/2=3655
	
set.seed(756)
count.nmds.m1.s2.noS = metaMDS(dist.m1.s2.noS, distance='bray', try=30, trymax=2000, k=3, trace=T,plot=T, binary=F)
#  Procrustes: rmse 0.001240137  max resid 0.004605783
	
count.nmds.m1.s1.noS
count.nmds.m1.s2.noS
	
pdf('R/graph/noSpike_m1_s12_nmds.pdf', width=12, height=16)
	
par(mfrow=c(3,2), cex=1.1)
	
stressplot(count.nmds.m1.s1.noS, main=paste('kelpie20200214, M1 (86 sites), S1, stress ', round(count.nmds.m1.s1.noS$stress,4), sep=''))
stressplot(count.nmds.m1.s2.noS, main=paste('Malaise 1 (86 sites), session 2, stress ', round(count.nmds.m1.s2.noS$stress,4), sep=''))
	
plot(count.nmds.m1.s1.noS, choice=c(1,2), main = 'session 1, Dim 1&2')
plot(count.nmds.m1.s1.noS, choice=c(1,3), main = 'session 1, Dim 1&3')
	
plot(count.nmds.m1.s2.noS, choice=c(1,2), main = 'session 2, Dim 1&2')
plot(count.nmds.m1.s2.noS, choice=c(1,3), main = 'session 2, Dim 1&3')
	
dev.off()
	
# . procrustes test .
proc1i = procrustes(count.nmds.m1.s1.noS, count.nmds.m1.s2.noS, symmetric = T)
	
summary(proc1i)
residuals(proc1i)
fitted(proc1i)
	
	
test = protest(count.nmds.m1.s1.noS, count.nmds.m1.s2.noS)
summary(test)
test$ss # signif
	
 

	
pdf('R/graph/noSpike_m1_s12_protest.pdf', height=10,width=10)
	
par(mfrow=c(2,2),tex=1.1)
	
plot(proc1i, choices=c(1,2), main='Dim 1&2') #,type = "text"
	
plot(proc1i, choices=c(1,3), main='kelpie20200214, M1, S1 (86 sites), Dim 1&2')
	
plot(proc1i, choices=c(2,3), main='Dim 2&3')
#points(proc1i, choices=c(2,3))
	
plot(test, kind=2, main=paste('sum of squares ', round(test$ss,4), ', significance ', test$signif, sep=''))
	
dev.off()
	










