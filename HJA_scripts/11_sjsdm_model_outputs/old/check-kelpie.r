# Yuanheng
# modified: Dec 24, 2019
# last modified: Feb 12, 2020
# check data (kelpie20191217), make clean dataset & apply PCA, iNext

rm(list=ls())
setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon/kelpie")
getwd()
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse'),library,character.only=T)
	
# ...................... read-in data ..........................
otu1 = read.csv('sample_by_species_table_F2308_minimap2_20191219_kelpie20191217.csv', header=T, sep=',')
otu2 = read.csv('sample_by_species_table_F2308_f0x2_minimap2_20191219_kelpie20191217.csv', header=T, sep=',')
# 'f0x2' different setting of 'minimap2'
	
print(c(dim(otu1), dim(otu2)))
	
# ............ correct & extract sampling info .................
# (M - malaise, S - session) 
otu1$site_trap_period
	
otu2$site_trap_period
	
site_info = read.csv('shotgun_samples_HJAndrews_Malaisetraps_20200108.csv', header=T, sep='\t')
str(site_info)
	
site_info$site
	
sort(table(site_info$sample))
hist(table(site_info$sample))
	
unique(site_info$sample)
	 
# HOBO-032-M1-S2-1 HOBO-032-M1-S2 HOBO-357-M1-S2

site_info$Correct_Site[site_info$sample=='HOBO-357-M1-S2'| site_info$sample=='HOBO-351-M1-S1'] 
site_info$Correct_Period[site_info$sample=='HOBO-357-M1-S2'| site_info$sample=='HOBO-351-M1-S1'] 
site_info$Correct_Trap[site_info$sample=='HOBO-032-M1-S2-1'| site_info$sample=='HOBO-032-M1-S2']
	
unique(site_info$Correct_Site)
hist(table(site_info$Correct_Site))
	

# ... set U1/2 to 1/2 ...
site_info$correct_trap_change = gsub('U1','1', site_info$Correct_Trap)
site_info$correct_trap_change = gsub('U2','2', site_info$correct_trap_change)
	


# ........ add correct site info to otu file .........
as.character(sort(unique(site_info$sample)))[1:5]
as.character(sort(unique(otu1$site_trap_period)))[1:5]
	
# . site info in 'otu1' & 'otu2' is same 
unique(otu1$site_trap_period) == unique(otu2$site_trap_period)
	
# . match site name of 'otu1' & 'site_info' .
otu1$sample = gsub('_','-', otu1$site_trap_period)
	
otu1$sample[1:5]
	
otu2$sample = gsub('_','-', otu2$site_trap_period)
	
as.character(sort(unique(site_info$sample))) == as.character(sort(unique(otu2
$sample)))
	
# 50, 68-70
as.character(sort(unique(site_info$sample)))[68:70]
as.character(sort(unique(otu2$sample)))[68:70]
	
as.character(sort(unique(otu1$sample))) == as.character(sort(unique(otu2$sample)))
	
# otu1 & out2 have identical '$sample'
site_compare = data.frame(name=as.character(rep('aaa',length=241)), otu.index=rep(0,length=241), site.index=rep(0,length=241), stringsAsFactors=F)
	
j=1
for (i in 1:length(unique(otu1$sample))) {
	a = as.character(sort(unique(otu1$sample))[i])
	
	if(match(a, unique(site_info$sample),nomatch=999) != 999) {
		site_compare$name[j]=a
		site_compare$otu.index[j] = match(a, sort(unique(otu1$sample)))
		site_compare$site.index[j] = match(a, sort(unique(site_info$sample)))
		j=j+1
	} else {
		print(a)}
	
}
	
# . there's 5 sites without corrected info from 'site_info' ......
#[1] "207746-M1-S1"
#[1] "268355-M2-S1"
#[1] "HOBO-060-M1-S1"
#[1] "HOBO-351-M1-S2"
#[1] "HOBO-353-M1-S2"
	
hist(table(site_compare$site.index))
	
which(site_compare$otu.index!=site_compare$site.index)
# ... 19 sites
# [1]  67  68 151 152 153 154 155 210 211 212 213 214 215 216 217 218 219 220 221
	

site_compare = subset(site_compare, name!= 'aaa')
	

otu1$correct_site = rep('Znot_sign', dim(otu1)[1])
otu1$correct_session = rep(0, dim(otu1)[1])
otu1$correct_trap = rep(0, dim(otu1)[1])
	
otu2$correct_site = rep('Znot_sign', dim(otu2)[1])
otu2$correct_session = rep(0, dim(otu2)[1])
otu2$correct_trap = rep(0, dim(otu2)[1])
	
for (i in 1:length(site_compare$name)) {
	k = site_compare$name[i]
	a = match(k, otu1$sample)
	b = match(k, otu2$sample)
	
	if (length(a)>1) print(c(k, a))
	if (length(b)>1) print(c(k, b))
	
	site = as.character(site_info$Correct_Site[match(k, site_info$sample)])
	session = site_info$Correct_Period[match(k, site_info$sample)]
	trap = site_info$correct_trap_change[match(k, site_info$sample)]
	
	otu1$correct_site[a] = site
	otu1$correct_session[a] = session
	otu1$correct_trap[a] = trap
	
	otu2$correct_site[b] = site
	otu2$correct_session[b] = session
	otu2$correct_trap[b] = trap
	
	print(i)
	rm(a,b,site,session,trap)
}
	
table(otu1$correct_site)
# 5 sites not signed
	
table(otu2$correct_site)
# 5 sites not signed
	
# ... replace not clear site info to NA ...
# . otu1 .
#HOBO-351-M1-S2 or HOBO-357-M1-S2
otu1$correct_site[otu1$correct_site=='Znot_sign' | otu1$correct_site=='HOBO-351-M1-S2 or HOBO-357-M1-S2'] = 'NULL'
	
otu1$correct_session[otu1$correct_site=='NULL'] = NA
otu1$correct_trap[otu1$correct_site=='NULL'] = NA
	
# . otu2 .
otu2$correct_site[otu2$correct_site=='Znot_sign' | otu2$correct_site=='HOBO-351-M1-S2 or HOBO-357-M1-S2'] = 'NULL'
	
otu2$correct_session[otu2$correct_site=='NULL'] = NA
otu2$correct_trap[otu2$correct_site=='NULL'] = NA
	
# .. save the edited otu file ..
write.table(otu1, 'correct_site-sample_by_species_table_F2308.csv', row.names=F, sep=',')
	
write.table(otu2, 'correct_site-sample_by_species_table_F2308_f0x2.csv', row.names=F, sep=',')
	

# HOBO-357_M1_S2 562701, 4895822 (otu2)
# 			     562701, 4895822 (otu1)
check1 = subset(otu1, otu1$site_trap_period=='HOBO-357_M1_S2')
check2 = subset(otu2, otu2$site_trap_period=='HOBO-357_M1_S2')
	
check1[1,888:890]
check2[1,888:889]
	
# ............... simple calculation ................
# -> species richness
spp.rich1 = otu1
spp.rich1[ ,4:893] = (otu1[ ,4:893]>0)*1
	
table(spp.rich1[,16])
spp.rich1[16,1:8]
	
rich1 = spp.rich1[,-(1:3)]
spp.rich1$count = rowSums(rich1)
rm(rich1)
	
spp.rich2 = otu2
spp.rich2[ ,4:889] = (otu2[ ,4:889]>0)*1
	
table(spp.rich2[,116])
spp.rich2[116,1:8]
	
rich2 = spp.rich2[,-(1:3)]
spp.rich2$count = rowSums(rich2)
rm(rich2)
	
g1 = ggplot(spp.rich1, aes(UTM_E, UTM_N)) + geom_point( aes(size=count))+  
 scale_radius(name = waiver(), breaks = waiver(), labels = waiver(),
  limits = NULL, range = c(1, 5), trans = "identity",
  guide = "legend") + ggtitle('minimap2 (default)') + theme(plot.title= element_text(hjust=.5))
	
spp.rich1$UTM_N[101:121]
#spp.rich1$UTM_E[101:121]
# there's a NA in coordination (row #112)
	
g2 = ggplot(spp.rich2, aes(UTM_E, UTM_N)) + geom_point(aes(size=count))+  
 scale_radius(name = waiver(), breaks = waiver(), labels = waiver(),
  limits = NULL, range = c(1, 4), trans = "identity",
  guide = "legend") + ggtitle('minimap2 (f0x2)') + theme(plot.title= element_text(hjust=.5))
	
pdf('R/graph/spp.rich.2data.pdf', width=6, height=7)
	
grid.arrange(g1,g2, nrow=2)
	
dev.off()
	
spp.rich2$UTM_N[101:121]
spp.rich2$UTM_E[101:121]
# there's a NA in coordination (row #112)
	
spp.rich1$count[112]
spp.rich2$count[112]
# only no coordinate info for row #112
	
# -> copy numbers
spp.copy1 = otu1[,1:3]
head(spp.copy1)
	
rich1 = otu1[,-(1:3)]
spp.copy1$cp.num = rowSums(rich1)
rm(rich1)
	
spp.copy2 = otu2[,1:3]
head(spp.copy2)
	
rich2 = otu2[,-(1:3)]
spp.copy2$cp.num = rowSums(rich2)
rm(rich2)
	
g1 = ggplot(spp.copy1, aes(UTM_E, UTM_N)) + geom_point( aes(size=cp.num))+  
 scale_radius(name = waiver(), breaks = waiver(), labels = waiver(),
  limits = NULL, range = c(1, 6), trans = "identity",
  guide = "legend") + ggtitle('minimap2 (default)') + theme(plot.title= element_text(hjust=.5))
	
g2 = ggplot(spp.copy2, aes(UTM_E, UTM_N)) + geom_point(aes(size=cp.num))+  
 scale_radius(name = waiver(), breaks = waiver(), labels = waiver(),
  limits = NULL, range = c(1, 3.2), trans = "identity",
  guide = "legend") + ggtitle('minimap2 (f0x2)') + theme(plot.title= element_text(hjust=.5))
	
pdf('R/graph/copy.number.2data.pdf', width=7, height=8)
	
grid.arrange(g1,g2, nrow=2)
	
dev.off()
	
# ... coordination
site.coord = spp.rich1[,1:3]
head(site.coord)
	
# ..................... transformation for PCA ....................
# Hellinger pre-transformation
decostand(data, 'hellinger')
	
# Chi-square transformation

	
# .......................... explore models ..............................
# NMSC

count.nmds = metaMDS(otu2[,-c(1:3)], distance='bray', trymax=2000, k=2)
count.nmds
	
count.nmds2 = metaMDS(otu2[,-c(1:3)], distance='bray', trymax=2000, k=2, previous.best=count.nmds)
	
count.nmds3 = metaMDS(otu2[,-c(1:3)], distance='bray', trymax=2000, k=2, previous.best=count.nmds2)
	

bestnmds()
	
wascores()
	
par(mfrow=c(1,2))
	
stressplot(count.nmds, main='Shepard')
	
gof = goodness(count.nmds)
plot(count.nmds, type='p', main='goodness of fit')
points(count.nmds, display='sites',cex=gof*2)
	



# initial configuration
# PCoA (principle coordinate analysis)
# Q mode -> relationship among objects; R mode -> among variables
?vegdist
	
head(names(spp.rich1))
otu1[1:5,1:10]
ls()
	
names(spp.copy1)
	
count.bray = vegdist(otu1[,-c(1:3)]) 
count.bray.pcoa = cmdscale(count.bray, k=dim(otu1)[1]-4, eig=T)
	
count.bray.pcoa$species = wascores(count.bray.pcoa$points[,1:2], otu1[,-c(1:3)])
	
# .. top 10 prevalent
names(spp.rich1)[1:5]
 table(spp.rich1$count)
	
spp.rich1$site_trap_period[spp.rich1$count>=81]
which(spp.rich1$count>=81)
	
spp.rich1$site_trap_period[spp.rich1$count<=11]
which(spp.rich1$count<=11)
	
?ordiplot
	
str(count.bray.pcoa)
head(count.bray.pcoa$points)
	
pdf('R/graph/PCoA.1.pdf', width=7, height=8)
#par(mfrow=c(1,2))
	
pl = ordiplot(count.bray.pcoa, type='none', main='PCoA (Bray-Curtis dissimilarity)')
points(pl, 'sites', pch=8)
	
abline(h=0, lty=3)
abline(v=0, lty=3)
	
# ... mark 7 sites with sparsest spp (richness/abundant) BLUE
# . richness
#points(count.bray.pcoa$points[which(spp.rich1$count<=11),1:2], pch=8,col= 'blue', cex=1.5)
	
table(spp.copy1$cp.num)
	
# . abundant
points(count.bray.pcoa$points[which(spp.copy1$cp.num<=4177),1:2], pch=8,col= 'blue', cex=1.5)
	
# ... mark 7 sites with richest spp (richness/abundant) RED
# . richness
#points(count.bray.pcoa$points[which(spp.rich1$count>=81),1:2], pch=8,col= 'red', cex=1.5)
	
# . abundant
points(count.bray.pcoa$points[which(spp.copy1$cp.num>=24869),1:2], pch=8,col= 'purple', cex=1.5)
	
# ... add species
points(pl, 'species', col='red', pch=20)
	
# ... mark most prevalent/abundant spp
spp.prevalent = data.frame(preva = colSums(spp.rich1[,-c(1,2,3,dim(spp.rich1)[2])]), cp.num = colSums(otu1[,-c(1,2,3)]))
row.names(spp.prevalent)
	
table(spp.prevalent$preva)
	
str(count.bray.pcoa)
	
# . prevalent
row.names(spp.prevalent)[which(spp.prevalent$preva>=94)]
count.bray.pcoa$species[which(spp.prevalent$preva>=94),]
	
points(count.bray.pcoa$species[which(spp.prevalent$preva>=94),1], count.bray.pcoa$species[which(spp.prevalent$preva>=94),2], pch=16, col='lightgreen')
	
# . abundant
table(spp.prevalent$cp.num)
	
row.names(spp.prevalent)[which(spp.prevalent$cp.num>=31837)]
count.bray.pcoa$species[which(spp.prevalent$cp.num>=400000),]
	
points(count.bray.pcoa$species[which(spp.prevalent$cp.num>=31837),1], count.bray.pcoa$species[which(spp.prevalent$cp.num>=31837),2], pch=1, col='green')
	
legend('bottomright', pch=c(8,8,1),col=c('purple','blue','green'),legend=c('sites with most cp.num', 'sites with least cp.num', 'most abundant spp'), bty='n')
	
dev.off()
	



ordiplot(scores(count.bray.pcoa)[,1:2], type='points')
	
text(count.bray.pcoa.spp[1:5,], rownames(count.bray.pcoa.spp)[1:5], col='red')
	




















