# Yuanheng
# create: Feb 17, 2020
# last modified:  , 2020
# apply NMDS for different sessions & mark them on geographic map
# with environmental data (in 'data_Feb_17/', 'data_Feb_25/')
# Feb 14, working with correct data
# derived from 
    
# ..............................
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)................................
# ............... explore functions ....................
data(package = .packages(all.available = TRUE))  # to list the data sets in all *available* packages.
data(list=c('varespec', 'varechem'))
str(varespecshape=hja,shape=hja,)
	
load(file='BCI.env.rda')
ls() # to see the name of the loaded file
	

# .................................................................
rm(list=ls())
setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon/kelpie")
getwd()
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse','scatterplot3d', 'gridBase','grid'),library,character.only=T)
	
# ...................... read-in data ..........................
# ('data_Feb_25' folder) 
otu.env1.noS = read.csv('data_Feb_25/sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
str(otu.env1.noS[,1:26])
	
unique(otu.env1.noS$trap)
unique(otu.env1.noS$period)
	
otu.env1.spike = read.csv('data_Feb_25/sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
str(otu.env1.spike[,1:26])
	
new.env1 = read.csv('data_Feb_17/biodiversity_site_info_vars_20200220.csv', header=T, sep=',', stringsAsFactors = F)
str(new.env1)
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
# 1173-26, more otus
	
names(new.env1)[3:24] == names(otu.env1)[2:23]
	
otu1 = otu.env1[,-c(2:23)]
dim(otu1)
	
otu1=dplyr::left_join(otu1, new.env1[,-1], by='SiteName')
dim(otu1)
names(otu1)
	
otu.env1[,2:23] = otu1[,895:916]
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
rm(otu1)
	
names(otu.env1.noS)[1:26] == names(otu.env1.spike)[1:26]
	
names(otu.env1.noS)[1:26]=c('site', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
names(otu.env1.spike)[1:26]=c('site', 'UTM_E','UTM_N','old.growth.str', 'yrs.disturb','point.ID','poly.ID','AGENCY','unit.log', 'name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T', 'precipitation','metre.road', 'metre.stream', 'yrs.disturb.min','hja','trap','session','site_trap_period')
	
otu.env1.noS$hja = as.character(otu.env1.noS$hja)
otu.env1.noS$yrs.disturb.level = as.character(otu.env1.noS$yrs.disturb.level)
otu.env1.noS$poly.ID = as.character(otu.env1.noS$poly.ID) 
otu.env1.noS$point.ID = as.character(otu.env1.noS$point.ID)
	
otu.env1.spike$hja = as.character(otu.env1.spike$hja)
otu.env1.spike$yrs.disturb.level = as.character(otu.env1.spike$yrs.disturb.level)
otu.env1.spike$poly.ID = as.character(otu.env1.spike$poly.ID) 
otu.env1.spike$point.ID = as.character(otu.env1.spike$point.ID)
	
names(otu.env1)[900:916]
# 1-26 site shape=hja,info, 27-916
names(otu.old1)[890:893]
# 4-893
sort(names(otu.env1)[27:916]) == sort(names(otu.old1)[4:893])
	
rm(otu.old1)
	
# num of sites in 'otu.env1' is different from in corrected files

env1 = data.frame(otu.env1[,1:26])
str(env1)
	
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
is.na(env1[,6:7])
	
par(mfrow=c(1,2))
hist(env1$point.ID)
hist(env1$poly.ID)
	
table(env1$AGENCY)
table(env1$UNIT_NAME)
	
table(env1$TREATMENT_)
hist(env1$e)
hist(env1$ht)
	
# 'o'->old-growth structural index; 'y','class'->years since disturbance(satellite) to 2018 ???; 'YrsSinceDi'->yrs since logging(to 2018) why 2018 ???; 'YrsSinceDisshape=hja,t'->combine 'y'&'YrsSinceDi', no info given 200; 'TREATMENT_'->logging treatment; 'e'->elevation; 'ht'->canopy height; 'hja'->not clear!!!
	
# .... re-name variables ....
rename(env1, 'elevation'='e','canopy.ht'=c'ht')
	
a=c('site','UTM_E','UTM_N','old.growth.str','yrs.disturb','point.ID','poly.ID','AGENCY','unit.log','name.log','log.yr','yrs.log.2018','log.treat','yrs.disturb.level', 'elevation','canopy.ht','min.T','max.T','precipitation','metre.road','metre.stream','yrs.disturb.min','hja','trap','session','site_trap_period')
	
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
data.frame(a,names(env1))
	
names(env1) = a
str(env1)
	
names(otu.env1)[1:26]=a
	
# ... save file with new col.names ...
write.table(otu.env1, 'data_Feb_17/sample_by_species_table_F2308_minimap2_20191219_kelpie20191217_tosa20200220.csv', row.names=F,sep=',')
	
# ............ sub-sets of data ............
dim(otu.env1.noS)
unique(otu.env1.noS$session)
table(otu.env1.noS$trap)
	
match(NA,otu.env1.noS$session)
match(NA,otu.env1.noS$session[182:dim(otu.env1.noS)[1]])
otu.env1.noS$trap[181:237]
	
otu.env1.noS$hja[otu.env1.noS$tr
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)ap=='NA']
	
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
	
unique(dataI.1.1.spike$trap)
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
unique(dataI.1.1.noS[,110])
	
# ...............................................................
# ................. NMDS with env.variable with/no spike......................
# . S1, OUT1 (88 sites) no spike .
distI1.noS = vegdist(wisconsin(dataI.1.1.noS[,-(1:26)]))
str(distI1.noS)
#88*87/2=3828
	
set.seed(8756)
count.nmdsI1.noS = metaMDS(distI1.noS, distance='bray', try=30, trymax=2000, k=2, trace=T,plot=T, binary=F)
#  Procrustes: rmse 0.0001431462  max resid 0.0008225903
	
set.seed(18756)  # set the seed so that the result can be repeated
count.nmdsI1_3.noS = metaMDS(distI1.noS, distance='bray', try=30, trymax=2000, k=3, trace=T,plot=T)
#  Procrustes: rmse 0.001676282  max resid 0.00992772 
	
count.nmdsI1.noS
# Stress:   shape=hja,  0.2692499
count.nmdsI1_3.noS
# Stress:     0.2102818
	

# ..... fit env.data .....   
fitI1_3.noS = envfit(count.nmdsI1_3.noS~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.1.1.noS, na.rm=T, choices=c(1:3))
	 
fitI1_3.noS 
	
fitI1_3.noS.2 = envfit(count.nmdsI1_3.noS~yrs.disturb.min +elevation+max.T +precipitation+min.T+ hja+canopy.ht, data=dataI.1.1.noS, na.rm=T, choices=c(1:3))
	
fitI1_3.noS.2
plot(count.nmdsI1_3.noS, choices=c(1,2))
plot(fitI1_3.2, choices=c(1,2)) 
	
a = data.frame(fitI1_3[[1]]$arrows)
range(scores(count.nmdsI1_3.noS))
range(a)
d=scatterplot3d(range(a), range(a), range(a),type="n") ;
#,xlab="Abundance of Species 1", ylab="Abundance of Species 2",zlab="Abundance of Species 3");
d$points3d(scores(count.nmdsI1_3), pch=19)
d$points3d(a, pch=19, col='red')
	shape=hja,
p2 <- d$xyz.convert(a[2,1],a[2,2],a[2,3])
p3 <- d$xyz.convert(0,0,0)
segments(p2$x,p2$y,p3$x,p3$y,lwd=2,col=2)
# ......................................................................
# ................. plotting .....................
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1",
  ylab="Abundance of Species 2",zlab="Abundance of Species 3"); d
d$points3d(5,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(3,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")
d$points3d(0,5,5); text(d$xyz.convert(0,5,5.5),labels="community C")




# ..... combine NMDS to geographic map .....
# ... s1, OTU1 (88 sites) no spike ...
names(dataI.1.1.noS)[1:26]
	
# . make a table for plotting (scores & goodness of fit of nmds & coordinates) .
stressplot(count.nmdsI1_3.noS)
	
gofIs1_3.noS = goodness(count.nmdsI1_3.noS)
range(gofIs1_3.noS)
	
nmds.geo.Is1_3.noS = data.frame(UTM_E = dataI.1.1.noS$UTM_E, UTM_N = dataI.1.1.noS$UTM_N, nmds1 = count.nmdsI1_3.noS$points[,1], nmds2 = count.nmdsI1_3.noS$points[,2],nmds3 =count.nmdsI1_3.noS$points[,3], good=gofIs1_3.noS, ele=dataI.1.1.noS$elevation, disturb=dataI.1.1.noS$yrs.disturb.min,hja=as.character(dataI.1.1.noS$hja),min.T=dataI.1.1.noS$min.T,preci=dataI.1.1.noS$precipitation,c.ht = dataI.1.1.noS$canopy.ht)
str(nmds.geo.Is1_3.noS)
	
c(range(nmds.geo.Is1$nmds1), range(nmds.geo.Is1$nmds2))
	
# . NMDS & hja, disturb .
plot(count.nmdsI1_3.noS, main=, type='n')
points(count.nmdsI1_3.noS, display='sites',cex=gofIs1_3.noS*50, pch=19)
	
pdf('R/graph/nmds_geo_nmds2_env_Is1II_0303_noSpike.pdf', width=18, height=11)
	
par(mfrow=c(2,3))
	
# . I. nmds with env.variables, plot 1:3
plot(count.nmdsI1_3.noS, type='n', pch=19, choices=c(1,2))
points(count.nmdsI1_3.noS, display='sites', pch=19)
plot(fitI1_3.noS, choices=c(1,2)) 
	
plot(count.nmdsI1_3.noS, type='n', pch=19, choices=c(1,3), xlim=c(min(nmds.geo.Is1_3.noS$nmds1)*1.1,max(nmds.geo.Is1_3.noS$nmds1)*1.1), ylim=c(min(nmds.geo.Is1_3.noS$nmds3)*1.3, max(nmds.geo.Is1_3.noS$nmds3)*1.3), main='no spike, session 1, trap 1')
points(count.nmdsI1_3.noS, display='sites', pch=19)
plot(fitI1_3.noS, choices=c(1,3)) 
	
plot(count.nmdsI1_3.noS, type='n', pch=19, choices=c(2,3), xlim=range(nmds.geo.Is1_3.noS$nmds2), ylim=c(min(nmds.geo.Is1_3.noS$nmds3)*1.25, max(nmds.geo.Is1_3.noS$nmds3)*1.25))
points(count.nmdsI1_3.noS, display='sites', pch=19)
plot(fitI1_3.noS, choices=c(3,2)) 
	
# . II. geo plots
plot.new()              ## suggested by @Josh
vps <- baseViewports()
#pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
#vp1 <-plotViewport(c(.1,.1,0,.1))
	
g.Is1.1 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=disturb, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('disturbance') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
g.Is1.3 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, shape=hja, size=c.ht, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('hja & canopy height') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='hja')
#g.Is1.2 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=c.ht, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('canopy height') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) 
	
vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
	
vp.BottomLeft <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.0)
	
#vp.BottomMiddle<- viewport(height=unit(.5, "npc"), width=unit(0.333, "npc"), 
#                           just=c("left","top"), 
#                           y=0.5, x=0.335)
	
print(g.Is1.1, vp=vp.BottomLeft)
#print(g.Is1.2, vp=vp.BottomMiddle)
print(g.Is1.3, vp=vp.BottomRight)
	
dev.off()
	
# . nmds & T, ele, rain .
g.Is1.1 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=ele, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.2 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=ele, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
g.Is1.3 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=preci, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.4 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=preci, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
g.Is1.5 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=min.T, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, (min)temperature') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.6 = ggplot(nmds.geo.Is1_3.noS, aes(UTM_E, UTM_N, size=min.T, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, (min)temperature') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
pdf('R/graph/nmds_geo_nmds_T_ele_preci_Is1II_0225_noSpike.pdf', width=11, height=15)
	
pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(.5,5,5,5), "null"))))
grid.text("kelpie 20200214, session 1, Malaise 1 (88 sites)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontfamily='sans'))
print(g.Is1.1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(g.Is1.2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
	
print(g.Is1.3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(g.Is1.4, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
	
print(g.Is1.5, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(g.Is1.6, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
	
dev.off()
	
# ................. NMDS with env.variable with/no spike......................
# . S2, OUT1 (91 sites) no spike .
dim(dataI.2.1.noS)
	
names(dataI.2.1.noS)[1:27]
	
distI2.noS = vegdist(wisconsin(dataI.2.1.noS[,-(1:26)]))
str(distI2.noS)
# 91*45=4095
	
set.seed(8756)
count.nmdsI2_3.noS = metaMDS(distI2.noS, distance='bray', try=30, trymax=6000, k=3, trace=T,plot=T) # , previous.best=cmdscale(distI2.noS
# Procrustes: rmse 0.001024561  max resid 0.004121664 
	
count.nmdsI2_3.noS
# Stress:     0.2380555 
	
stressplot(count.nmdsI2_3.noS)
	
fitI2_3.noS = envfit(count.nmdsI2_3.noS~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.2.1.noS, na.rm=T, choices=c(1:3))
	 
fitI2_3.noS
	
gofIs2_3.noS = goodness(count.nmdsI2_3.noS)
range(gofIs2_3.noS)
	
set.seed(333)  # set the seed so that the result can be repeated
count.nmdsI2_4.noS = metaMDS(distI2.noS, distance='bray', try=30, trymax=2000, k=4, trace=T,plot=T)
# Procrustes: rmse 0.0007636882  max resid 0.003287006
	
count.nmdsI2_4.noS
# Stress:     0.1907971 
	
# fit environmental variables
fitI2_4.noS = envfit(count.nmdsI2_4.noS~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.2.1.noS, na.rm=T, choices=c(1:4))
	 
fitI2_4.noS
	
stressplot(count.nmdsI2_4.noS)
	
gofIs2_4.noS = goodness(count.nmdsI2_4.noS)
range(gofIs2_4.noS)
	
# compare 4&3 dimension plots
pdf('R/graph/nmds_4_3_env_Is2_noSpike.pdf', width=18, height=18)
	
par(mfrow=c(3,3))
	
# 3 dim
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(1,2))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(1,2)) 
	
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(2,3), main=paste("kelpi 20200214, 3 mds, noSpike, session 2 (M1), stress ",round(count.nmdsI2_3.noS$stress,4), sep=''))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(2,3)) 
	
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(1,3))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(1,3)) 
	
# 4 dim
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(1,2))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(1,2)) 
	
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(2,3), main=paste("4 mds, noSpike, session 2 (M1, sites ",dim(dataI.2.1.noS)[1],"), stress ",round(count.nmdsI2_4.noS$stress,4), sep=''))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(2,3)) 
	
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(1,3))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(1,3)) 
	
# 4 dim
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(1,4))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(1,4)) 
	
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(2,4), main=paste("4 mds, noSpike, session 2 (M1, sites ",dim(dataI.2.1.noS)[1],"), stress ",round(count.nmdsI2_4.noS$stress,4), sep=''))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(2,4)) 
	
plot(count.nmdsI2_4.noS, type='n', pch=19, choices=c(3,4))
points(count.nmdsI2_4.noS, display='sites', pch=19)
plot(fitI2_4.noS, choices=c(3,4)) 
	
dev.off()
	
# ..... combine NMDS to geographic map .....
# ... s2, OTU1 (91 sites) no spike ...
nmds.geo.Is2_3.noS = data.frame(UTM_E = dataI.2.1.noS$UTM_E, UTM_N = dataI.2.1.noS$UTM_N, nmds1 = count.nmdsI2_3.noS$points[,1], nmds2 = count.nmdsI2_3.noS$points[,2],nmds3 =count.nmdsI2_3.noS$points[,3], good=gofIs2_3.noS, ele=dataI.2.1.noS$elevation, disturb=dataI.2.1.noS$yrs.disturb.min,hja=as.character(dataI.2.1.noS$hja),min.T=dataI.2.1.noS$min.T,max.T=dataI.2.1.noS$max.T,preci=dataI.2.1.noS$precipitation,c.ht = dataI.2.1.noS$canopy.ht, road=dataI.2.1.noS$metre.road, old=dataI.2.1.noS$old.growth.str)
str(nmds.geo.Is2_3.noS)
	
# nmds 3 dimensions
g1.dis = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=disturb, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('disturbance, hja') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g1.T = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=max.T, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('max Temperature, canopy') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g1.ele = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=ele, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g3.ele = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=ele, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g1.rain = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=preci, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g3.rain = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=preci, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
	
g2.road = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=road, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('distance to road') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) 
g3.road = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=road, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('distance to road') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) 
	
g2.old = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=old, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('old growth structure') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) 
	
g2.hja = ggplot(nmds.geo.Is2_3.noS, aes(UTM_E, UTM_N, size=c.ht, shape=hja, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('canopy height, hja') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) 
	
pdf('R/graph/nmds_geo_nmds3_12_env_Is2_noSpike.pdf', width=18, height=18)
	
par(mfrow=c(3,3))
	
# . I. nmds with env.variables, plot 1:3
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(1,2))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(1,2)) 
	
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(2,3), main=paste("3 mds, noSpike, session 2 (M1), stress ",round(count.nmdsI2_3.noS$stress,4), sep=''))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(2,3)) 
	
plot(count.nmdsI2_3.noS, type='n', pch=19, choices=c(1,3))
points(count.nmdsI2_3.noS, display='sites', pch=19)
plot(fitI2_3.noS, choices=c(1,3)) 
	
# . II. geo plots
plot.new()
vps <- baseViewports()
	
vp.MRight <- viewport(height=unit(.33, "npc"), width=unit(0.5, "npc"), just=c("left","top"), y=0.67, x=0.5)
	
vp.MLeft <- viewport(height=unit(.33, "npc"), width=unit(0.5, "npc"), just=c("left","top"), y=0.67, x=0.0)
	
vp.BLeft<- viewport(height=unit(.33, "npc"), width=unit(0.5, "npc"), just=c("left","top"), y=0.333, x=0.0)
	
vp.BRight<- viewport(height=unit(.33, "npc"), width=unit(0.5, "npc"), just=c("left","top"), y=0.333, x=0.5)
	
print(g1.dis, vp=vp.MLeft)
print(g1.T, vp=vp.MRight)
print(g2.hja, vp=vp.BLeft)
print(g2.old, vp=vp.BRight)
	
dev.off()
	
# plot II
pdf('R/graph/nmds_geo_nmds3_3_env_Is2_noSpike.pdf', width=18, height=18)
	
#par(mfrow=c(3,3))
	
# . II. geo plots
plot.new()
vps <- baseViewports()
	
print(g1.dis, vp=vp.MLeft)
print(g1.T, vp=vp.MRight)
print(g2.hja, vp=vp.BLeft)
print(g2.old, vp=vp.BRight)

pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(.5,5,5,5), "null"))))
grid.text("kelpie 20200214, session 2, Malaise 1 (91 sites), no spike", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontfamily='sans'))
print(g1.ele, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(g3.ele, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
	
print(g1.rain, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(g3.rain, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
	
print(g2.road, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(g3.road, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
	
dev.off()
	
# ...............................................................
# ................. NMDS with env.variable, with spike......................
# . S2, OUT1 (91 sites) with spike .
names(dataI.2.1.spike)[1:27]
	
distI2.S = vegdist(wisconsin(dataI.2.1.spike[,-(1:26)]))
str(distI2.S)
# 91*45 = 4095
	
set.seed(7756)
count.nmdsI2_3.S = metaMDS(distI2.S, distance='bray', try=30, trymax=2000, k=3, trace=T,plot=T)
#... Procrustes: rmse 0.00187576  max resid 0.00793357 
	
set.seed(11756)
count.nmdsI2_4.S = metaMDS(distI2.S, distance='bray', try=30, trymax=2000, k=4, trace=T,plot=T)
#... Procrustes: rmse 0.001843379  max resid 0.009314212
	
# .. fit environmental data ..
fitI2_4.S = envfit(count.nmdsI2_4.S~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.2.1.spike, na.rm=T, choices=c(1:4))
fitI2_4.S
	
fitI2_3.S = envfit(count.nmdsI2_3.S~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.2.1.spike, na.rm=T, choices=c(1:3))
fitI2_3.S
	
# compare 4&3 dimension plots
pdf('R/graph/nmds_4_3_env_Is2_Spike.pdf', width=18, height=18)
	
par(mfrow=c(3,3))
	
# 3 dim
plot(count.nmdsI2_3.S, type='n', pch=19, choices=c(1,2))
points(count.nmdsI2_3.S, display='sites', pch=19)
plot(fitI2_3.S, choices=c(1,2)) 
	
plot(count.nmdsI2_3.S, type='n', pch=19, choices=c(2,3), main=paste("kelpi 20200214, 3 mds, Spike, session 2 (M1), stress ",round(count.nmdsI2_3.S$stress,4), sep=''))
points(count.nmdsI2_3.S, display='sites', pch=19)
plot(fitI2_3.S, choices=c(2,3)) 
	
plot(count.nmdsI2_3.S, type='n', pch=19, choices=c(1,3))
points(count.nmdsI2_3.S, display='sites', pch=19)
plot(fitI2_3.S, choices=c(1,3)) 
	
# 4 dim
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(1,2))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(1,2)) 
	
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(2,3), main=paste("kelpi 20200214, 4 mds, Spike, session 2 (M1, sites ",dim(dataI.2.1.spike)[1],"), stress ",round(count.nmdsI2_4.S$stress,4), sep=''))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(2,3)) 
	
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(1,3))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(1,3)) 
	
# 4 dim
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(1,4))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(1,4)) 
	
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(2,4), main=paste("kelpi 20200214, 4 mds, Spike, session 2 (M1, sites ",dim(dataI.2.1.spike)[1],"), stress ",round(count.nmdsI2_4.S$stress,4), sep=''))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(2,4)) 
	
plot(count.nmdsI2_4.S, type='n', pch=19, choices=c(3,4))
points(count.nmdsI2_4.S, display='sites', pch=19)
plot(fitI2_4.S, choices=c(3,4)) 
	
dev.off()
	






# ...............................................................
# ................. NMDS with env.variable, with spike......................
# . S1, OUT1 (88 sites) with spike .
distI1.S = vegdist(wisconsin(dataI.1.1.spike[,-(1:26)]))
str(distI1.S)
#88*87/2=3828
	
set.seed(7756)
count.nmdsI1.S = metaMDS(distI1.S, distance='bray', try=30, trymax=2000, k=2, trace=T,plot=T)
#  Procrustes: rmse 0.0009470243  max resid 0.005603913
	
set.seed(11756)  # set the seed so that the result can be repeated
count.nmdsI1_3.S = metaMDS(distI1.S, distance='bray', try=30, trymax=2000, k=3, trace=T,plot=T)
#  Procrustes: rmse 0.001026232  max resid 0.006034556
	
count.nmdsI1.S
# Stress:     0.2744694
count.nmdsI1_3.S
# Stress:     0.2123643 
	
# ..... fit env.data .....   
fitI1_3.S = envfit(count.nmdsI1_3.S~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.1.1.spike, na.rm=T, choices=c(1:3))
	     
fitI1_3S.2 = envfit(count.nmdsI1_3.S~yrs.disturb.min +elevation+max.T +precipitation+min.T+ hja+canopy.ht, data=dataI.1.1.spike, na.rm=T, choices=c(1:3))
	
fitI1_3S.2
plot(count.nmdsI1_3.S, choices=c(1,2))
plot(fitI1_3.S, choices=c(1,2)) 
	
# .. plot stress ..
count.nmdsI1_3.S$stress
count.nmdsI1_3.noS$stress
	
pdf('R/graph/nmds_kelpi0214_with.no.spike.pdf',height=5,width=10)
	
par(mfrow=c(1,2))
	
stressplot(count.nmdsI1_3.S, main='Spike, kelpie0214, session 1, malaise 1')
text(.765,.8, paste('stess:',round(count.nmdsI1_3.S$stress, 4),sep=''))
stressplot(count.nmdsI1_3.noS, main='no Spike, kelpie0214, session 1, malaise 1')
text(.75,.8, paste('stess:',round(count.nmdsI1_3.noS$stress, 4),sep=''))
	
dev.off()
	
# .. make file for plotting ..
nmds.geo.Is1_3.S = data.frame(UTM_E = dataI.1.1.spike$UTM_E, UTM_N = dataI.1.1.spike$UTM_N, nmds1 = count.nmdsI1_3.S$points[,1], nmds2 = count.nmdsI1_3.S$points[,2],nmds3 =count.nmdsI1_3.S$points[,3], ele=dataI.1.1.spike$elevation, disturb=dataI.1.1.spike$yrs.disturb.min,hja=as.character(dataI.1.1.spike$hja),min.T=dataI.1.1.spike$min.T,preci=dataI.1.1.spike$precipitation,c.ht = dataI.1.1.spike$canopy.ht)
str(nmds.geo.Is1_3.S)
	
# . NMDS & hja, disturb .
	
pdf('R/graph/nmds_geo_nmds2_env_Is1II_0225_Spike.pdf', width=18, height=11)
	
par(mfrow=c(2,3))
	
# . I. nmds with env.variables, plot 1:3
plot(count.nmdsI1_3.S, type='n', pch=19, choices=c(1,2))
points(count.nmdsI1_3.S, display='sites', pch=19)
plot(fitI1_3.S, choices=c(1,2)) 
	
plot(count.nmdsI1_3.S, type='n', pch=19, choices=c(1,3), xlim=c(min(nmds.geo.Is1_3.S$nmds1)*1.1,max(nmds.geo.Is1_3.S$nmds1)*1.1), ylim=c(min(nmds.geo.Is1_3.S$nmds3)*1.3, max(nmds.geo.Is1_3.S$nmds3)*1.3), main='spike, session 1, trap 1')
points(count.nmdsI1_3.S, display='sites', pch=19)
plot(fitI1_3.S, choices=c(1,3)) 
	
plot(count.nmdsI1_3.S, type='n', pch=19, choices=c(2,3), xlim=range(nmds.geo.Is1_3.S$nmds2), ylim=c(min(nmds.geo.Is1_3.S$nmds3)*1.25, max(nmds.geo.Is1_3.S$nmds3)*1.25))
points(count.nmdsI1_3.S, display='sites', pch=19)
plot(fitI1_3.S, choices=c(3,2)) 
	
# . II. geo plots
plot.new()              ## suggested by @Josh
vps <- baseViewports()
#pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
#vp1 <-plotViewport(c(.1,.1,0,.1))
	
g.Is1.1 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=disturb, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('disturbance') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='disturb')
g.Is1.3 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, shape=hja, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('hja') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.2 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=c.ht, colour=nmds2)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('canopy height') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='bottom', legend.title=element_text(size=9)) 
	
vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.333, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.67)
	
vp.BottomLeft <- viewport(height=unit(.5, "npc"), width=unit(0.333, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.0)
	
vp.BottomMiddle<- viewport(height=unit(.5, "npc"), width=unit(0.333, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.335)
	
print(g.Is1.1, vp=vp.BottomLeft)
print(g.Is1.2, vp=vp.BottomMiddle)
print(g.Is1.3, vp=vp.BottomRight)
	
dev.off()
	
# . nmds & T, ele, rain .
g.Is1.1 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=ele, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.2 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=ele, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, elevation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
g.Is1.3 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=preci, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.4 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=preci, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, precipitation') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
g.Is1.5 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=min.T, colour=nmds1)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds1, (min)temperature') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
g.Is1.6 = ggplot(nmds.geo.Is1_3.S, aes(UTM_E, UTM_N, size=min.T, colour=nmds3)) + geom_point() + scale_colour_gradient(low='purple', high='orange') + ggtitle('nmds3, (min)temperature') + theme(plot.title= element_text(hjust=.5), text=element_text(family = "sans"), legend.position='right', legend.title=element_text(size=9)) #+labs(size='hja')
	
pdf('R/graph/nmds_geo_nmds_T_ele_preci_Is1II_0225_Spike.pdf', width=11, height=15)
	
pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(.5,5,5,5), "null"))))
grid.text("kelpie 20200214, session 1, Malaise 1 (88 sites), spike", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontfamily='sans'))
print(g.Is1.1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(g.Is1.2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
	
print(g.Is1.3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(g.Is1.4, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
	
print(g.Is1.5, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(g.Is1.6, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
	
dev.off()
	














# ..........................................................................
# ................... test if we need to do normalization 'envfit' ....................
mystd = function(x) { (x - mean(x, na.rm=T)) / sd(x, na.rm=T)}
	
a = mystd(dataI.1.1$yrs.disturb.min)
	
b = mystd(dataI.1.1$elevation)
	
par(mfrow=c(2,2))
	
hist(dataI.1.1$yrs.disturb.min)
hist(a)
hist(try1$yrs.disturb.min.1)
	
hist(dataI.1.1$elevation)
hist(b)
hist(try1$elevation.1)
	
str(dataI.1.1[,2:23])
	
data.env = dataI.1.1[,4:23]
	
str(data.env)
	shape=hja,
data.env[,21:36] = data.env[,-c(5:7,10)]
data.env[,21:36] = lapply(data.env[,21:36], mystd)
#try1 = lapply(data.env[,21:36], mystd)
	
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
names(data.env)[21:36] = paste('nor.',names(data.env)[-c(5:7,10,21:36)],sep='')
	
fitI1_3.nor = envfit(count.nmdsI1_3~nor.yrs.disturb.min+nor.canopy.ht+nor.elevation+max.T+nor.metre.road+nor.old.growth.str+nor.precipitation+nor.min.T+ nor.metre.stream+ nor.hja, data=data.env, na.rm=T, choices=c(1:3))
	
scores(fitI1, "vectors")

fitI1_3 = envfit(count.nmdsI1_3.noS~yrs.disturb.min+canopy.ht+elevation+max.T+metre.road+old.growth.str+precipitation+min.T+metre.stream+hja, data=dataI.1.1.noS, na.rm=T, choices=c(1:3))
	
scores(fitI1_3, "vectors")
	
pdf('R/graph/env.fit.compare.pdf', height=5, width=13)
	
par(mfrow=c(1,2))
	
plot(count.nmdsI1_3)
plot(fitI1_3.nor)
	
plot(count.nmdsI1_3.noS)
plot(fitI1_3shape=hja,)
	
nmds.geo.Is1_3.noS$hja = character(nmds.geo.Is1_3.noS$hja)
dev.off()
	
# ..........................................................................
