# Yuanheng
# create: Apr 27, 2020
# last modified:  Apr , 2020
# apply sj-sdm & hmsd 
# with data (in'data_Feb_25/', 'kelpie20200214')


# .................................................................
rm(list=ls())
setwd("/media/yuanheng/SD-64g2/files/Projects/Oregon")
getwd()
	
source('R/source/sjsdm_function.r')
	
lapply(c("ggplot2", "gridExtra",'vegan', 'labdsv','tidyverse','scatterplot3d', 'gridBase','grid', 'ggcorrplot'),library,character.only=T)
	
lapply(c('Hmsc','sjSDM', 'reticulate'),library, character.only=T)
	
# ...................... kelpie, lidar data ..............................
# (lidar data, load data formatted for sjsdm)
lidar.env = read.csv('kelpie/lidar_data/sjsdm_biodiversity_site_info_multispectral_2020-04-13.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
str(lidar.env)
	
# ('data_Feb_25' folder, load data formatted for sjsdm) 
otu.env1.noS = read.csv('kelpie/data_Feb_25/lidar_sample_by_species_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
otu.env1.spike = read.csv('kelpie/data_Feb_25/lidar_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA')
	
print(c(dim(otu.env1.spike), dim(otu.env1.noS)))
# 1173-26, more otus
	
# ......................... subsets of data ..................................
dataI.1.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M1' )
dataI.2.spike = subset(otu.env1.spike, session == 'S1' & trap == 'M2' )
print(c(dim(dataI.1.spike), dim(dataI.2.spike)))
	
dataII.1.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M1' )
dataII.2.spike = subset(otu.env1.spike, session == 'S2' & trap == 'M2' )
print(c(dim(dataII.1.spike), dim(dataII.2.spike)))
	
# .... if there's all zero OTUs ....
b = data.frame(otu=colnames(dataI.1.spike[,-c(1:36)]), zero=apply(dataI.1.spike[,-c(1:36)],2,sum)==0)
b$otu=as.character(b$otu)
	
dataI.1.spike2 = dplyr::select(dataI.1.spike, 1:36,b$otu[b$zero==F])
dataI.1.spike2[1:5,34:40]
	
########################################################################
# ............................. analyzing result .......................
########################################################################
# Max:
# * wrong "loss", the first value is the loglik of the model; second one is the regularization loss
# .... analyse the model ....
# s1, m1, spiked data
result<-readRDS("R/result/s-jSDM_result_s1_m1_spike2.RDS")
# "s-jSDM_result_otu_env1_spike.RDS"
	
str(result)
	
lrs = seq(-18, -1, length.out = 7);f = function(x) 2^x;lrs = f(lrs)
	
# ... look at the loss
result[[1]]$loss
	
loss.all<-data.frame(t(sapply(result, function(r) r$loss[1])))
loss.all<-data.frame(LogLik = unlist(loss.all))
str(loss.all)
loss.all$lrs<-lrs
	
loss.2 = data.frame(t(sapply(result, function(r) r$loss)))
loss.2 = unlist(loss.2)
loss.2<-data.frame(LogLik = loss.2[seq(1,length(loss.2),by=2)], unknown=loss.2[seq(2,length(loss.2),by=2)])
data.frame(loss.2,loss.all[,2])
	
plot(x=loss.all$lrs,y=loss.all$LogLik,ylab = "LogLik",xlab = "lrs")
	
# . cov <- result$sigma 
dim(result[[5]]$sigma)
# [1] 1147 1147
co.env.spp <- cov2cor(result[[4]]$sigma)
spp.names<-colnames(dataI.1.spike2[,-c(1:36)])
rownames(co.env.spp)<-spp.names
colnames(co.env.spp)<-spp.names
	
rownames(co.env.spp) <- 1:850
colnames(co.env.spp) <- 1:850
	
ggcorrplot(co.env.spp[1:400,1:400], hc.order = T, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,show.legend=T)   
# [1:100,1:100]
	
#a = apply(abind::abind(lapply(result[c(10,15)], function(r) cov2cor(r$sigma)), along = 0L), 2:3, sum)
#a[1:5,1]
	
#a1=abind::abind(cov2cor(result[[10]]$sigma), along = 0L)
#a1[1,1,1:5]
	
overall = apply(abind::abind(lapply(result, function(r) cov2cor(r$sigma)), along = 0L), 2:3, sum)
	
# sum of all models with diff lrs
pdf('R/graph/sjsdm_s1_m1_spike2/7lrs_spp_corr.pdf', height=20, width=20)
	
ggcorrplot(cov2cor(overall), hc.order = TRUE, outline.color = "white", insig = "blank",sig.level = 0.05, lab_size = 1,title=paste('lrs ',formatC(lrs[1],format='e',3),' to ',lrs[7],', sum of 7 models',sep=''))  
	
dev.off()
	
# . heatmap, don't understand .
for { # for lrs[4]
co.env.spp4<-cov2cor(result[[4]]$sigma)
range(co.env.spp4)
	
pdf(paste('R/graph/sjsdm_s1_m1_spike2/heatmap_lrs',round(lrs[1],digit=8),'_spp_corr.pdf', sep=''), height=20, width=20)
	
indices = heatmap(co.env.spp4) # [1:200,1:200] 
	
dev.off()
	
str(indices)
	
}
co.env.spp1<-cov2cor(result[[1]]$sigma)
co.env.spp7<-cov2cor(result[[7]]$sigma)
	
pdf('R/graph/sjsdm_s1_m1_spike2/spp_corr_heatmap_3lrs.pdf', height=5, width=15)
	
par(mfrow=c(1,3))
cols = (colorRampPalette(c("blue", "white", "red")))(10)
graphics::image(co.env.spp1[indices$rowInd, indices$rowInd], col = cols, main=paste('lrs ',formatC(lrs[1],format='e',3),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp1),4),' - ', max(co.env.spp1),sep=''), ylab='blue (min value) to red (max value)')
image(co.env.spp4[indices$rowInd, indices$rowInd], col = cols, main=paste('species correlation (850 spp), lrs ',formatC(lrs[4],format='e',3),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp4),4),' - ', max(co.env.spp4),sep=''))
image(co.env.spp7[indices$rowInd, indices$rowInd], col = cols, main=paste('lrs ',round(lrs[7],8),sep=''),xaxt = "n", yaxt = "n", xlab=paste('spp corr ',round(min(co.env.spp7),4),' - ', max(co.env.spp7),sep=''))
	
dev.off()
	
# ..... Polygon Drawing .....
# . species association .
for { # result[1], lrs 
number=10
lr_step=1
sigma = re_scale(result[[lr_step]]$sigma)[order(apply(dataI.1.spike2[,-c(1:36)], 2, sum)), order(apply(dataI.1.spike2[,-c(1:36)], 2, sum))]
# re_scale -> cov2cor
#select.otu = base::sample(1:1147,300)
#sigma = sigma[select.otu,select.otu]
sigmas = sigma[base::upper.tri(sigma)]
upper = order(sigmas, decreasing = TRUE)[1:number]
lower = order(sigmas, decreasing = FALSE)[1:number]
cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
summary(cuts)
	
to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
levels(cuts) = viridis::viridis(11)
cuts = as.character(cuts)
n = ncol(result[[lr_step]]$sigma) # [select.otu,select.otu]
lineSeq = 4.7
nseg = 100
	
pdf(paste('R/graph/sjsdm_s1_m1_spike2/spp_cov_lrs',formatC(lrs[7],format=NULL,1),'.pdf',sep=''), height=8, width=8)
	
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = 0, y = 5.7, pos = 3, xpd = NA, labels = paste("Penalty: ",formatC(lrs[lr_step],format=NULL,1),sep=''))
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
      #cat(counter, "\n")
  }
}
# ??? correlation plotted, why legend says 'covariance' 
	
OTU_log = log(sort(apply(dataI.1.spike2[, -c(1:36)], 2, sum)))
range(OTU_log)
OTU_log[1]=0 # ???
cuts = cut(OTU_log, breaks = 10)
cols = viridis::magma(10) #colfunc(5)
	
levels(cuts) = cols
sppnames2=paste("spp",1:20,sep = "")
abun=as.character(cuts)
sppsort=1:20
	
OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=dim(dataI.1.spike2[,-c(1:36)])[2]), sum = apply(dataI.1.spike2[,-c(1:36)], 2, sum))
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
  #curve_text(coords[i,3], label = spnames2[i],reverse = T,lineSeq =lineSeq +0.2, middle = TRUE, extend = 1.2, col = as.character(cuts[i]))
}
	
add_legend(viridis::viridis(11), angles = c(140,110),radius = 5.4)
text(cos(deg2rad(123))*(lineSeq+1), sin(deg2rad(123))*(lineSeq+1.2), labels = "covariance", pos = 2, xpd = NA)
	
add_legend(cols = cols, range = c(2, 850), angles = c(70,40),radius = 5.4)
# ??? range = c(2, 151), ???2:300
text(cos(deg2rad(53))*(lineSeq+1), sin(deg2rad(55))*(lineSeq+1.1), labels = "Sp. incidence", pos = 4, xpd = NA) 
text(cos(deg2rad(64))*(lineSeq+1.3), sin(deg2rad(62))*(lineSeq+1.1), labels = "sum of abun", pos = 4, xpd = NA) 
	
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
	
}
# ..... environmental effect .....
# Drawing parameter of OTU and environmental covariant
for { # result[[1]]
lr_step = 7
#formula = ~ elevation.scale+canopy.ht.scale+min.T.scale+max.T.scale+precipitation.scale+metre.road.scale+metre.stream.scale+yrs.disturb.min.scale,  
# 9 variables + intercept
effects= apply(result[[lr_step]]$beta[[1]], 1, function(o) sum(abs(o)))
turn_over = 1
n = ncol(result[[lr_step]]$sigma)# number of otus
max_effects= apply(result[[lr_step]]$beta[[1]],2, function(e) which.max(abs(e)))
turn_over = 1
effect_comb = data.frame(cbind(max_effects,sapply(1:n, function(i) result[[lr_step]]$beta[[1]][max_effects[i],i] )))
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
	
pdf(paste('R/graph/sjsdm_s1_m1_spike2/max_environ_lrs',formatC(lrs[lr_step],format=NULL,1),'.pdf',sep=''), height=8, width=8)
	
#Drawing figure parameter
#layout(matrix(c(1,2,3,rep(4,9)), 3,3 ,byrow = F), c(0.9,1.3,1.3), c(1,1,1))
par( mar = c(1,2,2.1,2)+0.1)
plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
text(x = -3.5, y = 5.3, pos = 3, xpd = NA, labels = paste("Penalty: ",formatC(lrs[lr_step],format='e',3), ', loss: ', loss.all[lr_step,1],sep=''), font = 2, cex = 1)
	
xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
angles = seq(0,360,length.out = n+1)[1:(n)] # for all otus
xx = cos(deg2rad(angles))*lineSeq
yy = sin(deg2rad(angles))*lineSeq
	
##inside circle
counter = 1
coords = cbind(xx, yy, angles)
#spnames=paste("sp",1:22,sep="")
	
for(i in 2:n) {
  for(j in 1:(i-1)){
    
      if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, species = T, lineSeq = 3.5, lwd = 1.3)
      counter = counter + 1
      #cat(counter, "\n")
  }
}
#for(i in 1:n) {
#curve_text(coords[i,3], label = effect_comb$name[i],reverse = T,lineSeq = 4, middle = TRUE, extend = 1.2, col = effect_comb$abun[i])
#}
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
evnames=c('intercept',"ele","canopy","min.T","max.T","preci","road","stream","disturb")
colourCount = length(unique(evnames))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
cols=getPalette(colourCount)
#cols = RColorBrewer::brewer.pal(12,"Set3")
coords = data.frame(cbind(xx, yy, angles))
effect_comb=effect_comb[,-c(3,4)]
effect_comb2 = effect_comb
effect_comb2[,2] = ff(effect_comb[,2])
# ff -> (x-min(x))/(max(x)-min(x))
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
  
  #names(habitat2)
  angleName = (from+to)/2
  if(angleName > 180) {reverse = TRUE} else {reverse = FALSE}
  ###environment variable text
  curve_text(angleName, label = evnames[i],reverse = reverse,lineSeq = 5.5, middle = TRUE, extend = 1.1, col = cols[i])
  ###environment variable bar
#  if(i == 3) polygon(x+0.4, y, xpd = NA,col = cols[i])
#  else if(i == 5) polygon(x, y+0.55, xpd = NA,col = cols[i])
#  else if(i == 7) polygon(x-0.6, y-0.2, xpd = NA,col = cols[i])
#  else if(i == 9) polygon(x, y-0.55, xpd = NA,col = cols[i])
  #else if(i == 12) polygon(x-0.2, y-0.5, xpd = NA,col = cols[i])
#  else 
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
	
abun=as.character(abun)
x = seq(-5.5,-3, length.out = 11)
for(i in 1:unique(length(abun))){
  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = unique(abun)[i], xpd = NA, border = NA)
  text(x= x[1]-0.2, y=-5.2, labels = "2", pos = 4, xpd = NA)
  text(x= x[10]-0.2, y=-5.2, labels = '850', pos = 4, xpd = NA)
}
text(x=-5.3, y=-5.3, labels = "Sp. incidence", pos = 4, xpd = NA)
	
dev.off()
	
}










