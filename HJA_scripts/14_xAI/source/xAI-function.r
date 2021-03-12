
# species covariance circular plot with min&max pairs of covariance
cov.circle.env = function (version.text, evnames, otu.text, effect_comb, otu.tbl) {
	OTU_log = log(sort(apply(otu.tbl, 2, sum))+.001)
	range(OTU_log)
	OTU_log[1]=0 
	cuts = cut(OTU_log, breaks = 10)
	cols = viridis::magma(10) 
	
	levels(cuts) = cols
	sppnames2=paste("spp",1:20,sep = "")
	abun=as.character(cuts)
	sppsort=1:20
	
	OTU_sort_abun <- data.frame(sppsort=rep(sppsort, len=ncol(otu.tbl)), sum = apply(otu.tbl, 2, sum))
	OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sum), ]
	OTU_sort_abun$abun<-abun
	OTU_sort_abun = OTU_sort_abun[order(OTU_sort_abun$sppsort), ]
	
	sppname3=seq(1:20); sppname3=as.character(sppname3)
	effect_comb$name <- rep(sppname3, len=nrow(effect_comb))
	effect_comb$abun <- OTU_sort_abun$abun
	effect_comb$abun<-as.character(effect_comb$abun)
	effect_comb_ind = order(effect_comb[,1], effect_comb[,2])
	effect_comb = effect_comb[effect_comb_ind,]

	n = ncol(otu.tbl)
	lineSeq = 3.5
	nseg = 100

	#Drawing figure parameter
	par( mar = c(1,2,2.1,2)+0.1)
	plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
	text(x = -1.7, y = 5.5, pos = 3, xpd = NA, labels = paste(version.text), font = 2, cex = 1)

	# print a circle
	xx = lineSeq*cos( seq(0,2*pi, length.out=nseg) )
	yy = lineSeq*sin( seq(0,2*pi, length.out=nseg) )
	polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
	angles = seq(0,360,length.out = n+1)[1:(n)]
	xx = cos(deg2rad(angles))*lineSeq
	yy = sin(deg2rad(angles))*lineSeq

	coords = cbind(xx, yy, angles)

	lineSeq = 4.0
	for(i in 1:n){
	  p1 = coords[i,]
	  x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
	  y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
	  segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = effect_comb$abun[i], lend = 1)
	}
	lineSeq = 3.5
	
	##outside circle 

	colourCount = length(unique(evnames))
	cols = sample(rainbow(colourCount))
	
	coords = data.frame(cbind(xx, yy, angles))
	effect_comb=effect_comb[,-c(3,4)]
	effect_comb2 = effect_comb
	effect_comb2[,2] = ff(effect_comb[,2])
	effect_comb2 = cbind(effect_comb2, effect_comb[,2])
	effect_comb2 = data.frame(effect_comb2)
			
	for(i in sort(unique(effect_comb2$max_effects))) {
	 # sub = coords %>% filter(effect_comb2$max_effects==i)
	  sub = coords[which(effect_comb2$max_effects==i),]
	  sub_eff = effect_comb2 %>% filter(max_effects==i)
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
	
# no need .............................................	
	###legend of bar  
# no need .............................................		

	#abun=as.character(abun)
	abun1 = as.character(viridis::magma(10))
	x = seq(-5.5,-3, length.out = 11)
	for(i in 1:unique(length(abun))){
	  rect(xleft = x[i], xright = x[i+1], ybottom = -5, ytop = -5+diff(x)[1], col = abun1[i], xpd = NA, border = NA)
	  text(x= x[1]-0.2, y=-5.2, labels = "2", pos = 4, xpd = NA)
	  text(x= x[10]-0.2, y=-5.2, labels = ncol(otu.tbl), pos = 4, xpd = NA)
	}
	text(x=-5.3, y=-5.3, labels = paste(otu.text), pos = 4, xpd = NA)
	
}

	
		
	






	
	
