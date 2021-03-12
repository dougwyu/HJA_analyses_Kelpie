# functions to make code of sjsdm graphs shorter
# Yuanheng Li, May 22, 2020
# last modified: 20201207

# variation partition plot with diff colors
vp.plot <- function(data, x1, x2, ind, textM, textL) {
    mycolor <- brewer.pal(count(unique(data[, ind]))[[1]], 'Paired')
    ternaryplot(data[, x1:x2], 
                pch = 20, 
                col = mycolor,
                main = textM, 
                labels = "outside", 
                grid = "dashed", 
                cex = 0.5, 
                prop_size = 1, 
                dimnames_position = "corner" # edge
                ) 
    grid_legend("topright", 
                pch=20, 
                col=mycolor, 
                labels = unique(data[, ind][[1]]),
                title = textL
                )
    }


# extract min & max spp pairs of correlation
extract.minmax.cor = function (otu, sigma) {
	spp.names = colnames(otu)
	corvalue=data.frame(cor=numeric(), rowname=character(), colname=character())
	for (j in 1:dim(sigma)[1]) {#
		print(j)
		a = data.frame(cor=sigma[,j], rowname=as.character(seq(1, dim(sigma)[1])), colname=as.character(rep(j,dim(sigma)[1])))
		
		corvalue = rbind(corvalue, a)
		rm(a)
	}
	
	corvalue1 = subset(corvalue, as.character(rowname)!=as.character(colname))
	corvalue1 = corvalue1[order(corvalue1$cor),]
	corvalue.minmax = corvalue1[c(1:20,(dim(corvalue1)[1]-19):dim(corvalue1)[1]),]
	corvalue.minmax = corvalue.minmax[seq(1,40,2), ]
		
	corvalue.minmax$spp1 = spp.names[corvalue.minmax$rowname]
	corvalue.minmax$spp2 = spp.names[corvalue.minmax$colname]
	my.output <<- list(corvalue, corvalue.minmax)
	
	return (str(my.output))
}
	

# species covariance circular plot with min&max pairs of covariance
cov.circle = function (version.text, otu.text, sigma, otu.tbl, result) {
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
	
	plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
	text(x = -2, y = 6, pos = 3, xpd = NA, labels = paste(version.text))
	#text(x = -6, y = 5.7, pos = 3, xpd = NA, labels = "A", font = 2, cex = 1.5)
	
	# print a circle
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
		
	sjSDM::add_legend(viridis::viridis(11), angles = c(140,110),radius = 5.4)
	text(cos(deg2rad(123))*(lineSeq+1), sin(deg2rad(123))*(lineSeq+1.2), labels = "covariance", pos = 2, xpd = NA)
		
	sjSDM::add_legend(cols = cols, range = c(2, ncol(result$sigma)), angles = c(70,40),radius = 5.4)
	text(cos(deg2rad(53))*(lineSeq+1), sin(deg2rad(55))*(lineSeq+1.1), labels = "low to high", pos = 4, xpd = NA) 
	text(cos(deg2rad(64))*(lineSeq+1.3), sin(deg2rad(62))*(lineSeq+1.1), labels = paste(otu.text), pos = 4, xpd = NA) 
		
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
		
	
	
}
	

# species covariance circle with max environmental variable
cov.circle.env = function (version.text, evnames, otu.text, result, effect_comb, otu.tbl) {
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
	
	sigma = re_scale(result$sigma)[effect_comb_ind, effect_comb_ind]
	
	sigmas = sigma[upper.tri(sigma)]
	number=10
	upper = order(sigmas, decreasing = TRUE)[1:number]
	lower = order(sigmas, decreasing = FALSE)[1:number]
	cuts = cut(sigmas, breaks = seq(-1,1,length.out = 12))
	to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
	levels(cuts) = viridis::viridis(11)
	cuts = as.character(cuts)
	n = ncol(sigma)
	lineSeq = 3.5
	nseg = 100
	
	
	#Drawing figure parameter
	par( mar = c(1,2,2.1,2)+0.1)
	plot(NULL, NULL, xlim = c(-5,5), ylim =c(-5,5),pty="s", axes = F, xlab = "", ylab = "")
	text(x = -1.7, y = 5.5, pos = 3, xpd = NA, labels = paste(version.text), font = 2, cex = 1)
		
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
	#formula = ~mean.NDVI.scale+mean.EVI.scale+ mean.green.scale + mean.bright.scale + mean.wet.scale
	
	colourCount = length(unique(evnames))
	cols = sample(rainbow(colourCount))
	
	coords = data.frame(cbind(xx, yy, angles))
	effect_comb=effect_comb[,-c(3,4)]
	effect_comb2 = effect_comb
	effect_comb2[,2] = ff(effect_comb[,2])
	effect_comb2 = cbind(effect_comb2, effect_comb[,2])
	effect_comb2 = data.frame(effect_comb2)
			
	for(i in sort(unique(max_effects))) {
	  sub = coords %>% filter(effect_comb2$max_effects==i)
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
	  ###environment variable bar
#	 if(i == 8) polygon(x-0.1, y, xpd = NA,col = cols[i])
	#  else if(i == 5) polygon(x, y+0.55, xpd = NA,col = cols[i])
	#  else if(i == 7) polygon(x-0.6, y-0.2, xpd = NA,col = cols[i])
	#  else if(i == 9) polygon(x, y-0.55, xpd = NA,col = cols[i])
	  #else if(i == 12) polygon(x-0.2, y-0.5, xpd = NA,col = cols[i])
#	 else 
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
	  text(x= x[10]-0.2, y=-5.2, labels = ncol(result$sigma), pos = 4, xpd = NA)
	}
	text(x=-5.3, y=-5.3, labels = paste(otu.text), pos = 4, xpd = NA)
	
	
}
	

# code from Cai Wang, edited by Doug
coef.figure <- function(summary.p, result, minsize,
                        maxsize = 200000, taxon = "all") {
  effect <- data.frame(summary.p$coefmat)
  effect$rownames <- rownames(effect)
  # effect <- rownames_to_column(effect, var = "rownames")
  effect <- separate(
    data = effect,
    col = rownames,
    into = c("species", "coef"),
    sep = " "
  )

  ifelse(taxon != "all",
    effect <- filter(effect, grepl(taxon, species)),
    effect <- effect
  )

  effect <- effect %>% filter(coef != "(Intercept)")
  effect <- effect %>% dplyr::select(-c(Z.value, Pr...z..))
  effect$coef <- as.factor(effect$coef)
  effect$max <- NA
  # str(effect$coef)
  n <- ncol(result$sigma) # number of species
  m <- ncol(data.frame(result$beta$env)) - 1
  for (i in 1:n) {
    max_effects <- as.numeric(
      apply(
        data.frame(effect[((m * i - (m - 1)):(m * i)), 1]), 2,
        function(e) which.max(abs(e))
      )
    )
    max_effects <- max_effects + ((i - 1) * m)
    effect$max[max_effects] <- "max"
  }
  effect$max[is.na(effect$max)] <- " "
  
  effect$group <- NA
  effect$group[grep("Diptera", effect$species)] <- "Diptera"
  effect$group[grep("Hymenoptera", effect$species)] <- "Hymenoptera"
  effect$group[grep("Araneae", effect$species)] <- "Araneae"
  effect$group[grep("Psocodea", effect$species)] <- "Psocodea"
  effect$group[grep("Hemiptera", effect$species)] <- "Hemiptera"
  effect$group[grep("Coleoptera", effect$species)] <- "Coleoptera"
  effect$group[grep("Neuroptera", effect$species)] <- "Neuroptera"
  effect$group[grep("Lepidoptera", effect$species)] <- "Lepidoptera"
  effect$group[grep("Lepidoptera", effect$species)] <- "Lepidoptera"
  effect$group[grep("Raphidioptera", effect$species)] <- "Raphidioptera"
  effect$group[grep("Orthoptera", effect$species)] <- "Orthoptera"
  effect$group[grep("Opiliones", effect$species)] <- "Opiliones"
  effect$group[grep("Plecoptera", effect$species)] <- "Plecoptera"
  # effect$coef <- as.factor(effect$coef)
  # levels species

  effect <- effect %>%
    tidyr::separate(species, into = c("OTU", "taxon"), sep = "__", remove = TRUE) %>%
    tidyr::separate(taxon, into = c("taxon", "size"), sep = "_size.", remove = TRUE) %>%
    # mutate(taxon = str_replace(taxon, "_sp", "_genus_sp")) %>%
    tidyr::separate(taxon, into = c("class", "order", "family", "genus", "species", "BOLD", "BOLDID")) %>%
    unite(OTU, c("OTU", "class", "order", "family", "genus", "species"), sep = "_", ) %>%
    select(-BOLD, -BOLDID) %>%
    mutate(size = as.numeric(size)) %>%
    rename(species = OTU) %>%
    filter(size > minsize & size < maxsize)

  effect$species <- as.factor(effect$species)
  effect$species <- factor(effect$species, levels = rev(levels(effect$species)))
  # effect$species=factor(effect$species,levels = levels(effect$species)[c(26:31,33:40,15:25,41:47,1:10,11:14,32)])
  # effect$species <- relevel(effect$species, "Wafl_Pele_Arde_Ardecine")

  # reorder levels coef
  # if (length((levels(effect$coef)))==13)  effect$coef=factor(effect$coef, levels = c('HSI2_Pond_area','HSI3_Pond_drying','HSI4_Water_quality','HSI5_Shade','HSI8_Pond_count','HSI10_Macrophytes','Inflow_present','Outflow_present','PC1','PC2','PC3','PC4','PC5'))
  #
  # #  else effect$coef=factor(effect$coef,levels = c('OS_area','cropN','cropP','cropK'))
  # effect <- effect %>% filter(effect$max == "max")
  t <- -(abs(effect$Estimate[which.max(abs(effect$Estimate))]))
  # t=-0.8
  # figure
  col <- brewer.pal(12, "Paired")
  # col <- sample(viridis::viridis(length(unique(effect_ind$max_effects))))
  p1 <- ggplot(
    effect,
    aes(
      x = species,
      y = Estimate,
      fill = group
    )
  ) +
    geom_bar(
      position = position_dodge(0.6),
      stat = "identity", width = 0.5, size = 0.5
    ) +
    scale_fill_manual(values = col) +
    guides(fill = guide_legend(reverse = F)) +
    xlab("species") +
    ylab("coef") +
    labs(fill = "Order") +
    coord_flip(expand = F) +
    geom_hline(aes(yintercept = 0),
      linetype = "dotted",
      size = 0.5
    ) +
    theme_classic() +
    facet_wrap(~coef, ncol = 5) +
    geom_text(aes(y = t, label = max),
      position = position_dodge(0.6),
      size = 1, fontface = "bold"
    ) +
    # geom_errorbar(aes(ymax = Estimate + Std.Err, ymin = Estimate - Std.Err), width = 0.3) +
    scale_y_continuous(limits = c(-1.1, 1)) +
    theme(axis.text.x = element_text(size = 2)) +
    theme(axis.text.y = element_text(size = 3))
  return(p1)
}


# list relation of spp & their envir factors
bar.coef = function(effect_ind, t_turn) {
    t= -(abs(effect_ind$beta[which.max(abs(effect_ind$beta))]))*t_turn
    #figure
    ggplot(effect_ind, aes(x = name, y = beta, fill=max_effects)) +
        geom_bar(stat="identity", width = 0.5)+
        #    position = position_dodge(0.4), 
        scale_fill_manual(values=col) +
        guides(fill = guide_legend(reverse=F)) +
        coord_flip(expand=F)  + 
        geom_hline(aes(yintercept = 0),linetype="dashed",size=1) +
        theme_classic()+ facet_wrap(~factors, ncol = 5) +
        xlab("otu") +
        ylab("beta") + 
        geom_text(aes(y= t, label = max),  size = 2, fontface = "bold")
    #    labs(fill="Species Index") 
    #    +geom_errorbar(aes(ymax = Estimate + Std.Err, ymin = Estimate - Std.Err), width = 0.3)
    #geom_errorbar(aes(ymax = Estimate + Std.Err, ymin = Estimate - Std.Err), width = 0.3)
}

