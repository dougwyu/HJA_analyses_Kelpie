### Map best AUC



wd <- here::here()
setwd(wd)
getwd()
# rm(wd)


# Jan 18, 2021

## Data set up for models on ADA

# code on sjsdm model run with DNN (for ADA cluster) - from Yuanheng
# to test how long it would take to run one model on the cluster

## I run this locally, save the data, copy it to ADA, then load the saved data file on ADA - Christian
## could also just take the minimap file straight from github online and then follow this script in ADA.

wd <- here::here()
wd # "J:/UEA/gitHRepos/HJA_analyses_Kelpie"  local
setwd(wd)

# lapply(c('tidyverse','reticulate','sjSDM','glue'), library, character.only=T) # removed here, conflicted

library(glue)
library(dplyr)


# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"

minocc = 5
abund = 'pa' # pa , qp    # !!! change accordingly
trap <- "M1"; period = "S1"
date.model.run = 20210119   # !!! change accordingly

outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('Kelpie_maps/{outputidxstatstabulatefolder}')

sjsdmV = '0.1.3.9000' # package version

# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')

# ..... load data ......

# what file am I loading?
basename(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')))

# when was it last modified?
file.mtime(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')))
# "2021-02-02 11:11:23 GMT"

alldata = read.csv(file.path(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
dim(alldata)
names(alldata)[1:10]

# ..... select trap, session .....
trap; period

alldata1 = subset(alldata, trap == 'M1' & period == 'S1')
dim(alldata1)

# select species with minocc
a = alldata1 %>% dplyr::select(contains('__'))
a = a[,which(vegan::specnumber(a, MARGIN=2)>=minocc)]
dim(a)

# join species back to data
alldata1 = cbind(dplyr::select(alldata1, -contains('__')), a)
dim(alldata1)

otu = dplyr::select(alldata1, contains('__'))
ori.XY = dplyr::select(alldata1, starts_with('UTM'))

otu <- (otu > 0)* 1

head(ori.XY)
head(otu)


load("Hmsc_CD/oregon_ada/results_sjSDM/best_auc.rdata")

str(auc, max.level = 1)



auc.train <- as.numeric(auc$train)
auc.test <- as.numeric(auc$test)

plot(auc.train, auc.test)
boxplot(auc.train, auc.test)

# join species to AUC values, and map

otu[1:10, 1:10]

# m <- matrix(c(1,0,0,0,1,1,1,1,0,0,1,1,1,0,0), ncol = 3)
# m
# x <- c(1,2,3)
# m * x
# sapply(1:ncol(m), function(i) m[,i] * x[i] )

auc_spp <- sapply(1:ncol(otu), function(i) otu[,i] * auc.test[i])
auc_spp[1:10, 1:10]

# convert to 0 to NA for metrics only where presnet
auc_spp[auc_spp == 0] <- NA

#av.auc <- rowSums(auc_spp, na.rm = T)/rowSums(auc_spp>0, na.rm = T)
av.auc <- apply(auc_spp, 1, mean, na.rm = T)
sd.auc <- apply(auc_spp, 1, sd, na.rm = T)

plot(av.auc, sd.auc/av.auc)

# auc.df <- data.frame(rbind(cbind(ori.XY, value = av.auc, type = "mean auc"), 
#                            cbind(ori.XY, value = sd.auc, type = "auc sd")))
# head(auc.df)


library(ggplot2)
library(cowplot)
# 
# ggplot(auc.df, aes(x = UTM_E, y = UTM_N, size = value))+
#   geom_point()+
#   facet_wrap(~type)

p1 <- ggplot(cbind(ori.XY, av.auc), aes(x = UTM_E, y = UTM_N, size = av.auc))+
  geom_point()+
  guides(size = guide_legend(title = "Average auc"))


p2 <- ggplot(cbind(ori.XY, sd.auc), aes(x = UTM_E, y = UTM_N, size = sd.auc))+
  geom_point()+
  guides(size = guide_legend(title = "SD auc"))

plot_grid(p1, p2)
ggsave2("Hmsc_CD/local/plots/auc_location1.png", width = 250, height = 150, units = "mm")


## Species with lowest AUC

cbind(colSums(otu), auc.test)

png("Hmsc_CD/local/plots/auc_prevalence.png")
plot(colSums(otu), auc.test, xlab = "prevalence")
dev.off()
cor.test(colSums(otu), auc.test, method = "spearman")

cbind(ori.XY, pa = otu[, which.min(auc.test)])

size = 3

m1 <- ggplot(cbind(ori.XY, pa = otu[, which.min(auc.test)]), aes(x = UTM_E, y = UTM_N, shape = as.factor(pa)))+
  geom_point(size = size)+
  #scale_color_manual(values = c("grey20", "darkred"))+
  scale_shape_manual(values = c(1, 16), name = "Presence/Absence")


m2 <- ggplot(cbind(ori.XY, pa = otu[, which.max(auc.test)]), aes(x = UTM_E, y = UTM_N, shape = as.factor(pa)))+
  geom_point(size = size)+
  scale_shape_manual(values = c(1, 16), name = "Presence/Absence")
  #scale_color_manual(values = c("grey20", "darkred"), name = "Presence/Absence")


## lowest 10% of AUC
dat10 <- cbind(ori.XY, pa = rowSums(otu[, which(auc.test < quantile(auc.test, probs = 0.1, na.rm = T))])>0 )

q1 <- ggplot(dat10, aes(x = UTM_E, y = UTM_N, shape = pa))+
  geom_point(size = size)+
  scale_shape_manual(values = c(1, 16), name = "Presence/Absence")
  #scale_color_manual(values = c("grey20", "darkred"), name = "Presence/Absence")

dat90 <- cbind(ori.XY, pa = rowSums(otu[, which(auc.test > quantile(auc.test, probs = 0.9, na.rm = T))])>0 )

q2 <- ggplot(dat90, aes(x = UTM_E, y = UTM_N, shape = pa))+
  geom_point(size = size)+
  scale_shape_manual(values = c(1, 16), name = "Presence/Absence")
  #scale_color_manual(values = c("grey20", "darkred"), name = "Presence/Absence")

plot_grid(m1, m2, q1, q2, labels = c("min AUC", "max AUC", "< 10 percentile AUC", "> 90 percentile AUC"))
ggsave2("Hmsc_CD/local/plots/auc_location2.png", width = 250, height = 150, units = "mm")
