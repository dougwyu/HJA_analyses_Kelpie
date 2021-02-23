# Feb 22, 2021
# train = 's1m1', test = 's1m2'

```{r setup}
# setwd('/media/yuanheng/SD-64g2/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/13_sjsdm_prediction')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan','pROC', 'gridExtra','ggeffects','cowplot')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.3.9000’
	

```


```{r set-names}
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
trap <- "M1"; period = "S1"
date.model.run = 20210221   # according to tuning result
abund = 'pa'		# 'qp','pa'  !!! change accordingly
formula.env = 'newVars'		# 'envDNN.newVars', 'changed.envDNN'
minocc = 5
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	
```


```{r load-data}
# ..... load data ......
alldata = read.csv(here(outputpath, glue('sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_FSL_qp.csv')), header=T, sep=',', stringsAsFactors = F, na.strings='NA')
dim(alldata)
names(alldata)[110:125]
# 124 Var
	
alldata1 = alldata %>% filter(period == period[1]) #& trap == trap[1]
dim(alldata1); unique(alldata1$period)
	
minocc
a = alldata1 %>% select(contains('__'))
a = a[,which(vegan::specnumber(a, MARGIN=2)>=minocc)]
dim(a)
	
alldata1 = cbind(select(alldata1, -contains('__')), a)
dim(alldata1)
	
envvars = alldata1 %>% dplyr::select(!contains("__"), UTM_E, UTM_N, -starts_with("nor")) %>% 
	mutate(uniqueID = paste(SiteName, trap, period, sep = "_"),
         elevation_m = elevation_f * 0.3048, ## convert to metres
         canopyHeight_m = canopyHeight_f * 0.3048,
         lg_DistStream = log(distToStream_m + 0.001),
         lg_DistRoad = log(distToRoad_m + 0.001),
         lg_YrsDisturb = log(YrsSinceDist + 0.001),
         lg_cover2m_max = log(l_Cover_2m_max + 0.001),
         lg_cover2m_4m = log(l_Cover_2m_4m + 0.001),
         lg_cover4m_16m = log(l_Cover_4m_16m + 0.001)) %>% mutate(clearcut = factor(clearcut), insideHJA = factor(insideHJA))
str(envvars)
	

# ..... make tables (test->m2, train->m1) .....
ttindex = which(alldata1$trap %in% 'M1')
	
oldVars <- c("insideHJA", "elevation_f", "canopyHeight_f", "minT_annual", "precipitation_mm", "distToRoad_m", "distToStream_m", "YrsSinceDist", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_Cover_2m_max", "l_Cover_2m_4m", "l_Cover_4m_16m", "l_p25", "l_p95", "l_rumple")
	
newvars <- c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")
	
ori.env = select(envvars, newvars)
str(ori.env)
	
ori.env.train = ori.env[ttindex,]
ori.env.test = ori.env[-ttindex,]
str(ori.env.test)
	
ori.XY = select(alldata1, starts_with('UTM'))
ori.XY.train = ori.XY[ttindex,]
ori.XY.test = ori.XY[-ttindex,]
str(ori.XY.train)
	
otu = select(alldata1, contains('__'))
otu.train = otu[ttindex,]
dim(otu.train)
otu.test = otu[-ttindex,]
	
# .. env data
scale.env.train.all = select(ori.env.train, -'insideHJA') %>% scale() 
str(scale.env.train.all)
	
scale.env.train = data.frame(scale.env.train.all) %>% add_column(insideHJA=as.factor(ori.env.train$insideHJA), .before=names(ori.env.train)[2])
# insideHJA factor !!!
str(scale.env.train)
	
dd.env.scaler = data.frame(t(data.frame(env.mean = attr(scale.env.train.all, "scaled:center"), env.sd = attr(scale.env.train.all, "scaled:scale"))))
str(dd.env.scaler)
	
rm(scale.env.train.all)
	
scale.env.test = as.data.frame(do.call(rbind, apply(select(ori.env.test, -'insideHJA'), 1, function(x){(x-dd.env.scaler['env.mean',])/dd.env.scaler['env.sd',]} ) )) %>% add_column(insideHJA=as.factor(ori.env.test$insideHJA), .before=names(ori.env.test)[2])
str(scale.env.test)
	
# .. spatial data
XY.train.all = scale(ori.XY.train)
str(XY.train.all)
	
XY.train = data.frame(XY.train.all)
str(XY.train)
	
dd.xy.scaler = data.frame(t(data.frame(sp.mean = attr(XY.train.all, "scaled:center"), sp.sd = attr(XY.train.all, "scaled:scale"))))
str(dd.xy.scaler)
base::rownames(dd.xy.scaler)
	
rm(XY.train.all)
	
XY.test = as.data.frame(do.call(rbind, apply(ori.XY.test, 1, function(x){(x-dd.xy.scaler['sp.mean',])/dd.xy.scaler['sp.sd',]} ) ))
str(XY.test)
	

```


```{r save-data}
datapath=glue('../12_sjsdm_general_model_outputs')
save(XY.train,XY.test,scale.env.train,scale.env.test,otu.train,otu.test, file=here(datapath,'source',glue('forada_data_{period}_m1m2_min{minocc}_{date.model.run}_{formula.env}.rdata')))
	
```


```{r explore-graph}
# ..... explore m1,m2 ......
exp = select(alldata1, 'UTM_E', 'UTM_N', 'trap', 'SiteName') %>% bind_cols(spprich=rowSums(otu>0))
str(exp)
	
which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'])
which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'] & exp$trap=='M1')

# pdf(here(outputpath,'prediction_outputs','sjsdm-graph', sjsdmVfolder,glue('coordinate_s1_m1m2_min{minocc}.pdf')), height=10,width=10)
	
p1 = ggplot(exp, aes(UTM_E, UTM_N),shape=trap) + geom_point(aes(colour=trap,shape=trap)) 
p2 = ggplot(exp, aes(UTM_E, UTM_N),colour=spprich) + geom_point(aes(colour=spprich)) + geom_point(colour='white') + geom_point(data=filter(exp,trap=='M2'), aes(colour=spprich)) + scale_colour_gradient2(low='black', mid='red',high='blue',midpoint=max(exp$spprich)/2) + ggtitle('M2,test') + theme(legend.position='none')
p3 = ggplot(exp, aes(UTM_E, UTM_N),colour=spprich) + geom_point(aes(colour=spprich)) + geom_point(colour='white') + geom_point(data=exp[which(exp$SiteName %in% exp$SiteName[exp$trap=='M2'] & exp$trap=='M1'),], aes(colour=spprich)) + scale_colour_gradient2(low='black', mid='red',high='blue',midpoint=max(exp$spprich)/2) + ggtitle('M1, explain')
	
grid.arrange(p1,p2,p3, nrow=2, layout_matrix= rbind(c(1), c(2,3)), heights=c(.4,.6),widths=c(.44,.56))
	
dev.off()
	

```
