# Feb 05, 2021
# code on sjsdm model run with DNN (for ADA cluster)



# setup
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.3.9000’
	



# set-names
# ....... folder structure .......
# bioinfo structure
samtoolsfilter = "F2308" # F2308 filter only
samtoolsqual = "q48"
minimaprundate = 20200929
kelpierundate = 20200927
primer = "BF3BR2"
	
minocc = 5
abund = 'qp'		# 'qp','pa'  !!! change accordingly
trap <- "M1"; period = "S1"
date.model.run = 20210125   # !!! change accordingly
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	


# load data
# ... training data
s.otu.train = as.matrix(read.table(here(PATHtoDATA, glue('otu.train.{abund}.csv')), sep=','))
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
scale.env.train = read.table(here(PATHtoDATA, 'scale.env.train.csv'), header=T, sep=',')
str(scale.env.train)
	
XY.train = read.table(here(PATHtoDATA, 'scale.XY.train.csv'), header=T, sep=',')
str(XY.train)
	
# ... test data
s.otu.test = as.matrix(read.table(here(PATHtoDATA, glue('otu.test.{abund}.csv')), sep=','))
attr(s.otu.test, 'dimnames') = NULL
str(s.otu.test)
	
scale.env.test = read.table(here(PATHtoDATA, 'scale.env.test.csv'), header=T, sep=',')
str(scale.env.test)
	
XY.test =read.table(here(PATHtoDATA, 'scale.XY.test.csv'), header=T, sep=',')
str(XY.test)
	



# calculate-auc
# set variables
formula.env = 'envDNN'		# change accordingly !!!
lambda.env = seq(0,.3, length.out=4)	# .1
alpha.env = seq(.7,1, length.out=4)		# .9
lambda.sp = seq(0,1, length.out=7)	# .1
alpha.sp =  seq(0,1, length.out=7)	# .5 
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.1,.5, length.out=3) # .3
sample.bio = seq(0,1,length.out=11)
	
# hiddenN=1; lambda.envN=2; alpha.envN=2; lambda.spN=2; alpha.spN=2;dropN=2
	
tuning.dd = data.frame(lambda.env = numeric(), alpha.env = numeric(), lambda.sp = numeric(), alpha.sp = numeric(), lambda.bio = numeric(), alpha.bio = numeric(), drop = numeric(), hidden = character(), loglike = numeric(), loss=numeric(), AUC.explain=numeric(), AUC.test=numeric())
	
for (hiddenN in 1:2) {
for (lambda.envN in 1:4) {
for (alpha.envN in 1:4) {
	for (lambda.spN in 1:7) {
		for (alpha.spN in 1:7) {
			for (dropN in 1:3) {
	
		tryCatch({
				model.train = readRDS(here(PATHTOsjsdmRDS, glue('s-jSDM_tuning.model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
			
				lambda.bioN = model.train$random$lambda.bioN
				alpha.bioN = model.train$random$alpha.bioN
				 
				tdd = data.frame(lambda.env = lambda.env[lambda.envN], alpha.env = alpha.env[alpha.envN], lambda.sp = lambda.sp[lambda.spN], alpha.sp = alpha.sp[alpha.spN], lambda.bio = sample.bio[lambda.bioN], alpha.bio = sample.bio[alpha.bioN], drop = drop[dropN], hidden = as.character(hiddenN), loglike = logLik(model.train$model), loss= model.train$model$history[length(model.train$model$history)], AUC.explain=.1, AUC.test=.1)
			
				for (pred in 1:2) {
				# 1 -> 'test'
					newdd = scale.env.test ; newsp = XY.test; otudd = s.otu.test
					if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train }
					
					pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train$model, newdata=newdd, SP=newsp)) , along = -1L), 2:3, mean)
					attr(pred.dd, 'dimnames') = NULL
					
					if (pred==1) {
						otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
						pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
						otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
					}
					
					otudd.pa = (otudd>0)*1
					roc.dd = lapply(1:dim(otudd)[2], function(i) roc(otudd.pa[,i], pred.dd[,i], direction = "<", quiet=T))
					auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) str_split(auc(i), ':')), function(x) x[[1]][1] )))
					
					if (pred==2) {tdd$AUC.explain=auc.mean}
					if (pred==1) {tdd$AUC.test=auc.mean}
				}
			
				tuning.dd = rbind(tuning.dd, tdd)
				print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN))
			
				rm(model.train, tdd)
			
				write.table(tuning.dd, file=here(PATHtoRESULT, glue('manual_tuning_sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')), row.names=F, sep=',')
	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			}
		}
	}
}	}}
	


