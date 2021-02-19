# Feb 02, 2021
# code on sjsdm model run with DNN (on laptop)



# ...r setup
# setwd('/media/yuanheng/SD-64g2/Downloads/backup2/HJA_analyses_Kelpie/HJA_scripts/12_sjsdm_general_model_outputs')
	
pacman::p_load('tidyverse','here','conflicted','reticulate','sjSDM','glue','vegan','pROC')
	
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')
	
here()
packageVersion('sjSDM')
#[1] ‘0.1.3.9000’
	



# ...r set-names
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
date.model.run = 20210215   # !!! change accordingly
	
formula.env = 'envDNN-cv'			# 'envDNN-newgis'
	
outputidxstatstabulatefolder = glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")
outputpath = glue('../../Kelpie_maps/{outputidxstatstabulatefolder}')
	
sjsdmV = '0.1.3.9000' # package version
	
# names for graph
sjsdmVfolder = glue('sjsdm-{sjsdmV}')
	


# ... r read-data
# ...... from rdata ........
cvn = 'cv1'				# cv2 , cv3
select.percent = .7		# .7
load( here('source',glue('YL_ada_{formula.env}_{select.percent}_{period}_{trap}_{minocc}_{abund}'),glue('YL_ada_{formula.env}_{select.percent}_{period}_{trap}_{minocc}_{abund}_{cvn}.rdata')))
	
# .. generate 's.' data
s.otu.train = as.matrix(otu.train.cv1)		# otu.train.cv2 , otu.train.cv3
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)
	
s.otu.valid = as.matrix(otu.valid.cv1)		# otu.valid.cv2 , otu.valid.cv3
attr(s.otu.valid, 'dimnames') = NULL
str(s.otu.valid)
	
scale.env.train = scale.env.train.cv1		# scale.env.train.cv2 , scale.env.train.cv3
scale.env.valid = scale.env.valid.cv1		# scale.env.valid.cv2 , scale.env.valid.cv3 
XY.train = XY.train.cv1			# XY.train.cv2 , XY.train.cv3
XY.valid = XY.valid.cv1			# XY.valid.cv2 , XY.valid.cv3
	
str(s.otu.train); abund
str(scale.env.valid); formula.env
	


# ... r model
# make sure input of sjsdm are numeric matrix
str(s.otu.train)
str(scale.env.train)
	
# set variables
formula.env 
lambda.env = seq(0,1, length.out=11)	# .1
alpha.env = seq(0,1, length.out=11)		# .9
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.1,.4, length.out=4) # .3
sample.spbio = seq(0,1,length.out=11)
	
#4^2*7^2*2*3
# 10*11*2*4+2*4
	
set.seed(55)
lseed = sample(1:5000, (10*11*2*4+2*4))
ccc = 0
	
for (lambda.envN in 1:11) {
	for (alpha.envN in 1:11) {
		for (dropN in 1:4) {
			for (hiddenN in 1:2) {
				
				if (lambda.envN == 1 & alpha.envN >1) {next}
				ccc = ccc + 1
				set.seed(lseed[ccc])
				lambda.bioN = sample(1:11,1); alpha.bioN = sample(1:11,1)
				lambda.spN = sample(1:11,1); alpha.spN = sample(1:11,1)
				
#				print(c(ccc, lambda.bioN, alpha.bioN,lambda.spN, alpha.spN))
				
				model.train = sjSDM(Y = s.otu.train,
				  env = DNN(data=scale.env.train, formula = ~.,
				  hidden=hidden[[hiddenN]], lambda = lambda.env[lambda.envN], alpha = alpha.env[alpha.envN], activation=acti.sp, dropout=drop[dropN], bias=T),
				  
				  biotic = bioticStruct(lambda=sample.spbio[lambda.bioN], alpha=sample.spbio[alpha.bioN], on_diag=F, inverse = FALSE),
				  
				  spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=sample.spbio[lambda.spN], alpha=sample.spbio[alpha.spN]),
				  
				  learning_rate = 0.003, # 0.003 recommended for high species number 
				  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 3000L # 150L, 5000L
				 )
				saveRDS(list(model=model.train, random=data.frame('lambda.bioN'=lambda.bioN, 'alpha.bioN'=alpha.bioN, 'lambda.spN'=lambda.spN, 'alpha.spN'=alpha.spN)), here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS',glue('cv_{date.model.run}'), glue('s-jSDM_tuning.cv_{period}_{trap}_{abund}_{cvn}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
	
			}
		}
	}
}
	

# ... r manually-tune
formula.env 
lambda.env = seq(0,1, length.out=11)	# .1
alpha.env = seq(0,1, length.out=11)		# .9
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.1,.4, length.out=4) # .3
sample.spbio = seq(0,1,length.out=11)
	
tuning.dd = data.frame(lambda.env = numeric(), alpha.env = numeric(),lambda.sp = numeric(), alpha.sp = numeric(), lambda.bio = numeric(), alpha.bio = numeric(), drop = numeric(), hidden = character(), loglike = numeric(), loss= numeric(), AUC.train=numeric(), AUC.valid=numeric(), ll.train=numeric(), ll.valid=numeric(), nagel.train=numeric(), nagel.valid=numeric(),plr.train=numeric(),plr.valid=numeric())
# plr->positive likelihood rate
	
# hiddenN=1; lambda.envN=2; alpha.envN=2; lambda.spN=2; alpha.spN=2;dropN=2
	
for (lambda.envN in 1:11) {
	for (alpha.envN in 1:11) {
		for (dropN in 1:4) {
			for (hiddenN in 1:2) {
						
		tryCatch({
			model.train = readRDS(here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS', glue('s-jSDM_tuning.cv_{period}_{trap}_{abund}_{cvn}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
#			plot(model.train$model$history)

			lambda.bioN = model.train$random$lambda.bioN; alpha.bioN = model.train$random$alpha.bioN
			lambda.spN = model.train$random$lambda.spN; alpha.spN = model.train$random$alpha.spN
			
			tdd = data.frame(lambda.env = lambda.env[lambda.envN], alpha.env = alpha.env[alpha.envN], lambda.sp = sample.spbio[lambda.spN], alpha.sp = sample.spbio[alpha.spN], lambda.bio = sample.spbio[lambda.bioN], alpha.bio = sample.spbio[alpha.bioN], drop = drop[dropN], hidden = as.character(hiddenN), loglike = logLik(model.train$model), loss= model.train$model$history[length(model.train$model$history)], AUC.train=.1, AUC.valid=.1,ll.train=.1, ll.valid=.1,nagel.train=.1, nagel.valid=.1,plr.train=.1,plr.valid=.1)
			
			for (pred in 1:2) {
				# 1 -> 'test'
				# AUC
				newdd = scale.env.valid ; newsp = XY.valid; otudd = s.otu.valid
				if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train }
				
				pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train$model, newdata=newdd, SP=newsp)) , along = -1L), 2:3, mean)
				attr(pred.dd, 'dimnames') = NULL
				
				if (pred==1) {
					otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<nrow(otudd))*1 )
					pred.dd = pred.dd[ ,which(otudd[nrow(otudd),] == 1)]
					otudd = otudd[1:(nrow(otudd)-1), which(otudd[nrow(otudd),] == 1)]
				}
				
				otudd.pa = (otudd>0)*1
				auc.mean = mean(sapply(1:ncol(otudd), function(i) as.numeric(pROC::roc(otudd.pa[,i], pred.dd[,i],direction='<')$auc)))
				
				if (pred==2) {tdd$AUC.train=auc.mean}
				if (pred==1) {tdd$AUC.valid=auc.mean}
				rm(auc.mean)
				
				# ll, nagel & plr for spp 
				rsq = data.frame(ll=rep(.1, length.out=ncol(pred.dd)), nagel=rep(.1, length.out=ncol(pred.dd)), plr=rep(.1, length.out=ncol(predy.dd)))
				for (j in 1:ncol(pred.dd)) { 
					p = pred.dd[,j]; y = otudd.pa[,j]
					loglikP = sum( log( p*y + (1-p)*(1-y) ) )
					loglikN = sum( log( mean(p)*y + (1-mean(p))*(1-y) ) )
					rsq$nagel[j] = (1-exp(2/length(p)*(loglikN-loglikP))) / (1-exp(2/length(p)*loglikN))
					rsq$ll[j] = loglikP
					tp = sum(p*y); fp = sum(p*(1-y)); fa = sum((1-p)*y); ta = sum((1-p)*(1-y))
					rsq$plr[j] = tp/(tp+fa)/fp*(fp+ta)
				}
				if (pred==2) {tdd$ll.train=mean(rsq$ll); tdd$nagel.train=mean(rsq$nagel); tdd$plr.train=mean(rsq$plr)}
				if (pred==1) {tdd$ll.valid=mean(rsq$ll); tdd$nagel.valid=mean(rsq$nagel); tdd$plr.valid=mean(rsq$plr)}
			}
			
		tuning.dd = rbind(tuning.dd, tdd)
		print(c(hiddenN, lambda.envN, alpha.envN, lambda.spN, alpha.spN, lambda.bioN, alpha.bioN, dropN))
		
		rm(model.train, tdd)
		
		write.table(tuning.dd, file=here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('manual.cv_tuning.sjsdm_{period}_{trap}_{abund}_min{minocc}_{cvn}_{formula.env}_{date.model.run}.csv')), row.names=F, sep=',')
	
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		
}}}}
	

