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
select.percent = .8		# .7
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
				  step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L # 150L, 5000L
				 )
				saveRDS(list(model=model.train, random=data.frame('lambda.bioN'=lambda.bioN, 'alpha.bioN'=alpha.bioN, 'lambda.spN'=lambda.spN, 'alpha.spN'=alpha.spN)), here(outputpath,'sjsdm_general_outputs',sjsdmVfolder,'sjsdm-model-RDS', glue('s-jSDM_tuning.cv_{period}_{trap}_{abund}_{cvn}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_hidden{hiddenN}_{dropN}_{date.model.run}.RDS')) )
	
			}
		}
	}
}
	


