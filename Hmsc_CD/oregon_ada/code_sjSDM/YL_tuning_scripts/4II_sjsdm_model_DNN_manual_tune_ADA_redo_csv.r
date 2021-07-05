
## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()

# Jan 18, 2021
# Jan 25, 2021
# code on sjsdm model run with DNN (for ADA cluster)


## REDO lambdaENV == 1 prediction csv results were overwritten
## EXTRACT AUC from rds


## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
getwd() # always run sub from oregon_ada

options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] 0.1.3.9000

# lapply(c('tidyverse','reticulate','sjSDM','glue','vegan'), library, character.only=T)


library(glue)
library(pROC)

abund = 'pa'		# 'qp','pa'  # is overwritten by data file loading.
minocc = 5
trap <- "M1"; period = "S1"
date.model.run = 20210125   # !!! change accordingly

### r load data### 
# otu.train <- read.csv(paste0("data/otu.train.", abund, ".csv"), header = F)
# scale.env.train <- read.table("data/scale.env.train.csv", sep = " ", header = T)
# XY.train <- read.table("data/scale.XY.train.csv", sep = ",", header = T)

# head(otu.train)
# str(otu.train)
# 
# head(scale.env.train)
# str(scale.env.train)
# 
# head(XY.train)

# from here
# "J:\UEA\gitHRepos\HJA_analyses_Kelpie\Hmsc_CD\oregon_ada\data\otu.train.pa.csv"
# "J:\UEA\gitHRepos\HJA_analyses_Kelpie\Hmsc_CD\oregon_ada\data\otu.train.qp.csv"
# "J:\UEA\gitHRepos\HJA_analyses_Kelpie\Hmsc_CD\oregon_ada\data\scale.env.train.csv"
# "J:\UEA\gitHRepos\HJA_analyses_Kelpie\Hmsc_CD\oregon_ada\data\scale.XY.train.csv"

# data creatd here and saved: Hmsc_CD|oregon_ada|code_sjSDM|sjSDM_trial|ada_sjsdm1_data_setup.r
# OLD VARS
load(paste0("data/yuanheng_mod_data_", abund, ".rdata"))
# s.otu.train,scale.env.train, XY.train, s.otu.test, scale.env.test, XY.test, abund


# NEW VARS
# load(paste0("data/yuanghen_mod_data_newVars_", abund, ".rdata"))


### r model-DNN.env### 
# make sure input of sjsdm are numeric matrix
s.otu.train = as.matrix(s.otu.train)
attr(s.otu.train, 'dimnames') = NULL
str(s.otu.train)

names(scale.env.train)

# set variables
formula.env = 'envDNN'
lambda.env = seq(0,.3, length.out=4)	# .1
alpha.env = seq(.7,1, length.out=4)		# .9
lambda.sp = seq(0,1, length.out=7)	# .1
alpha.sp =  seq(0,1, length.out=7)	# .5 
hidden = list(c(50L,50L,10L), c(25L,25L,10L))
acti.sp = 'relu'
drop = seq(.1,.5, length.out=3) # .3
sample.bio = seq(0,1,length.out=11)

# no models
4*7*7*2*3  # each lambda

# data storage
# 4*7*7*2*3 * 642
# time
# 4*7*7*2*3*2 * .5 # 162 models - 1 hour

tuning.dd = data.frame(lambda.env = numeric(),
                       alpha.env = numeric(),
                       lambda.sp = numeric(),
                       alpha.sp = numeric(),
                       lambda.bio = numeric(), 
                       alpha.bio = numeric(),
                       drop = numeric(), 
                       hidden = character(),
                       loglike = numeric(), 
                       loss= numeric(), 
                       AUC.explain=numeric(), 
                       AUC.test=numeric())


# hiddenN <- dropN <- alpha.spN <- lambda.spN <- alpha.envN <- lambda.envN <- 1

# lambda.envN = 1		# 1,2,3,4 four jobs in AD

# testing
# for( lambda.envN in 1){
#   for (alpha.envN in 1) {
#     for (lambda.spN in 1) {
#       for (alpha.spN in 1) {
#         for (dropN in 1) {
#           for (hiddenN in 1) {
            
for(lambda.envN in seq_along(lambda.env)){
  for (alpha.envN in seq_along(alpha.env)) {
    for (lambda.spN in seq_along(lambda.sp)) {
      for (alpha.spN in seq_along(alpha.sp)) {
        for (dropN in seq_along(drop)) {
          for (hiddenN in seq_along(hidden)) {
            
            # lambda.bioN = sample(1:11,1)
            # alpha.bioN = sample(1:11,1)
            # print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN, hiddenN))
            # 
            
            # Old vars            
            mod <- readRDS(file.path("results_sjSDM", paste0("tuning_YL_", abund),
                                     glue::glue('s-jSDM_tuning_model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}.RDS')
            ))

            # new vars            
            # mod <- readRDS(file.path("results_sjSDM", paste0("tuning_YL_", abund, "_newVars"),
            #                          glue::glue('s-jSDM_tuning_model__newVars_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}.RDS')
            # ))
            
            
            # contains:
            # list(model=model.train, random=data.frame('lambda.bioN'=lambda.bioN, 'alpha.bioN'=alpha.bioN))
            
            model.train <- mod$model
            lambda.bioN <- mod$random$lambda.bioN
            alpha.bioN <- mod$random$alpha.bioN
            
            # Do testing and save results in data frame
            tdd = data.frame(lambda.env = lambda.env[lambda.envN], 
                             alpha.env = alpha.env[alpha.envN], 
                             lambda.sp = lambda.sp[lambda.spN], 
                             alpha.sp = alpha.sp[alpha.spN], 
                             lambda.bio = sample.bio[lambda.bioN],
                             alpha.bio = sample.bio[alpha.bioN], 
                             drop = drop[dropN],
                             hidden = as.character(hiddenN), 
                             loglike = logLik(model.train), 
                             loss= model.train$history[length(model.train$history)], 
                             AUC.explain=.1, AUC.test=.1)
            
            for (pred in 1:2) {
              
              # 1 -> 'test'
              newdd = scale.env.test ; newsp = XY.test; otudd = s.otu.test
              
              if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train}
              
              pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train, newdata=newdd, SP=newsp)),
                                           along = -1L), 2:3, mean)
              
              attr(pred.dd, 'dimnames') = NULL
              
              if (pred==1) {
                otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
                pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
                otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
              }
              
              otudd.pa = (otudd>0)*1
              roc.dd = lapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i], direction = "<", quiet=T))
              auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) stringr::str_split(pROC::auc(i), ':')), function(x) x[[1]][1] )))
              
              if (pred==2) {tdd$AUC.explain=auc.mean}
              if (pred==1) {tdd$AUC.test=auc.mean}
            }
            
            tuning.dd = rbind(tuning.dd, tdd)
            
            print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN))
            
            rm(model.train, tdd)
            
            # old vars
            write.table(tuning.dd,
                        file = file.path("results_sjSDM",paste0("tuning_YL_", abund),
                                         glue::glue('manual_tuning_sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')),
                        row.names=F, sep=',')
            
            # new vars
            # write.table(tuning.dd, 
            #             file = file.path("results_sjSDM",paste0("tuning_YL_", abund, "_newVars"), 
            #                              glue::glue('manual_tuning_sjsdm_newVars_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')),
            #             row.names=F, sep=',')
            
          }
        }
      }
    }
  }
}


