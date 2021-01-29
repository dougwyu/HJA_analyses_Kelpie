
## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# dir()

# Jan 18, 2021
# Jan 25, 2021
# code on sjsdm model run with DNN (for ADA cluster)

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
getwd() # always run sub from oregon_ada

options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’

# lapply(c('tidyverse','reticulate','sjSDM','glue','vegan'), library, character.only=T)

library(glue)
library(pROC)

abund = 'qp'		# 'qp','pa'  # is overwritten by data file loading.
minocc = 5
trap <- "M1"; period = "S1"
date.model.run = 20210125   # !!! change accordingly

resFolder <- "results_sjSDM/tuning_YL2"

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
load("data/yuanghen_mod_data.rdata")
# s.otu.train,scale.env.train, XY.train, s.otu.test, scale.env.test, XY.test, abund


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


# hiddenN <- dropN <- alpha.spN <- lambda.spN <- alpha.envN <- 1

lambda.envN = 3		# 1,2,3,4 four jobs in AD

# testing
# for (alpha.envN in 1:2) {
#   for (lambda.spN in 1) {
#     for (alpha.spN in 1) {
#       for (dropN in 1) {
#         for (hiddenN in 1) {

for (alpha.envN in seq_along(alpha.env)) {
  for (lambda.spN in seq_along(lambda.sp)) {
    for (alpha.spN in seq_along(alpha.sp)) {
      for (dropN in seq_along(drop)) {
        for (hiddenN in seq_along(hidden)) {

          lambda.bioN = sample(1:11,1)
          alpha.bioN = sample(1:11,1)
          print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN, hiddenN))
          
          model.train = sjSDM(
            Y = s.otu.train,
            env = DNN(data=scale.env.train, formula = ~.,
                      hidden=hidden[[hiddenN]],
                      lambda = lambda.env[lambda.envN],
                      alpha = alpha.env[alpha.envN],
                      activation=acti.sp,
                      dropout=drop[dropN],
                      bias=T),
            biotic = bioticStruct(lambda=sample.bio[lambda.bioN],
                                  alpha=sample.bio[alpha.bioN], 
                                  on_diag=F, inverse = FALSE),
            spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, 
                             lambda=lambda.sp[lambda.spN], 
                             alpha=alpha.sp[alpha.spN]),
            learning_rate = 0.003, # 0.003 recommended for high species number 
            step_size = NULL, iter = 150L, 
            family=stats::binomial('probit'), 
            sampling = 5000L, # 150L, 5000L
            device = "gpu"
            )
          
          saveRDS(
            list(model=model.train, random=data.frame('lambda.bioN'=lambda.bioN, 'alpha.bioN'=alpha.bioN)),
            file.path(resFolder,
                      glue::glue('s-jSDM_tuning_model_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}.RDS')
                      )
            )
          
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
# 
#             if (pred==1) {
#               otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
#               pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
#               otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
#             }
 
             otudd.pa = (otudd>0)*1

#             roc.dd = lapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i]))
#             auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) stringr::str_split(pROC::auc(i), ':')), function(x) x[[1]][1] )))
            
            auc <- tryCatch(expr = sapply(1:ncol(pred.dd), function(i) Metrics::auc(otudd.pa[,i],pred.dd[,i])),
                     error = function(err){ return(NA) })
            
            auc.mean <- mean(auc, na.rm= T)

            if (pred==2) {tdd$AUC.explain=auc.mean}
            if (pred==1) {tdd$AUC.test=auc.mean}
          }
          
          tuning.dd = rbind(tuning.dd, tdd)
          
          print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN))
          
          rm(model.train, tdd)
          
          write.table(tuning.dd, 
                      file = file.path(resFolder, 
                  glue::glue('manual_tuning_sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}_lambdaE{lambda.envN}.csv')),
                      row.names=F, sep=',')
          
          }
      }
    }
  }
}


# ..... for test run .....
# lambda.env = .1
# alpha.env = .9
# lambda.sp = .1
# alpha.sp = .5 
# hidden = c(25L,25L,10L)
# acti.sp = 'relu'
# drop = .3
# 
# model.train = sjSDM(Y = s.otu.train,
#                     env = DNN(data=scale.env.train, formula = ~.,
#                               hidden=hidden, lambda = lambda.env, alpha = alpha.env, activation=acti.sp, dropout=drop, bias=T),
#                     
#                     biotic = bioticStruct(lambda=lambda.sp, alpha=alpha.sp, on_diag=F, inverse = FALSE),
#                     
#                     spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp, alpha=alpha.sp),
#                     
#                     learning_rate = 0.003, # 0.003 recommended for high species number 
#                     step_size = NULL, iter = 5L, family=stats::binomial('probit'), sampling = 50L # 150L, 5000L
# )
# # ..... for test run .....
# 
# model.train = readRDS(here(tempsjsdmpath,'results','sjsdm-model-RDS', sjsdmVfolder, glue('s-jSDM_tuned.model_{period}_{trap}_{abund}_min{minocc}_{date.model.run}_{formula.env}.RDS')) )
# 
# plot(model.train$history)





### r model-save-less ###
#hiddenN=1; lambda.envN=2; alpha.envN=2; lambda.spN=2; alpha.spN=2;dropN=2

# set variables (from 'model-...' chunk)
# tuning.dd = data.frame(lambda.env = numeric(), alpha.env = numeric(),  lambda.sp = numeric(), alpha.sp = numeric(), lambda.bio = numeric(), alpha.bio = numeric(), drop = numeric(), hidden = character(), loglike = numeric(), loss= numeric(), AUC.explain=numeric(), AUC.test=numeric())
# 
# lambda.envN = 1		# 1,2,3,4 four jobs in AD
# hiddenN = 1			# 1,2 4*2 jobs
# 
# for (alpha.envN in 1:3) {
#   for (lambda.spN in 1:3) {
#     for (alpha.spN in 1:3) {
#       for (dropN in 1:3) {
#         
#         lambda.bioN = sample(1:11,1)
#         alpha.bioN = sample(1:11,1)
#         
#         model.train = sjSDM
#         (Y = s.otu.train,
#           env = DNN(data=scale.env.train, formula = ~.,
#                                       hidden=hidden[[hiddenN]], lambda = lambda.env[lambda.envN], alpha = alpha.env[alpha.envN], activation=acti.sp, dropout=drop[dropN], bias=T),
#                             
#          biotic = bioticStruct(lambda=sample.bio[lambda.bioN], alpha=sample.bio[alpha.bioN], on_diag=F, inverse = FALSE),
#          spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp[lambda.spN], alpha=alpha.sp[alpha.spN]),
#          learning_rate = 0.003, # 0.003 recommended for high species number 
#          step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L # 150L, 5000L
#         )
#         
#         tdd = data.frame(lambda.env = lambda.env[lambda.envN], alpha.env = alpha.env[alpha.envN], lambda.sp = lambda.sp[lambda.spN], alpha.sp = alpha.sp[alpha.spN], lambda.bio = sample.bio[lambda.bioN], alpha.bio = sample.bio[alpha.bioN], drop = drop[dropN], hidden = as.character(hiddenN), loglike = logLik(model.train), loss= model.train$history[length(model.train$history)], AUC.explain=.1, AUC.test=.1)
#         
#         for (pred in 1:2) {
#           # 1 -> 'test'
#           newdd = scale.env.test ; newsp = XY.test; otudd = otu.test
#           if (pred==2) { newdd = NULL; newsp = NULL; otudd = otu.train }
#           
#           pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train, newdata=newdd, SP=newsp)) , along = -1L), 2:3, mean)
#           attr(pred.dd, 'dimnames') = NULL
#           
#           if (pred==1) {
#             otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
#             pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
#             otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
#           }
#           
#           otudd.pa = (otudd>0)*1
#           roc.dd = lapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i]))
#           auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) stringr::str_split(pROC::auc(i), ':')), function(x) x[[1]][1] )))
#           
#           if (pred==2) {tdd$AUC.explain=auc.mean}
#           if (pred==1) {tdd$AUC.test=auc.mean}
#         }
#         
#         tuning.dd = rbind(tuning.dd, tdd)
#         print(c(lambda.envN, alpha.envN, lambda.spN, alpha.spN, dropN))
#         
#         rm(model.train, tdd)
#         
#         write.table(tuning.dd, file=here(outputpath, 'sjsdm_general_outputs', sjsdmVfolder, 'DNN_tune', glue('manual_tuning.sjsdm_{period}_{trap}_{abund}_min{minocc}_{formula.env}_{date.model.run}.csv')), row.names=F, sep=',')
#         
#         
#       }
#     }
#   }
# }
# 
# 
# 
# 
