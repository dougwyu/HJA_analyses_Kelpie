## best models

# from manual tuning 

## sJSDM on ADA
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada
library(dplyr)

resFolder <-"code_sjSDM/r20210217/results"
if(!dir.exists(resFolder)) stop("res folder doesn't exist")

# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

### 1. Get data  #####




## 2. tuning results #####
qp <- read.csv("results/manual_tuning_sjsdm_5CV_M1S1_mean_AUC_pa_min_5_nSteps_1000.csv")
qp.best <- qp[which.max(qp$AUC.test),, drop = T]
qp.best

pa <- read.csv("results_sjSDM/tuning_YL_pa_newVars/manual_tuning_sjsdm_newVars_S1_M1_pa_min5_envDNN_20210125.csv")
pa.best <- pa[which.max(pa$AUC.test),,drop = T]

rm(qp, pa)

hidden <- list(c(50L,50L,10L), c(25L,25L,10L))


model_cpu_best <- sjSDM(Y = otu.pa.minocc,
                      env = DNN(data = env.vars.new, formula = ~.,
                                hidden = hidden[[pa.best$hidden]],
                                activation = "relu", 
                                bias = TRUE, 
                                lambda = qp.best$lambda.env,
                                alpha = qp.best$alpha.env,
                                dropout = qp.best$drop),
                      spatial = linear(data = Sp.data.scale, formula = ~0+UTM_E:UTM_N,
                                    lambda = qp.best$lambda.sp,
                                    alpha = qp.best$alpha.sp),
                      biotic = bioticStruct(lambda = qp.best$lambda.bio, alpha = qp.best$alpha.bio),
                      learning_rate=0.003,
                      se = TRUE, family=binomial("probit"),
                      sampling = 5000L, iter = 150L,
                      device = "cpu")




model_gpu_best <- sjSDM(Y = otu.pa.minocc,
                      env = DNN(data = env.vars.new, formula = ~.,
                                hidden = hidden[pa.best$hidden],
                                activation = "relu", 
                                bias = TRUE, 
                                lambda = qp.best$lambda.env,
                                alpha = qp.best$alpha.env,
                                dropout = qp.best$drop),
                      spatial = linear(data = Sp.data.scale, formula = ~0+UTM_E:UTM_N,
                                       lambda = qp.best$lambda.sp,
                                       alpha = qp.best$alpha.sp),
                      biotic = bioticStruct(lambda = qp.best$lambda.bio, alpha = qp.best$alpha.bio),
                      learning_rate=0.003,
                      se = TRUE, family=binomial("probit"),
                      sampling = 5000L, iter = 150L,
                      device = "gpu")


save(model_cpu_best, model_cpu_best, file = file.path(resFolder, "best_models.rdata"))

