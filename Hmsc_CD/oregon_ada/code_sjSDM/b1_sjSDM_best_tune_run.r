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

resFolder <-"results_sjSDM"
if(!dir.exists(resFolder)) stop("res folder doesn't exist")


# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# try basic model on oregon data
source("code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train
rm(P)

# Data reduction - reduce species
raretaxa <- which(colSums(otu.pa.csv > 0) < 10)
length(raretaxa)
otu.pa.minocc <- as.matrix(otu.pa.csv[, -raretaxa]) # reduced species
rm(raretaxa)


# scale all togher for now... just for test
env.vars.scale <- env.vars %>%
  mutate(across(where(is.numeric), scale))
# head(env.vars.scale)

env.vars.old <- env.vars.scale[, c("insideHJA", "elevation_m", "canopyHeight_m", "minT_annual", "precipitation_mm", "lg_DistRoad", "lg_DistStream", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m", "l_p25", "l_p95", "l_rumple")]

# new vars
env.vars.new <- env.vars.scale[, c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")]


# XFormula1 <- as.formula(~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo)
# check names
# all(all.vars(XFormula1) %in% names(env.vars.scale))


# spatial data here:
# spatial data here:
head(Sp.data)
# scale spatial coords
Sp.data.scale <- scale(Sp.data[, c("UTM_E", "UTM_N")])
head(Sp.data.scale)
is.matrix(Sp.data.scale)


## tuning results
qp <- read.csv("results_sjSDM/tuning_YL_qp_newVars/manual_tuning_sjsdm_newVars_S1_M1_qp_min5_envDNN_20210125.csv")
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

