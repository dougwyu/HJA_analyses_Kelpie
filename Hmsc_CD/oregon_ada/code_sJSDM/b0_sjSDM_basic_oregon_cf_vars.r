
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


# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# try basic model on oregon data
source("code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train
rm(P, otu.qp.csv)

# Data reduction - reduce species
raretaxa <- which(colSums(otu.pa.csv > 0) < 10)
length(raretaxa)
otu.pa.minocc <- as.matrix(otu.pa.csv[, -raretaxa]) # reduced species
rm(raretaxa)


# scale all togher for now... just for test
env.vars.scale <- env.vars %>%
  mutate(across(where(is.numeric), scale))
head(env.vars.scale)

env.vars.old <- env.vars.scale[, c("insideHJA", "elevation_m", "canopyHeight_m", "minT_annual", "precipitation_mm", "lg_DistRoad", "lg_DistStream", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m", "l_p25", "l_p95", "l_rumple")]

# new vars
env.vars.new <- env.vars.scale[, c("be10", "tri", "slope", "Nss", "Ess", "ht", "ht.r250", "ht.r500", "ht.r1k", "cov2_4", "cov2_4.r250", "cov2_4.r500", "cov2_4.r1k", "cov4_16", "cov4_16.r250", "cov4_16.r500", "cov4_16.r1k", "be500", "mTopo", "cut.r1k.pt", "insideHJA", "minT_annual", "maxT_annual", "precipitation_mm", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "B1_20180717", "B2_20180717", "B3_20180717", "B4_20180717", "B5_20180717", "B6_20180717", "B7_20180717", "B10_20180717", "B11_20180717", "NDVI_20180717", "EVI_20180717", "B_20180717", "G_20180717", "W_20180717", "l_p25", "l_rumple")]


length(env.vars.new)
length(env.vars.old)

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

# this is what goes into the model
# mm <- model.matrix(XFormula1, data = env.vars.scale)
# head(mm)

# model <- sjSDM(Y = otu.pa.minocc,
#                env = linear(data = env.vars.scale,
#                             formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 +
#                               cov4_16 + cov4_16.r500 + mTopo), # linear model on env covariates
#                spatial = linear(data = Sp.data.scale, formula = ~0+UTM_E:UTM_N), # interactions of coordinates
#                se = TRUE, family=binomial("probit"), sampling = 1000L, iter = 100L,
#                device = "gpu")
# 
# # model summary
# summary(model)
# 
# # coefficients
# # coef(model)
# 
# imp = importance(model)
# print(imp)
# # 
# pdf("results_sjSDM/trial_oregon_sJSDM_importance.pdf")
# plot(imp)
# dev.off()
# 
# an = anova(model)
# print(an)
# 
# pdf("results_sjSDM/trial_oregon_sJSDM_anova.pdf")
# plot(an)
# dev.off()
# 
# save(model, file ="results_sjSDM/oregon_trial.rdata")

## Full tuning parameters
# lambda.env = seq(0,.3, length.out=4)	# .1
# alpha.env = seq(.7,1, length.out=4)		# .9
# lambda.sp = seq(0,1, length.out=7)	# .1
# alpha.sp =  seq(0,1, length.out=7)	# .5 

## tuning not included in sjSDM_cv()
# hidden = list(c(50L,50L,10L), c(25L,25L,10L))
# acti.sp = 'relu'
# drop = seq(.1,.5, length.out=3) # .3
# sample.bio = seq(0,1,length.out=11)

# sjSDM_cv
results_old_vars = sjSDM(Y = otu.pa.minocc,
                         env = linear(data = env.vars.old,
                                   formula = ~.),
                         spatial = linear(Sp.data.scale, ~0 + UTM_E:UTM_N),
                         biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
                         se = TRUE,
                         learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
                         family = stats::binomial("probit"), # for both p/a and quasiprob data, default
                         device = "gpu",
                         iter = 100L, # 2L
                         sampling = 1000L # default is 5000L
)



# ncores
# NULL, # or 10L for small models, run this many sjsdm models at once (1 per CPU core). 10L is too much for large models because the 10 CPUs try to run on only the 2 GPUs available on ada: 5/GPU

# sjSDM_cv
results_new_vars <- sjSDM(Y = otu.pa.minocc,
                          env = linear(data = env.vars.new, 
                                       formula = ~.),
                          spatial = linear(Sp.data.scale, ~0 + UTM_E:UTM_N),
                          biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' 
                          se = TRUE,
                          learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
                          family = stats::binomial("probit"), # for both p/a and quasiprob data, default
                          device = "gpu",
                          iter = 100L, # 2L
                          sampling = 1000L # default is 5000L
)


summ.old <- summary(results_old_vars)
summ.new <- summary(results_new_vars)


# save results
save(results_new_vars, results_old_vars, summ.old, summ.new, file ="results_sjSDM/oregon_trial_cf_vars.rdata")
