
## sJSDM on ADA

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada


# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# try basic model on oregon data
source("code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train

# Data reduction - reduce species
raretaxa <- which(colSums(Y.train.pa > 0) < 10)
length(raretaxa)
Y.train.pa_min10 <- as.matrix(Y.train.pa[, -raretaxa]) # reduced species
rm(raretaxa)

XFormula1 <- as.formula(~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo)
# check names
all(all.vars(XFormula1) %in% names(X.train))
str(X.train)


# spatial data here:
head(S.train)
xy.scale <- scale(S.train[,c("UTM_E", "UTM_N")])
# head(xy.scale)

# this is what goes into the model
# mm <- model.matrix(XFormula1, data = env.vars)
# head(mm)

# make and run model
model <- sjSDM(Y = Y.train.pa_min10,
               env = linear(data = env.vars, 
                            formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo), # linear model on env covariates
               spatial = linear(data = xy.scale, formula = ~0+UTM_E:UTM_N), # interactions of coordinates
               se = TRUE, family=binomial("probit"), sampling = 100L,
               device = "gpu")

# model summary
summary(model)

# coefficients
# coef(model)

# imp = importance(model)
# print(imp)
# 
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
tune_results = sjSDM_cv(Y = Y.train.pa_min10,
                        env = linear(data = env.vars, 
                                     formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo),
                        spatial = linear(xy.scale, ~0 + UTM_E:UTM_N),
                        biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
                        tune = "random", # random steps in tune-parameter space
                        learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
                        family = stats::binomial("probit"), # for both p/a and quasiprob data, default
                        CV = 5L, #  5L for 5-fold cross validation, nrow(as.matrix(otu.data)) for LOOCV
                        tune_steps = 60L, # 20L is default - is this the number of tunes for random search??? 
                        alpha_cov = seq(0, 1, 0.2), # species regularisation
                        alpha_coef = seq(0, 1, 0.2), # env reg
                        alpha_spatial = seq(0, 1, 0.2), # spatial reg
                        lambda_cov = seq(0, 0.1, 0.01), 
                        lambda_coef = seq(0, 0.1, 0.01),
                        lambda_spatial = seq(0,0.1,0.01),
                        device = "gpu",
                        n_cores = 2L, # 
                        n_gpu = 2L,# spread over this many GPU cards
                        iter = 100L, # 2L
                        sampling = 5000L # default is 5000L
)

# ncores
# NULL, # or 10L for small models, run this many sjsdm models at once (1 per CPU core). 10L is too much for large models because the 10 CPUs try to run on only the 2 GPUs available on ada: 5/GPU


# visualize tuning and best points:
pdf(file = "results_sjSDM/oregon_tune_trial.pdf")
plot(tune_results, perf = "logLik")
dev.off()
# green points in the best plot are (close to) the best lambda and alpha values

# save best lambdas and alphas
best <- plot(tune_results, perf = "logLik")

# save results
save(model, tune_results, best, file ="results_sjSDM/oregon_trial_tune.rdata")