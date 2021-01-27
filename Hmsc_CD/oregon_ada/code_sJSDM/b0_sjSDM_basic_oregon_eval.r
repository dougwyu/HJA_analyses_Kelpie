
## sJSDM on ADA

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada




# local
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# try basic model on oregon data
source("code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train

# Data reduction - reduce species
raretaxa <- which(colSums(otu.pa.csv > 0) < 10)
length(raretaxa)
Y.pa.min10 <- as.matrix(otu.pa.csv[, -raretaxa]) # reduced species
rm(raretaxa)

XFormula1 <- as.formula(~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo)
# check names
all(all.vars(XFormula1) %in% names(env.vars))
# str(X.train)


# spatial data here:
head(S.train)
xy.scale <- scale(S.train[,c("UTM_E", "UTM_N")])
# head(xy.scale)

# this is what goes into the model
# mm <- model.matrix(XFormula1, data = env.vars)
# head(mm)

# do testing and training data
test.size <- floor(nrow(Y.pa.min10)*0.2)
train.size <- nrow(Y.pa.min10) - test.size

set.seed(99)
test.rows <- sample(1:nrow(Y.pa.min10), size = test.size)
train.rows <- setdiff(1:nrow(Y.pa.min10), test.rows)

# sort(c(test.rows, train.rows))

Y.test.pa.min10 <- Y.pa.min10[test.rows, ]
Y.train.pa.min10 <- Y.pa.min10[train.rows, ]

X.test <- env.vars[test.rows,]
X.train <- env.vars[train.rows,]

xy.train <- xy.scale[train.rows,]
xy.test <- xy.scale[test.rows,]

# make and run model
model <- sjSDM(
  Y = Y.train.pa.min10,
  env = linear(data = X.train,
               # linear model on env covariates
               formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo),
  spatial = linear(data = xy.train, formula = ~0+UTM_E:UTM_N), # interactions of coordinates
  se = TRUE, family=binomial("probit"), 
  sampling = 100L,
  device = "gpu")

# model summary
summary(model)

## S3 method for class 'sjSDM'
pred.link.train <- predict(model, newdata = X.train, SP = xy.train, type = "link")
pred.raw.train <- predict(model, newdata = X.train, SP = xy.train, type = "raw")

pred.link.test <- predict(model, newdata = X.test, SP = xy.test, type = "link")
pred.raw.test <- predict(model, newdata = X.test, SP = xy.test, type = "raw")


## testing AUC - Yuanghen code.
for (pred in 1:2) {
  
  # 1 -> 'test'
  newdd = X.test; newsp = xy.test; otudd = Y.test.pa.min10
  
  if (pred==2) { newdd = NULL; newsp = NULL; otudd = Y.train.pa.min10}
  
  predL <- lapply(1:3, function(i) predict(model, newdata=newdd, SP=newsp))
  pred.dd = apply(abind::abind(predL, along = -1L), 2:3, mean)
  attr(pred.dd, 'dimnames') = NULL
  
  if (pred==1){
    otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
    pred.dd = pred.dd[ , which(otudd[dim(otudd)[1],] == 1)]
    otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
  }
  
  otudd.pa = (otudd>0)*1
  roc.dd = lapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i]))
  auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) stringr::str_split(pROC::auc(i), ':')), function(x) x[[1]][1] )))
  
  if (pred==2) {AUC.explain <- auc.mean}
  if (pred==1) {AUC.test <- auc.mean}
}


save(model, 
     pred.link.train, pred.link.test, pred.raw.train, pred.raw.test,
     predL, AUC.explain, AUC.test,
     file ="results_sjSDM/oregon_eval.rdata")



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

# sjSDM_cv
# 
# tune_results = sjSDM_cv(Y = Y.train.pa_min10,
#                         env = linear(data = env.vars, 
#                                      formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo),
#                         spatial = linear(xy.scale, ~0 + UTM_E:UTM_N),
#                         biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
#                         tune = "random", # random steps in tune-parameter space
#                         learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
#                         family = stats::binomial("probit"), # for both p/a and quasiprob data, default
#                         CV = 5L, #  5L for 5-fold cross validation, nrow(as.matrix(otu.data)) for LOOCV
#                         tune_steps = 60L, # 20L is default - is this the number of tunes for random search??? 
#                         alpha_cov = seq(0, 1, 0.2), # species regularisation
#                         alpha_coef = seq(0, 1, 0.2), # env reg
#                         alpha_spatial = seq(0, 1, 0.2), # spatial reg
#                         lambda_cov = seq(0, 0.1, 0.01), 
#                         lambda_coef = seq(0, 0.1, 0.01),
#                         lambda_spatial = seq(0,0.1,0.01),
#                         device = "gpu",
#                         n_cores = 2L, # 
#                         n_gpu = 2L,# spread over this many GPU cards
#                         iter = 100L, # 2L
#                         sampling = 5000L # default is 5000L
# )
# 
# # ncores
# # NULL, # or 10L for small models, run this many sjsdm models at once (1 per CPU core). 10L is too much for large models because the 10 CPUs try to run on only the 2 GPUs available on ada: 5/GPU
# 
# 
# # visualize tuning and best points:
# pdf(file = "results_sjSDM/oregon_tune_trial.pdf")
# plot(tune_results, perf = "logLik")
# dev.off()
# # green points in the best plot are (close to) the best lambda and alpha values
# 
# # save best lambdas and alphas
# best <- plot(tune_results, perf = "logLik")
# 
# # save results
# save(model, tune_results, best, file ="results_sjSDM/oregon_trial_tune.rdata")