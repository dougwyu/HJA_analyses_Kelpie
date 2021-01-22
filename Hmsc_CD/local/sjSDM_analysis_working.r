

wd <- here::here()
setwd(wd)
getwd()
rm(wd)

# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# get model input data
source("Hmsc_CD/oregon_ada/code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train

# Data reduction - reduce species
raretaxa <- which(colSums(Y.train.pa > 0) < 10)
length(raretaxa)
Y.train.pa_min10 <- as.matrix(Y.train.pa[, -raretaxa]) # reduced species
rm(raretaxa)



load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_trial.rdata")
load("results_sjSDM/oregon_trial.rdata")
str(model, max.level = 1)

head(link_pred)
summary(link_pred)

summary(raw_pred)

link_pred[1:10,1:10]

(binomial(link = "probit")$linkinv(raw_pred))[1:10,1:10]


# pred data is link_pred
# obs data is Y.train.pa_min10


m1 <- matrix(c(1,5,3,7,3,4,5,5,7,8,8,10), ncol = 3)
m2 <- matrix(c(2,0,3,0,3,4,10,5,10,5,7,10), ncol = 3)

m1
m2

colSums(m1 * m2)

prob.cm <- function(pred_data, obs_data){
  
  TP <- colSums(pred_data * obs_data)
  FP <- colSums(pred_data * (1 - obs_data))
  TN <- colSums((1 - pred_data) * (1 - obs_data))
  FN <- colSums((1 - pred_data) * obs_data)
  
  return(cbind(TP=TP, FP=FP, TN=TN, FN=FN))
  
}

cm <- prob.cm(link_pred, Y.train.pa_min10)
head(cm)

cm.bin <- prob.cm(link_pred>0, Y.train.pa_min10)

head(cm.bin)

colSums(link_pred>0)


#### Tuning data

load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_trial_tune.rdata")

str(model, max.level = 1)
str(tune_results, max.level = 1)
head(tune_results$short_summary) # 60 steps is 60 iterations in full summary...  so 60 random choices of tuning grid. 
# average AUC, etc.


# tune_results - list of 60. one for eah step

# each step has 5 cross valication runs
str(tune_results$tune_results[[1]], max.level = 1) # cross validation results


str(tune_results$tune_results[[1]][[1]], max.level = 1) # cross validation results

head(tune_results$summary) # here 60 iterations, times 5 for cross validation = 60*5
summary(tune_results$summary) 

tune_results$settings

max(tune_results$short_summary$AUC_test)
tune_results$short_summary$AUC_test

summary(tune_results)

