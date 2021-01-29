

wd <- here::here()
setwd(wd)
getwd()
rm(wd)

# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# get model input data
setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
source("code_sjSDM/S1_read_data.r")
setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")

# gets env.vars, out.pa.csv, otu.qp.csv, S.train

# Data reduction - reduce species
raretaxa <- which(colSums(Y.train.pa > 0) < 10)
length(raretaxa)
Y.train.pa.min10 <- as.matrix(Y.train.pa[, -raretaxa]) # reduced species
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


# tuning results - best tune: based on min logLik
best

# tune results short summary
head(tune_results$short_summary) # 60 steps is 60 iterations in full summary...  so 60 random choices of tuning grid. 
# average AUC, etc.
tail(tune_results$short_summary)


# best test AUC
head(tune_results$short_summary[order(tune_results$short_summary$AUC_test, decreasing = T),])
max(tune_results$short_summary$AUC_test)

# best loglik
head(tune_results$short_summary[order(tune_results$short_summary$logLik),])


# tune_results - list of 60. one for each step
# each step has 5 cross valication runs
str(tune_results$tune_results[[1]], max.level = 1) # cross validation results

str(tune_results$tune_results[[1]][[1]], max.level = 1) # cross validation results
str(tune_results$tune_results[[1]][[2]], max.level = 1) # cross validation results

# should be same parameter from tuning grid
tune_results$tune_results[[1]][[1]]$pars
tune_results$tune_results[[1]][[2]]$pars

head(tune_results$summary) # here 60 iterations, times 5 for cross validation = 60*5
summary(tune_results$summary) 

tune_results$settings # tune_samples is the tuning grid


## if you want to get species scores... need to go into CV results
str(tune_results$tune_results[[1]][[1]], max.level = 1) # cross validation results

# parameters from tuning grid
tune_results$tune_results[[1]][[1]]$pars

# test and train predictions here
pred.test <- tune_results$tune_results[[1]][[1]]$pred_test
tune_results$tune_results[[1]][[1]]$pred_train

# indices of data
tune_results$tune_results[[1]][[1]]$indices

# so training response will be:
Y.test <- Y.train.pa.min10[tune_results$tune_results[[1]][[1]]$indices, ]
Y.train <- Y.train.pa.min10[-tune_results$tune_results[[1]][[1]]$indices, ]

pROC::roc(Y.test[,1], tune_results$tune_results[[1]][[1]]$pred_test[,1])

# remove empty sites
Y.test = rbind(Y.test, count=(base::colSums(Y.test)>0 & base::colSums(Y.test)<dim(Y.test)[1])*1 )
pred.test = pred.test[ , which(Y.test[dim(Y.test)[1],] == 1)]
Y.test <- Y.test[1:(dim(Y.test)[1]-1), which(Y.test[dim(Y.test)[1],] == 1)]
mean(sapply(1:dim(Y.test)[2], function(i) pROC::roc(Y.test[,i], pred.test[,i], quiet = TRUE)$auc))

mean(sapply(1:dim(Y.test)[2], function(i) pROC::roc(Y.test[,i], pred.test[,i], quiet = TRUE, direction = "<")$auc))

sapply(1:dim(Y.test)[2], function(i) pROC::roc(Y.test[,i], pred.test[,i], quiet = TRUE)$auc)
sapply(1:dim(Y.test)[2], function(i) pROC::roc(Y.test[,i], pred.test[,i], quiet = TRUE, direction = "<")$auc)

# what are predictions for species with all 0 at sites?
# ind <- which(colSums(Y.test)==0)
# pred.test[, ind]

Pred.test <- pred.test
save(Y.test, Pred.test, file = "Hmsc_CD/local/test_example.rdata")
Y.test <- unname(Y.test[,1:4])
row.names(Y.test) <- NULL
Pred.test[,1:5]

dput(Y.test)
dput(round(Pred.test[,1:4],5))

#######  evaluation 2nd run working #####
## from b0_sjSDM_basic_oregon_eval.r
load("results_sjSDM/oregon_eval.rdata")

# # first, how long does it take?  ## get input data from above ....eval.r --  needs pytorch - better on gpu
# system.time(
#   tmp <- predict(model, newdata = X.train, SP = xy.train, type = "link")
# )
# 

single.a <- abind::abind(predL, along = -1L) # bind to array. no predicts: no sites: no species

# get mean of all predictions for each species, each site
pred.dd = apply(single.a, 2:3, mean)
attr(pred.dd, 'dimnames') = NULL

otudd = Y.test.pa.min10
pred.dd <- pred.link.test
  
otudd = obs_data
pred.dd <- pred_data


if (pred==1){ # testing, otud is test - remove species with all 0 counts at site
    otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
    pred.dd = pred.dd[ , which(otudd[dim(otudd)[1],] == 1)]
    otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
  }
  
  otudd.pa = (otudd>0)*1
  # do roc for each species i = 1
  
  tmp <- pROC::roc(otudd.pa[,i], pred.dd[,i])
  table(otudd.pa[,i])
  pROC::auc(tmp)
  str(tmp)
  tmp$auc
  
  lapply(1:5, function(i) pROC::roc(otudd.pa[,i], pred.dd[,i], quiet = FALSE))
  
  roc.dd = lapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i], quiet = FALSE))
  
  # auc for all species
  auc.all <- sapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i], quiet = TRUE)$auc)
  
  roc.mn <- mean(auc.all)
  
  # mean of all species
  auc.mean = mean(as.numeric(sapply(lapply(roc.dd, function(i) stringr::str_split(pROC::auc(i), ':')), function(x) x[[1]][1] )))

auc.mean == roc.mn  
  
if (pred==2) {AUC.explain <- auc.mean}
if (pred==1) {AUC.test <- auc.mean}


## Single predictrion results

table(Y.test.pa.min10[,i])

otudd = Y.test.pa.min10
pred.dd <- pred.link.test

otudd = Y.train.pa.min10
pred.dd <- pred.link.train


otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
pred.dd = pred.dd[ , which(otudd[dim(otudd)[1],] == 1)]
otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]

otudd.pa = (otudd>0)*1

auc.single.pred <- sapply(1:dim(otudd)[2], function(i) pROC::roc(otudd.pa[,i], pred.dd[,i], quiet = TRUE)$auc)
  
singl.mn <- mean(auc.single.pred)
