

wd <- here::here()
setwd(wd)
getwd()
rm(wd)


#### Tuning data
load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_trial_cf_vars.rdata")

str(summ.old, max.level = 1)

summ <- summ.new
summ <- summ.old

str(summ.new$coefs, max.level = 1)

str(summ.new$coefs$env, max.level = 1)

summ.new$coefs$env[[7]][1:5, 1:5]

var.sig <- function(summ, p = 0.05, ...){
  
  nSig <- sum(summ$coefmat[, 4] < p, na.rm = T)
  
  print(nSig)
  
  if(nSig >0) {
    
  sig.vars <- rownames(summ$coefmat)[summ$coefmat[, 4] < 0.05]
  
  var.df <- data.frame(rowN = sig.vars) %>%
    tidyr::separate(col = rowN, into = c("species", "predictor"),
                    remove = FALSE, sep = " ") %>%
    filter(predictor != "(Intercept)")
  return(var.df)
  
  barplot(sort(table(var.df$predictor)))

  }
}

old.df <-var.sig(summ.old)
new.df <- var.sig(summ.new)

png("Hmsc_CD/local/plots/cf_vars.png")
par(mfrow = c(1,2))
var.sig(summ.old, main = "Old")
var.sig(summ.new, main = "new")
dev.off()





# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")

#### Tuning data
load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_trial_tune_cv2.rdata")

str(tune_results_new_vars, max.level = 1)

# tune results short summary
head(tune_results_new_vars$short_summary) # 60 steps is 60 iterations in full summary...  so 60 random choices of tuning grid. 

# tune_results
# each step has 5 cross valication runs
str(tune_results_new_vars$tune_results[[1]], max.level = 1) # cross validation results
str(tune_results_new_vars$tune_results[[1]][[1]], max.level = 1) # cross validation results

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
