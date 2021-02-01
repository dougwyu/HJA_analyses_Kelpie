
wd <- here::here()
setwd(wd)
getwd()
# rm(wd)

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
Y.pa.min10 <- as.matrix(Y.train.pa[, -raretaxa]) # reduced species
rm(raretaxa)
# spatial data here:
head(S.train)
xy.scale <- scale(S.train[,c("UTM_E", "UTM_N")])
# head(xy.scale)

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


load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_eval.rdata")
str(model, max.level = 1)


# m1 <- matrix(c(1,5,3,7,3,4,5,5,7,8,8,10), ncol = 3)
# m2 <- matrix(c(2,0,3,0,3,4,10,5,10,5,7,10), ncol = 3)
# 
# m1
# m2
# 
# colSums(m1 * m2)

pred_data <- pred.link.test
obs_data <- Y.test.pa.min10

# x <- pred_data[,58]
# y <- obs_data[,58]
# diff(tapply(x, y, mean, na.rm = T))



prob.cm <- function(pred_data, obs_data, model, p, k, metrics){

  # p - number of predictors
  # k - no estimated parameters
  
  if(nrow(pred_data) != nrow(obs_data)) stop("")
  
  if(missing(metrics)) metrics <- c("TPR", "TPR", "TNR", "FNR", "Prev", "Acc", "Tjur", "AUC", "R2", "adjR2")
  
  TP <- colSums(pred_data * obs_data) # sum of probability of presences at True presences
  FP <- colSums(pred_data * (1 - obs_data)) # sum of prob of presences at True absences
  TN <- colSums((1 - pred_data) * (1 - obs_data)) # sum of probability of absence at true absences
  FN <- colSums((1 - pred_data) * obs_data) # sum of probability of absences at true presences
  
  TPR <- tryCatch(expr = TP / (TP + FN),
                  error = function(err){ return(NA) })                          # True Positive Rate / Sensitivity
  
  FPR <- tryCatch(expr = FP / (FP + TN),
                  error = function(err){ return(NA) })                         # False Positive Rate
  
  TNR <- tryCatch(expr = TN / (FP + TN),
                  error = function(err){ return(NA) })                         # True Negative Rate / Specificity
  
  FNR <- tryCatch(expr = FN / (TP + FN),
                  error = function(err){ return(NA) })                        # False Negative Rate
  
  Prevalence <- tryCatch(expr = (TP + FN) / (TP + FP + TN + FN),
                         error = function(err){ return(NA) })
  
  Accuracy <- tryCatch(expr = (TP + TN) / (TP + FP + TN + FN),
                       error = function(err){ return(NA) })
  
  
  Tjur <- tryCatch(expr = mapply(function(x,y) {
    tj <- base::diff(tapply(x, y, mean, na.rm = T)) # difference of average predicted values at 1 and 0
    if(!length(tj) > 0) tj <- NA
    return(tj)
    }, asplit(pred_data,2), asplit(obs_data,2)), error = function(err){ return(NA)})
  
  AUC <- tryCatch(expr = sapply(1:ncol(pred_data), function(i) Metrics::auc(obs_data[,i],pred_data[,i])),
                  error = function(err){ return(NA) })
  
  # lawson et al 2014 Appendix S9
  loglikP <- colSums( log( pred_data * obs_data + (1-pred_data)*(1-obs_data) ) )
  loglikN <- colSums( log( mean(pred_data)*obs_data + (1-mean(pred_data))*(1-obs_data) ) ) # NUll intercept only, (so mean)
  rsq <- (loglikP-loglikN)/(1-loglikN)
  
  # Guisan
  adjR2 <- 1-( (n-1)/(n-p) ) * (1-rsq)
  
  # Nakagawa 
  rsq <- 1 - (loglikN/loglikP)^(2/n)
  
  adjR2 <- rsq/(1 - loglikN^(2/n))
  
  # cm = cbind(TP=TP, FP=FP, TN=TN, FN=FN),
  res <- cbind(TPR=TPR, FPR=FPR, TNR=TNR, FNR=FNR, 
               Prev=Prevalence, Acc=Accuracy, 
               Tjur=Tjur, 
               AUC=AUC, 
               R2 = rsq, adjR2 = adjR2)
  
  # subset, if required
  res <- res[,metrics]
  
  return(res)
}

# probabilistic confusion matrix
cm <- prob.cm(pred.link.test, Y.test.pa.min10, p = 10)
head(cm)
apply(cm, 2, mean, na.rm = T)

auc.all <- prob.cm(pred.link.test, Y.test.pa.min10, metrics = "AUC")

# binary preditction with threshold for 'classic' confusion matrix
threshold <- 0.5
cm.bin <- prob.cm(pred.link.test>0.5, Y.test.pa.min10)
head(cm.bin)
apply(cm.bin, 2, mean, na.rm = T)


## AUC direction

lapply(1:5, function(i) pROC::roc(obs_data[,i], pred_data[,i], direction = "<", quiet = FALSE))
lapply(1:5, function(i) pROC::roc(obs_data[,i], pred_data[,i], direction = "<", quiet = FALSE)$auc)


tryCatch(expr = mapply(function(x,y) Metrics::auc(y,x),asplit(pred_data,2)[1:5], asplit(obs_data,2)[1:5]),
                error = function(err){ return(NA) })






