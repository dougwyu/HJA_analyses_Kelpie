### AUC functions


# testing data: 17 sites (20%) with response and predicted values from pa cross validation, 4 species
Y.test <- structure(c(1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
                      1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 
                      0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 
                      0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(17L, 4L))

Pred.test <- structure(c(0.93892, 0.66146, 0.47845, 0.84288, 0.3483, 0.415, 
                         0.97109, 0.91469, 0.86533, 0.76702, 0.6164, 0.91962, 0.78906, 
                         0.91927, 0.73108, 0.83802, 0.92296, 0.60071, 0.13734, 0.29069, 
                         0.35382, 0.11483, 0.17761, 0.3793, 0.15075, 0.41787, 0.33493, 
                         0.24288, 0.14474, 0.19351, 0.33881, 0.13977, 0.11995, 0.12558, 
                         0.53059, 0.85161, 0.84674, 0.95095, 0.81294, 0.99165, 0.34469, 
                         0.85878, 0.65415, 0.42095, 0.92685, 0.89672, 0.88219, 0.85416, 
                         0.97768, 0.97595, 0.91641, 0.34734, 0.18854, 0.10646, 0.27799, 
                         0.28179, 0.48508, 0.11032, 0.24834, 0.14395, 0.18885, 0.21009, 
                         0.18454, 0.16242, 0.21792, 0.26658, 0.19304, 0.14045), 
                       .Dim = c(17L, 4L), .Dimnames = list(NULL, NULL))

Y.test
Pred.test

# calculate AUC values - pROC
sapply(1:ncol(Pred.test), function(i) as.numeric(pROC::roc(Y.test[,i],Pred.test[,i])$auc))
# Metrics package
sapply(1:ncol(Pred.test), function(i) Metrics::auc(Y.test[,i],Pred.test[,i]))
# second value is different

# pROC::auc() requires direction to be specified
rocList <- lapply(1:ncol(Pred.test), function(i) pROC::roc(Y.test[,i],Pred.test[,i], direction = "<"))

## In the second species, the average predicted values at presences happen to be lower than the average values for absences
## Can see this as the equivalent of the Wilcoxon test - see plots

par(mfcol = c(2,2))

for(i in 1:4){
  
  plot(rocList[[i]])
  text(1, 0.8, paste("direction:", rocList[[i]]["direction"]))
  boxplot(Pred.test[,i] ~ Y.test[,i])
  
}


# specify direction in pROC - see 'direction' in helpfile
# "<": if the predictor values for the control [0] group are lower or equal than the values of the case [1] group  
sapply(1:ncol(Pred.test), function(i) as.numeric(pROC::roc(Y.test[,i],Pred.test[,i], direction = "<")$auc))

# Metrics package
# from help file: where the smallest values correspond to the observations most believed to be in the negative class and the
# largest values indicate the observations most believed to be in the positive class
sapply(1:ncol(Pred.test), function(i) Metrics::auc(Y.test[,i],Pred.test[,i]))



