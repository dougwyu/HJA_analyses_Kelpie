
## get best run AUC to plot over species... 

options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada



pa <- read.csv("results_sjSDM/tuning_YL_pa_newVars/manual_tuning_sjsdm_newVars_S1_M1_pa_min5_envDNN_20210125.csv")

pa[which.max(pa$AUC.test),]

pa.best <- pa[which.max(pa$AUC.test),,drop = T]
pa.best

sprintf("%f10", seq(0,0.3, length.out=4))
sprintf("%f10",pa.best$lambda.env)

lambda.envN <- which(round(seq(0,0.3, length.out=4),5) == round(pa.best$lambda.env,5))
alpha.envN <- which(round(seq(.7,1, length.out=4),5) == round(pa.best$alpha.env,5))
lambda.spN <- which(round(seq(0,1, length.out=7),5) == round(pa.best$lambda.sp,5))
alpha.spN <- which(round(seq(0,1, length.out=7),5) == round(pa.best$alpha.sp,5))
hiddenN <- pa.best$hidden
dropN <- which(round(seq(.1,.5, length.out=3),5) == round(pa.best$drop,5))

abund = 'pa'		# 'qp','pa'  # is overwritten by data file loading.
minocc = 5
trap <- "M1"; period = "S1"
date.model.run = 20210125 
formula.env = 'envDNN'


best.mod <- glue::glue('s-jSDM_tuning_model__newVars_{period}_{trap}_{abund}_min{minocc}_{formula.env}_lambdaE{lambda.envN}_{alpha.envN}_{lambda.spN}_{alpha.spN}_hidden{hiddenN}_{dropN}.RDS')

best.mod

rm(lambda.envN, alpha.envN, lambda.spN, alpha.spN, hiddenN, dropN, trap, minocc, period, date.model.run, formula.env)

## load data 
load(paste0("data/yuanheng_mod_data_newVars_", abund, ".rdata"))

# load best model
mod.res <- readRDS(file.path("results_sjSDM/tuning_YL_pa_newVars", best.mod))

model.train <- mod.res$model

# best aucs

auc <- list()

for (pred in 1:2) {
  
  # 1 -> 'test'
  newdd = scale.env.test ; newsp = XY.test; otudd = s.otu.test
  
  if (pred==2) { newdd = NULL; newsp = NULL; otudd = s.otu.train}
  
  pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model.train, newdata=newdd, SP=newsp)),
                               along = -1L), 2:3, mean)
  
  attr(pred.dd, 'dimnames') = NULL

  # use all species, if no records, then NA is returned in the right place by tryCatch
  # if (pred==1) {
  #   otudd = rbind(otudd, count=(base::colSums(otudd)>0 & base::colSums(otudd)<dim(otudd)[1])*1 )
  #   pred.dd = pred.dd[ ,which(otudd[dim(otudd)[1],] == 1)]
  #   otudd = otudd[1:(dim(otudd)[1]-1), which(otudd[dim(otudd)[1],] == 1)]
  # }
  
 
   otudd.pa = (otudd>0)*1
 auc.dd = lapply(1:dim(otudd)[2], function(i) {
    
    tryCatch({
      pROC::roc(otudd.pa[,i], pred.dd[,i], direction = "<", quiet=T)$auc},
      error = function(err){ return(NA)}
      )
   }
    )
  
 auc[[pred]] <- auc.dd

  
  # if (pred==2) {tdd$AUC.explain=auc.mean}
  # if (pred==1) {tdd$AUC.test=auc.mean}
}


names(auc) <- c("test", "train")


save(auc, model.train, best.mod, s.otu.train, s.otu.test, file = "results_sjSDM/best_auc.rdata")