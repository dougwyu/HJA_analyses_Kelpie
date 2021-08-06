#### ecocopula analysis

## copied initially from oregon_ada/code/P1_ecocopula_analysis_v3.r

options(echo=TRUE) # if you want see commands in output file

# Qsne 2020-01-24 is now installed on ADA as an environment module
# It can be accessed by executing the following command
# # module add qsne/2020-01-24
# 
# For accessing in R 3.6.2  you have to run the following command inside the R
# 
# module add R/3.6.2
# R
source('code_sjSDM/qsne.r')

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
# setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")


baseFolder <- "code_sjSDM/r20210716b"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"


## load species AUC resutls for filtering
load(file.path(resFolder, "sp_test_results.rdata")) # # eval.results, sp.res.test, sp.res.train

# load model data - for filename parameters
load(file.path(resFolder, paste0("modelData_",abund,".rdata")))
rm(env.vars, k, noSteps, vars, device, iter, sampling, otuenv)
# otu.pa.csv, otu.qp.csv

## Mean AUC per species (and other eval metrics)
str(sp.res.test, max.level = 1)
head(sp.res.test$auc)

## Filter species by auc
auc.filt <- 0.70
# threshold for presence absence data
# tr <- 0.5

# how many species after AUC filter?
sum(sp.res.test$auc > auc.filt, na.rm = T)

# extrapolated predictions
# load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, ".rdata")))
# pred.mn, pred.sd, 

# clamp predictions
load(file.path(resFolder, paste0("sjSDM_predictions_", "M1S1_", "min", minocc, "_", varsName, "_", abund, "_clamp", ".rdata")))
# pred.mn.cl, pred.sd.cl

dim(pred.mn.cl)

## filter for species performance
# pred.in <- pred.mn[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
# dim(pred.in)

# clamp pred
# pred.in.cl <- pred.mn.cl[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]

## Full data set
Xmat <- pred.mn.cl[,sp.res.test$auc > auc.filt & !is.na(sp.res.test$auc)]
# [1] 247743  rows


# half data set by pca
# x.rda <- vegan::rda(pred.in.cl, scale = TRUE)
# Xmat <- x.rda$Ybar[,1:20]

## reduced data set trial
# Xmat <- pred.in.cl[1:10000,]

# pa version
# Xmat <- (pred.mod[indNA2, ] >= tr)*1

dim(Xmat)
#Xmat[1:10, 1:10]

# perplexity <- 30
# 3 * perplexity < nrow(Xmat) - 1
# 
# N <- nrow(Xmat)
# N^(1/2)


# Max
(nrow(Xmat) - 1)/3


## qsne function
# qsne <- function(X, dims=2, perplexity=0, hess_rank=10, kld_ival=10, pca=F, max_iter=1000, tol=1e-6, Y_init=NA, verbose=F)

# qsne_res <- qsne(Xmat, dims = 2, perplexity = perplexity, max_iter=1000)

dims <- 2
perplexity <- 50
hess_rank <- 10
kld_ival <- 10
max_iter <- 1000
tol <- 0.0001
nCores <- parallel::detectCores()-1
pY <- file.path(resFolder, "qsne_out.txt")

pX <- tempfile("X")
write.table(Xmat, pX, sep="\t", row.names=F, col.names=F)

cmd <- sprintf("/gpfs/software/ada/qsne/2020-01-24/bin/qsne -d %d -p %s -m %d -K %d -k %d -T %d -o %s %s", dims, perplexity, hess_rank, kld_ival, max_iter, nCores, pY, pX) # pYi pE  pO
# run
system(cmd)

# system("/gpfs/software/ada/qsne/2020-01-24/bin/qsne -d2 -p30 -o code_sjSDM/qsne_out.txt -O code_sjSDM/qsne_obj.txt -S code_sjSDM/qsne_ent.txt code_sjSDM/iris.txt")

# save(qsne_res, file = file.path(resFolder, "qsne_results_cl_pca20.rdata"))


## again with perplexity tuning.
# perplexity <- "30:300"
# 
# qsne_p_300 <- qsne(Xmat, dims = 2, perplexity = perplexity, max_iter=1000)
# save(qsne_p_300, file = file.path(resFolder, "qsne_results_cl_p300_pca20.rdata"))