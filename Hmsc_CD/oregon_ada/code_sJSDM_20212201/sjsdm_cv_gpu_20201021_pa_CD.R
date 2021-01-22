# set some general options
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library("sjSDM")
packageVersion("sjSDM")


# # arguments
# args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
# print(args)
# # test if there is at least one argument: if not, return an error
# if (length(args) < 3) {
#   stop("At least three arguments must be supplied", call.=FALSE)
# }
# 

# packages
library(sjSDM)
library(here)
library(tidyverse)
library(fs)
library(glue)

# set variables
rundate <- 20201119 # run date
minocc <- 5 # minimum occupancy (incidence) per OTU
envvar <- "gismslidar" # gismslidarmin, gismslidar, gis, ms, lidar, mslidar

abund <- "pa" # pa is 0/1 data, qp is quasiprob data
# chk sjsdm_cv() code to see if have chosen DNN to fit the env covariates

resultsfolder <- glue("results_{rundate}_{minocc}minocc_{envvar}_{abund}_loocv")
dir_create(resultsfolder) # create results/ directory, but only if the results/ directory does not already exist
datafolder <- glue("data_{rundate}_{minocc}minocc_{envvar}")
  # paste0 alternative syntax
  # datafolder <- paste0("data_", rundate, "_", minocc, "minocc_", envvar)

# read in data
# env data:  scale.env, should be scaled
scale.env <- read_csv(here(datafolder, "scale.env.csv"))

# species data: otu.pa.csv, otu.qp.csv, qp == quasiprob, pa == 0/1
# check value of abund for qp or pa data
otu.data <- read_csv(here(datafolder, glue("otu.{abund}.csv")))

# XY data: XY.csv, should be scaled
XY <- read_csv(here(datafolder, "XY.csv"))

# sjSDM_cv
tune_results = sjSDM_cv(
  Y = as.matrix(otu.data),
  env = linear(as.matrix(scale.env)),
  # env = DNN(as.matrix(scale.env)),
  spatial = linear(XY, ~0 + UTM_E:UTM_N),
  biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
  tune = "random", # random steps in tune-parameter space
  learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
  family = stats::binomial("probit"), # for both p/a and quasiprob data, default
  CV = nrow(as.matrix(otu.data)), # 5L for 5-fold cross validation, nrow(as.matrix(otu.data)) for LOOCV
  tune_steps = 60L, # 20L is default
  alpha_cov = seq(0, 1, 0.1),
  alpha_coef = seq(0, 1, 0.1),
  lambda_cov = seq(0, 0.1, 0.001), 
  lambda_coef = seq(0, 0.1, 0.001),
  n_cores = 2L, # NULL, # or 10L for small models, run this many sjsdm models at once (1 per CPU core). 10L is too much for large models because the 10 CPUs try to run on only the 2 GPUs available on ada: 5/GPU
  n_gpu = 2L,# spread over this many GPU cards
  iter = 150L, # 2L
  sampling = 5000L # default is 5000L
  )
# https://rscs.uea.ac.uk/ada/using-ada/jobs/gpu
# each gpu node has 2 NVIDIA Quadro P5000 cards, each with 2560 GPU cores


# save cross-validation results
saveRDS(tune_results, here(resultsfolder, glue("sjsdm_tune_results_HJA_{rundate}.RDS")))

# visualize tuning and best points:
# pdf(file = here("results", "best_20200830.pdf"))
pdf(file = here(resultsfolder, glue("best_{rundate}.pdf")))
plot(tune_results, perf = "logLik")
dev.off()
# green points in the best plot are (close to) the best lambda and alpha values

# save best lambdas and alphas
best <- plot(tune_results, perf = "logLik")
saveRDS(best, here(resultsfolder, glue("sjsdm_tune_results_HJA_{rundate}_bestonly.RDS")))

# copy datafolder into results folder
dir_copy(datafolder, resultsfolder)

# END






# loocv with pa data commented out.  will be running either only qp or only pa data, which is set above in the abund variable
# ####
# # now run with pa data
# abund <- "pa" # "qp" # pa is 0/1 data, qp is quasiprob data
# 
# resultsfolder <- glue("results_{rundate}_{minocc}minocc_{envvar}_{abund}_loocv")
# dir_create(resultsfolder) # create results/ directory, but only if the results/ directory does not already exist
# 
# # species data:  otu.data
# otu.data <- read_csv(here(datafolder, glue("otu.data.{abund}.csv")))
# # or if i don't have a pa dataset,  i convert qp to 0/1 data
#   # otu.data.pa.csv <- otu.data.qp.csv
#   # otu.data.pa.csv[otu.data.pa.csv > 0] <- 1
#   # otu.data <- otu.data.pa.csv
# 
# # sjSDM_cv
# tune_results = sjSDM_cv(
#   Y = as.matrix(otu.data),
#   env = linear(as.matrix(scale.env)), 
#   spatial = linear(XY, ~0 + UTM_E:UTM_N),
#   biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
#   tune = "random", # random steps in tune-parameter space
#   learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
#   family = stats::binomial("probit"), # for both p/a and quasiprob data, default
#   CV = nrow(as.matrix(otu.data)), # 5L for 5-fold cross validation, nrow(Y) for LOOCV
#   tune_steps = 60L, # 20L is default
#   alpha_cov = seq(0, 1, 0.1),
#   alpha_coef = seq(0, 1, 0.1),
#   lambda_cov = seq(0, 0.1, 0.001), 
#   lambda_coef = seq(0, 0.1, 0.001),
#   n_cores = 2L, # NULL, # or 10L for small models, run this many sjsdm models at once (1 per CPU core). 10L is too much for large models because the 10 CPUs try to run on only the 2 GPUs available on ada: 5/GPU
#   n_gpu = 2L,# spread over this many GPU cards
#   iter = 150L, # 2L
#   sampling = 5000L # default is 5000L
# )
# # https://rscs.uea.ac.uk/ada/using-ada/jobs/gpu
# # each gpu node has 2 NVIDIA Quadro P5000 cards, each with 2560 GPU cores
# # link = "probit", # not used in sjSDM_cv
# 
# 
# # save cross-validation results
# saveRDS(tune_results, here(resultsfolder, glue("sjsdm_tune_results_HJA_{rundate}.RDS")))
# 
# # visualize tuning and best points:
# # pdf(file = here("results", "best_20200830.pdf"))
# pdf(file = here(resultsfolder, glue("best_{rundate}.pdf")))
# plot(tune_results, perf = "logLik")
# dev.off()
# # green points in the best plot are (close to) the best lambda and alpha values
# 
# # save best lambdas and alphas
# best <- plot(tune_results, perf = "logLik")
# saveRDS(best, here(resultsfolder, glue("sjsdm_tune_results_HJA_{rundate}_bestonly.RDS")))
# 
# # copy datafolder into results folder
# dir_copy(datafolder, resultsfolder)
# 



# # sjSDM test code
# community <- simulate_SDM(sites = 100, species = 10, env = 5)
# Env <- community$env_weights
# Occ <- community$response
# Sp <- matrix(rnorm(800), 100, 2) # spatial coordinates (no effect on species occurrences)
# 
# model <- sjSDM(Y = Occ, 
#                env = linear(data = Env, 
#                             formula = ~0 + X1*X2 + X3 + X4),
#                spatial = linear(data = Sp, 
#                                 formula = ~0 + X1:X2), 
#                se = TRUE,
#                device = "gpu"
# )
# summary(model)
# saveRDS(model, file = "model.RDS")
# Rsquared(model)
# 
# # model2 <- readRDS("model.RDS"); objects()
# # Rsquared(model2)
# 
# 
# imp <- importance(model)
# print(imp)
# pdf(file = "imp.pdf")
# plot(imp)
# dev.off()
# 
# 
# an = anova(model, cv = FALSE)
# print(an)
# pdf(file = "an.pdf")
# plot(an)
# dev.off()
