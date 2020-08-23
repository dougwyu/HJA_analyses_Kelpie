# set some general options
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/b042/sjSDM_env/bin/python")

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

# read in data
# env data:  scale.env1
scale.env1 <- read_csv(here("data", "scale.env1.csv"))

# species data:  otu.data
otu.data <- read_csv(here("data", "otu.data.csv"))

# XY data: XY
XY <- read_csv(here("data", "XY.csv"))

# sjSDM_cv
tune_results = sjSDM_cv(
  Y = as.matrix(otu.data),
  env = linear(as.matrix(scale.env1)), 
  spatial = linear(XY, ~ 0 + UTM_E:UTM_N),
  biotic = bioticStruct(on_diag=FALSE), 
  tune = "random", # random steps in tune-parameter space
  CV = nrow(as.matrix(otu.data)), # 5L for 5-fold cross validation, nrow(Y) for LOOCV
  tune_steps = 60L, # 20L is default
  alpha_cov = seq(0, 1, 0.1),
  alpha_coef = seq(0, 1, 0.1),
  lambda_cov = seq(0, 0.1, 0.001), 
  lambda_coef = seq(0, 0.1, 0.001),
  n_cores = 2L, # NULL, # or 10L for small models, meaning that 10 sjSDM jobs are started (1 per CPU core), and distributed evenly over the GPUs: 5/GPU
  n_gpu = 2,# ada GPU queues have 2 GPUs
  iter = 150L, # 2L
  sampling = 5000L # default is 5000L
  )
# https://rscs.uea.ac.uk/ada/using-ada/jobs/gpu
# each gpu node has 2 NVIDIA Quadro P5000 cards, each with 2560 GPU cores
# link = "probit", # not used in sjSDM_cv


# save cross-validation results
saveRDS(tune_results, here("results", "sjsdm_tune_results_HJA_20200823.RDS"))

# visualize tuning and best points:
pdf(file = here("results", "best_20200823.pdf"))
plot(tune_results, perf = "logLik")
dev.off()


# # read in cross-validation results
# tune_results <- readRDS(here::here("results", "sjsdm_tune_results_HJA_20200823.RDS"))
# best = plot(tune_results, perf = "logLik") # perf = "AUC"

# # run sjSDM model
# model <-  sjSDM(
#   Y = as.matrix(otu.data),
#   iter = 100L,
#   learning_rate = 0.01, # 0.01 default, but 0.003 recommended for high species number
#   link = "probit", # for both p/a and quasiprob data
#   env = linear(data = as.matrix(scale.env1),
#                formula = ~.,
#                # formula = ~ elevation.scale + canopy.ht.scale +
#                #   min.T.scale + max.T.scale + precipitation.scale +
#                #   metre.road.scale + metre.stream.scale +
#                #   yrs.disturb.min.scale,
#                lambda = best[["lambda_coef"]],
#                alpha = best[["alpha_coef"]]
#   ),
#   spatial = linear(
#     data = XY,
#     formula = ~0 + UTM_E:UTM_N, # ~0 removes intercept from the spatial model
#     lambda = best[["lambda_spatial"]],
#     alpha = best[["alpha_spatial"]]
#   ),
#   biotic = bioticStruct(
#     lambda = best[["lambda_cov"]],
#     alpha = best[["alpha_cov"]],
#     on_diag = FALSE
#   )
# )


# # summary(model)
# # calculate post-hoc p-values:
# p = getSe(model)
# (summary.p=summary(p))
# plot(model$history) # check iter
# 
# #save result
# result = list(beta = coef(model), 
#               sigma = getCov(model), 
#               history = model$history, 
#               p = p, 
#               logLik=logLik(model)
# )

# saveRDS(model, here("results", "sjsdm_model_HJA_20200608.RDS"))
# saveRDS(result, here("results", "sjsdm_results_HJA_20200608.RDS"))
# result <- readRDS(here::here("results", "sjsdm_results_HJA_20200608.RDS"))
# model <- readRDS(here::here("results", "sjsdm_model_HJA_20200608.RDS"))

# #VP for test
# imp = importance(model)
# plot(imp)
# Rsquared2(model, individual = FALSE)
# # Rsquared2(model)
# 
# an <- anova(model, cv = FALSE)
# # print(an)
# plot(an, percent = FALSE)
# plot(an)

# an_cv <- anova(model, cv = 5L)
# # print(an)
# plot(an_cv)
# plot(an_cv, percent = FALSE)



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
