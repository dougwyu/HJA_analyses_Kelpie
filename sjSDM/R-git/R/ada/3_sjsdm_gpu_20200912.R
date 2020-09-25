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

# rm(list=ls())

# packages
library(sjSDM)
library(here)
library(tidyverse)
library(fs)
library(glue)
library(RColorBrewer)
dir_ls()

rundate <- 20200922 # sjsdm_cv run date
envvar <- "gismslidar" # gismslidar, mslidar, gis, ms, lidar
abund <- "qp" # "qp" # pa is 0/1 data, qp is quasiprob data

minocc <- 5 # minimum occupancy (incidence) per OTU, value from dataprep.Rmd

resultsfolder <- glue("results_{rundate}_{minocc}minocc_{envvar}_{abund}_loocv")
resultsfolder
datafolder <- glue("data_{rundate}_{minocc}minocc_{envvar}")
datafolder

# read in data
# env data:  scale.env1
scale.env1 <- read_csv(here(resultsfolder, datafolder, "scale.env1.csv"))

# species data:  otu.data
# comment in the dataset that i want to use. qp == quasiprob, pa == 0/1
otu.data <- read_csv(here(resultsfolder, datafolder, glue("otu.data.{abund}.csv")))

# XY data: XY
XY <- read_csv(here(resultsfolder, datafolder, "XY.csv"))


# read in cross-validation output from resultsfolder
best <- readRDS(here(resultsfolder, 
                     glue("sjsdm_tune_results_HJA_{rundate}_bestonly.RDS")))
best
# # lambda is the regularization strength
# # sjSDM supports l1 (lasso) and l2 (ridge) regularization:
# # alpha weights the lasso or ridge penalty:
# # - alpha = 0 --> pure lasso
# # - alpha = 1.0 --> pure ridge
# green points in the best plot are (close to) the best lambda and alpha values


# # run sjSDM model
model <-  sjSDM(
  Y = as.matrix(otu.data),
  iter = 150L,
  learning_rate = 0.003, # 0.01 default, 0l002 or 0.003 recommended for high species number, try a few values and choose the one with the smoothest model history
  family = stats::binomial("probit"), # for both p/a and quasiprob data, default
  env = linear(data = as.matrix(scale.env1),
               formula = ~.,
               # formula = ~ elevation.scale + canopy.ht.scale +
               #   min.T.scale + max.T.scale + precipitation.scale +
               #   metre.road.scale + metre.stream.scale +
               #   yrs.disturb.min.scale,
               lambda = best[["lambda_coef"]],
               alpha = best[["alpha_coef"]]
  ),
  spatial = linear(
    data = XY,
    formula = ~0 + UTM_E:UTM_N, # ~0 removes intercept from the spatial model
    lambda = best[["lambda_spatial"]],
    alpha = best[["alpha_spatial"]]
  ),
  biotic = bioticStruct(
    lambda = best[["lambda_cov"]],
    alpha = best[["alpha_cov"]],
    on_diag = FALSE,
    inverse = FALSE # inverse=TRUE is 'better' but much slower
  ),
  device = "gpu"
)


# summary(model)
# calculate post-hoc p-values:
p <- getSe(model)
summary.p <- summary(p)
# str(summary.p)

# model history
pdf(file = here(resultsfolder, glue("model_history_{rundate}.pdf")))
plot(model$history) # check iter
dev.off()

#save result
result = list(beta = coef(model),
              sigma = getCov(model),
              history = model$history,
              p = p,
              logLik=logLik(model)
              )

# VP for test
imp <- importance(model)
pdf(file = here(resultsfolder, glue("importance_{rundate}.pdf")))
plot(imp)
dev.off()


saveRDS(model, here(resultsfolder, glue("sjsdm_model_HJA_{rundate}.RDS")))
saveRDS(result, here(resultsfolder, glue("sjsdm_result_HJA_{rundate}.RDS")))
saveRDS(summary.p, here(resultsfolder, glue("sjsdm_summary.p_HJA_{rundate}.RDS")))
saveRDS(imp, here(resultsfolder, glue("sjsdm_imp_HJA_{rundate}.RDS")))
# model <- readRDS(here(resultsfolder, glue("sjsdm_model_HJA_{rundate}.RDS")))
# result <- readRDS(here(resultsfolder, glue("sjsdm_result_HJA_{rundate}.RDS")))
# summary.p <- readRDS(here(resultsfolder, glue("sjsdm_summary.p_HJA_{rundate}.RDS")))
# imp <- readRDS(here(resultsfolder, glue("sjsdm_imp_HJA_{rundate}.RDS")))


# Rsquared2(model, individual = FALSE)
# Rsquared2(model)

# does eight runs, each ~22 secs, total 3 mins
an <- anova(model, cv = FALSE)
pdf(file = here(resultsfolder, glue("anova_{rundate}.pdf")))
plot(an, percent = FALSE)
dev.off()
dir_ls(resultsfolder)



# END



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
