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

# # read in cross-validation output from sjsdm_cv_gpu_20200821.R
resultsfolder <- "results_20200824_275OTU_loocv"
tune_results <- readRDS(here(resultsfolder, 
                             "sjsdm_model_HJA_20200823.RDS"))
best = plot(tune_results, perf = "logLik")

# # run sjSDM model
# # lambda is the regularization strength
# # sjSDM supports l1 (lasso) and l2 (ridge) regularization:
  # # alpha weights the lasso or ridge penalty:
  # # - alpha = 0 --> pure lasso
  # # - alpha = 1.0 --> pure ridge
# green points in the best plot are (close to) the best lambda and alpha values
model <-  sjSDM(
  Y = as.matrix(otu.data),
  iter = 150L,
  learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
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
pdf(file = here(resultsfolder, "model_history.pdf"))
plot(model$history) # check iter
dev.off()

#save result
result = list(beta = coef(model),
              sigma = getCov(model),
              history = model$history,
              p = p,
              logLik=logLik(model)
              )

saveRDS(model, here(resultsfolder, "sjsdm_model_HJA_202080824.RDS"))
saveRDS(result, here(resultsfolder, "sjsdm_result_HJA_202080824.RDS"))
# model <- readRDS(here(resultsfolder, "sjsdm_model_HJA_202080824.RDS"))
# result <- readRDS(here(resultsfolder, "sjsdm_result_HJA_202080824.RDS"))

# #VP for test
imp <- importance(model)
pdf(file = here(resultsfolder, "importance.pdf"))
plot(imp)
dev.off()

# Rsquared2(model, individual = FALSE)
# Rsquared2(model)

an <- anova(model, cv = FALSE)
# # print(an)
pdf(file = here(resultsfolder, "anova.pdf"))
plot(an, percent = FALSE)
dev.off()

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
