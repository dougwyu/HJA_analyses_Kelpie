
## sJSDM on ADA

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
getwd() # where the .sub is run from (not where it is)


# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# trial
set.seed(42)

# simulate data
community <- simulate_SDM(sites = 100, species = 10, env = 3)
Env <- community$env_weights
Occ <- community$response
SP <- matrix(rnorm(200, 0, 0.3), 100, 2) # spatial coordinates (no effect on species occurences)

# head(community)
head(Env)
head(Occ)
head(SP)

Env <- as.data.frame(Env)
head(Env)

form1 <- as.formula(~V1+V2+V3)
form1

# data frame works within linear() but formula object doesn't/

# make and run model
model1 <- sjSDM(Y = Occ, 
               env = linear(data = Env, formula = ~V1+V2+V3), # # data frame works within linear()
               spatial = linear(data = SP, formula = ~0+X1:X2), # interactions of coordinates
               se = TRUE, family=binomial("probit"), sampling = 25L,
               device = "gpu")

summary(model1)

# try without linear()
model2 <- sjSDM(Y = Occ, 
               env = as.matrix(Env),   # just matrix of covariates..intercept included default, uses all covars,ie ~.
               spatial = linear(data = SP, formula = ~0+X1:X2), # interactions of coordinates
               se = TRUE, family=binomial("probit"), sampling = 25L,
               device = "gpu")

summary(model2)

# coefficients
coef(model1)

imp <- importance(model1)
print(imp)

pdf("results_sjSDM/trial_sJSDM_importance.pdf")
plot(imp)
dev.off()

an <- anova(model1)
print(an)

pdf("results_sjSDM/trial_sJSDM_anova.pdf")
plot(an)
dev.off()


save(model1, model2, file = "results_sjSDM/sjDM_trial_mods.rdata")