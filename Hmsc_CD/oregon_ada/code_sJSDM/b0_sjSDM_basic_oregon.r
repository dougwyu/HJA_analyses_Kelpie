
## sJSDM on ADA

## trial run 
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() # always run sub from oregon_ada. # should be /gpfs/home/hsp20azu/oregon_ada
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")



# try basic model on oregon data
source("code_sjSDM/S1_read_data.r")
# gets env.vars, out.pa.csv, otu.qp.csv, S.train
rm(P, otu.qp.csv)

# Data reduction - reduce species
raretaxa <- which(colSums(otu.pa.csv > 0) < 10)
length(raretaxa)
otu.pa.minocc <- as.matrix(otu.pa.csv[, -raretaxa]) # reduced species
rm(raretaxa)


# scale all togher for now... just for test
env.vars.scale <- env.vars %>%
  mutate(across(where(is.numeric), scale))
head(env.vars.scale)

# XFormula1 <- as.formula(~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo)
# check names
# all(all.vars(XFormula1) %in% names(env.vars.scale))


# spatial data here:
# spatial data here:
head(Sp.data)
# scale spatial coords
Sp.data.scale <- scale(Sp.data[, c("UTM_E", "UTM_N")])
head(Sp.data.scale)
is.matrix(Sp.data.scale)

# this is what goes into the model
# mm <- model.matrix(XFormula1, data = env.vars.scale)
# head(mm)

# make and run model
model <- sjSDM_cv(Y = otu.pa.minocc,
               env = linear(data = env.vars.scale, 
                            formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht +
                              ht.r500 + cov4_16 + cov4_16.r500 + mTopo), # linear model on env covariates
               spatial = linear(data = Sp.data.scale, formula = ~0+UTM_E:UTM_N), # interactions of coordinates
               se = TRUE, family=binomial("probit"), sampling = 1000L, iter = 500L,
               device = "gpu")


# model summary
# summary(model)

# coefficients
# coef(model)


# model summary
summary(model)

# coefficients
coef(model)

imp = importance(model)
print(imp)

pdf("results_sjSDM/trial_oregon_sJSDM_importance_k40.pdf")
plot(imp)
dev.off()

an = anova(model)
print(an)

pdf("results_sjSDM/trial_oregon_sJSDM_anova_k40.pdf")
plot(an)
dev.off()


save(model, file ="results_sjSDM/oregon_trial_k40.rdata")


