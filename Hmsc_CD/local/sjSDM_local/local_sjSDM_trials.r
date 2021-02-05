# sjSDM::install_sjSDM(version = "cpu")

## sjSDM trial on laptop


library(sjSDM)
packageVersion("sjSDM")
# [1] ‘0.1.3.9000’
getwd() 

# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie")

# https://theoreticalecology.github.io/s-jSDM/
# https://github.com/TheoreticalEcology/s-jSDM

# try basic model on oregon data
source("Hmsc_CD/local/L1_read_data_v3.r")
rm(P, otu.sp.csv)

# Data reduction - reduce species to only those with more than minocc across all sites
minocc <- 25
raretaxa <- which(colSums(otu.pa.csv > 0) < minocc)
# no of species remaining
ncol(otu.pa.csv) - length(raretaxa)
otu.pa.minocc <- as.matrix(otu.pa.csv[, -raretaxa]) # reduced species
rm(raretaxa)

# sites without species?
sum(rowSums(otu.pa.minocc) == 0)

## scale data - all for test run
env.vars.scale <- env.vars %>%
  mutate(across(where(is.numeric), scale))

XFormula1 <- as.formula(~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + cov4_16 + cov4_16.r500 + mTopo)
# check names
all(all.vars(XFormula1) %in% names(env.vars))

# this is what goes into the model
# mm <- model.matrix(XFormula1, data = env.vars.scale)
# head(mm)

# spatial data here:
head(Sp.data)
# scale spatial coords
Sp.data.scale <- Sp.data %>%
  mutate(across(where(is.numeric), scale))
head(Sp.data.scale)


# make and run model
model <- sjSDM(Y = otu.pa.minocc,
               env = linear(data = env.vars.scale, 
                            formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + 
                              cov4_16 + cov4_16.r500 + mTopo), # linear model on env covariates
               spatial = linear(data = Sp.data.scale, formula = ~0+UTM_E:UTM_N), # interactions of coordinates
               se = TRUE, family=binomial("probit"), sampling = 100L, iter = 50L,
               device = "cpu")

# model summary
summary(model)

# coefficients
coef(model)

imp <- importance(model)
print(imp)

plot(imp)

an <- anova(model)
print(an)

plot(an)

save(model, an, imp, otu.pa.minocc, Sp.data, env.vars, file = "Hmsc_CD/local/sjSDM_local/trial29_res.rdata")


## tuning local#
tune_results = sjSDM_cv(Y = otu.pa.minocc,
                        env = linear(data = env.vars.scale, 
                                     formula = ~be10+B11_median+mean.EVI+insideHJA + Ess + ht + ht.r500 + 
                                       cov4_16 + cov4_16.r500 + mTopo),
                        spatial = linear(Sp.data.scale, ~0 + UTM_E:UTM_N),
                        biotic = bioticStruct(on_diag = FALSE, inverse = FALSE), # inverse=TRUE is 'better' but much slower
                        tune = "random", # random steps in tune-parameter space
                        learning_rate = 0.003, # 0.01 default, 0.003 recommended for high species number
                        family = stats::binomial("probit"), # for both p/a and quasiprob data, default
                        CV = 5L, #  5L for 5-fold cross validation, nrow(as.matrix(otu.data)) for LOOCV
                        tune_steps = 5L, # 20L is default - is this the number of tunes for random search??? 
                        alpha_cov = seq(0, 1, 0.2), # species regularisation
                        alpha_coef = seq(0, 1, 0.2), # env reg
                        alpha_spatial = seq(0, 1, 0.2), # spatial reg
                        lambda_cov = seq(0, 0.1, 0.01), 
                        lambda_coef = seq(0, 0.1, 0.01),
                        lambda_spatial = seq(0,0.1,0.01),
                        device = "cpu",
                        n_cores = 2L, # 
                        iter = 100L, # 2L
                        sampling = 1000L # default is 5000L
)

save(tune_results, file = "Hmsc_CD/local/sjSDM_local/trial29_tune.rdata")



load("Hmsc_CD/local/sjSDM_local/trial29_res.rdata")

mod.summ <- summary(model)

str(mod.summ, max.level = 1)

# coeffecients for all evn vars and species
mod.summ$coefs

mod.summ$coefmat[1:10, 1:4]
mod.summ$P[1:10, 1:10]

## extract ssignificant vars
29*11 # no sp * no coeff

sum(mod.summ$coefmat[, 4] < 0.05)

sig.vars <- rownames(mod.summ$coefmat)[mod.summ$coefmat[, 4] < 0.05]

var.df <- data.frame(rowN = sig.vars) %>%
  tidyr::separate(col = rowN, into = c("species", "predictor"),
                  remove = FALSE, sep = " ") %>%
  filter(predictor != "(Intercept)")

var.df
barplot(sort(table(var.df$predictor)))


