
# Jan 18, 2021
# code on sjsdm model run with DNN (for ADA cluster) - from Yuanheng
# to test how long it would take to run one model on the cluster


options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
library(sjSDM)
packageVersion("sjSDM")

getwd() # should be /gpfs/home/hsp20azu/oregon_ada
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")

# one model to see how long it takes in ADA
load("data/yuanghen_mod_data.rdata")
# s.otu.train,scale.env.train, XY.train,  s.otu.test, scale.env.test, XY.test


# set variables
formula.env = 'envDNN'
lambda.env = .1
alpha.env = .5
lambda.sp = .1
alpha.sp = .9 
hidden = c(50L,50L,10L)
acti.sp = 'relu'
drop = .3

model.train = sjSDM(Y = s.otu.train,
                    env = DNN(data=scale.env.train, formula = ~.,
                              hidden=hidden, lambda = lambda.env, alpha = alpha.env, activation=acti.sp, dropout=drop, bias = T),
                    biotic = bioticStruct(lambda=lambda.env, alpha=alpha.env, on_diag=F, inverse = FALSE),
                    spatial = linear(data=XY.train, ~0+UTM_E*UTM_N, lambda=lambda.sp, alpha=alpha.sp),
                    learning_rate = 0.003, # 0.003 recommended for high species number 
                    step_size = NULL, iter = 150L, family=stats::binomial('probit'), sampling = 5000L 
)

names(model.train)

saveRDS(model.train, file = "results_sjSDM/ada_mod1.rds")



