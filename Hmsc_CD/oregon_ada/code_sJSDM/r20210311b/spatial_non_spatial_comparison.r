#### Plot differences in spatial vs non spatial model ####


## RUNS ON LOCAL

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
wd <- here::here()
wd
setwd(file.path(wd, "Hmsc_CD"))
dir()

getwd()

resFolder <-"oregon_ada/code_sjSDM/r20210311b/results"
load(file.path(resFolder, "sp_results.rdata"))


# mean species evaluation metrics:
sp.mn.test
## here, 79 species - M1 and M2 shared

sapply(sp.mn.test, mean) # average across all mean metrics per species
# same as results here
apply(eval.results, 2, mean) # average across avereaged metrics for all species across 5 folds

## combine into single data frame 
spc.df <- cbind(stack(sp.mn.test), spatial = "yes")
head(spc.df)
rm(eval.results, sp.mn.train, sp.mn.test, sp.res.test, sp.res.train)

## load non spatial model results and combine
load(file.path(resFolder, "nsp_results.rdata"))

spc.df <- rbind(spc.df, cbind(stack(sp.mn.test), spatial = "no"))

head(spc.df)

unique(spc.df$ind)
table(spc.df[,c("ind", "spatial")])

library(ggplot2)


tapply(spc.df$values, spc.df$ind, range)


# spc.df %>%
#   group_by(ind)%>%
#   mutate(valN = (values-min(values))/(max(values)-min(values))) %>%
#   summarise(mx = max(valN),
#             min = min(valN))
# 
# spc.df %>%
#   group_by(ind)%>%
#   mutate(valN = (values-min(values))/(max(values)-min(values))) %>%
#   ggplot(aes(x = ind, y = valN, fill = spatial))+
#   geom_boxplot()+
#   ylab("Normalised metric values")


pdf("local/plots/spatial_vs_non_spatial_models_sjSDM_pa_M1M2_shared_79spp_vars2.pdf", width = 10, height = 6)
spc.df %>%
  ggplot(aes(x = ind, y = values, fill = spatial))+
  geom_boxplot()+
  facet_wrap(~ind, ncol = 6, scales = "free")
dev.off()