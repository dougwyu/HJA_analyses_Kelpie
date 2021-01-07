#### ecocopula analysis

getwd()

wd <- here::here()
setwd(wd)

library(dplyr)
library(mvabund)
library(ecoCopula)
library(here)
library(glue)
library(vegan)
library(ggplot2)
library(RColorBrewer)

## Get data from 
# "Hmsc_CD/oregon_ada/code/S1_read_data_cd.r" 
load("Hmsc_CD/oregon_ada/data/allData_pre_selection.rdata")

env.csv <- env.vars %>%
  select(clearcut, insideHJA, oldGrowthIndex, elevation_m, dem500, tri.pt, canopyHeight_m, precipitation_mm, mean.NDVI, mean.EVI, mean.green, mean.wet, l_p25, l_rumple, B1_median, B4_median, lg_DistStream, lg_DistRoad, lg_YrsDisturb, lg_cover2m_max, lg_cover2m_4m, lg_cover4m_16m)

source("Hmsc_CD/local/ecoCopula_plot_fn.r")

## do the models
# mvabund models, mod0 and mod1
# abund <- "pa"
# glue("otu.{abund}.csv")

# make a pa matrix as mvabund object
otu.pa <- mvabund(otu.pa.csv)
otu.qp <- mvabund(otu.qp.csv)

is.mvabund(otu.pa)
is.mvabund(otu.qp)
str(otu.pa)

# with no env covariates - PA data, no offset. (Only for raw, ab data)
mod0 <- manyglm(otu.pa ~ 1, 
                    family = "negative.binomial", # family = binomial("cloglog") 
                    data = env.csv) 

plot(mod0) # chk residuals

# do ordination
mod0.ord <- cord(mod0)
str(mod0.ord, max.level = 1)

# Site scores: join ordination axes to data frame of predictors, coordinates
site_res <- data.frame(mod0.ord$scores, env.csv, XY.csv)
str(site_res, max.level=1)

# Species scores
sp_res <- data.frame(mod0.ord$loadings, species = colnames(otu.pa))
str(sp_res)

cor.preds <- colnames(env.csv)[sapply(env.csv, is.numeric)]

(mod0.cor <- cor(mod0.ord$scores, site_res[,cor.preds]))

alpha <- 4*0.95

# plot factor 1 and 2
ggplot() + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = Factor1 * alpha, 
                   yend = Factor2 * alpha), 
               data = sp_res, 
               size = .1) +
  geom_point(aes(x = Factor1, y = Factor2, 
                 color = dem500, 
                 size = exp(lg_YrsDisturb), 
                 shape = as.factor(insideHJA)),
             data = site_res) + 
  scale_color_gradientn(colours = brewer.pal(n = 10, name = "RdYlBu")) +
  theme_classic() + # or "PuOr" or "RdYlBu"
  labs(title = "null model pa") +
  xlab("Factor 1") +
  ylab("Factor 2")


# Plot over coords
ggplot() + 
    geom_point(aes(x = UTM_E, y = UTM_N, 
                   color = Factor1, 
                   size = exp(lg_YrsDisturb), # YrsSinceDist, l_rumple
                   shape = as.factor(insideHJA)
    ), 
    data = site_res) +
    scale_color_gradientn(colours = brewer.pal(n = 10, name = "PuOr")) +
    theme_classic() + # or "PuOr" or "RdYlBu"
    labs(title = glue("HJA, XY plot, Factor1, pa")) +
    xlab("UTM_E") +
    ylab("UTM_N")

# factor 1 mainly inside HJA

mod01 <- manyglm(otu.pa ~ dem500, 
                family = "negative.binomial", # family = "binomial"
                data = env.csv) # , show.time = "all" # ??

mod01.ord <- cord(mod01)
site_res01 <- data.frame(mod01.ord$scores, env.csv, XY.csv)
sp_res01 <- data.frame(mod01.ord$loadings, species = colnames(otu.pa))
(mod01.cor <- cor(mod01.ord$scores, site_res01[,cor.preds]))

plot_factors(alpha, "pa", "mod01", sp_res01, site_res01)


# with 3 predictors, + offset
mod02 <- manyglm(otu.pa ~ dem500 + insideHJA + lg_YrsDisturb, 
                    family = "negative.binomial", # family = "binomial"
                    data = env.csv) # , show.time = "all" # ??
plot(mod02)

#AIC(mod02)


## anova(mod1, nBoot = 100)
# Time elapsed: 0 hr 5 min 48 sec
# Analysis of Deviance Table
# 
# Model: otu.pa ~ elevation_m + insideHJA + lg_YrsDisturb + offset(log(offset))
# 
# Multivariate test:
#   Res.Df Df.diff   Dev Pr(>Dev)   
# (Intercept)       87                          
# elevation_m       86       1 550.6     0.01 **
#   insideHJA         85       1 715.5     0.01 **
#   lg_YrsDisturb     84       1 532.9     0.01 **

# do ordination
mod02.ord <- cord(mod02)

# Site scores, species scores
site_res02 <- data.frame(mod02.ord$scores, env.csv, XY.csv)
sp_res02 <- data.frame(mod02.ord$loadings, species = colnames(otu.pa))

plot_factors(alpha, "pa", "mod02", sp_res02, site_res02)
plot_xy_factor1("pa", "mod02", site_res02)

# get correlations with selected predictors

(mod02.cor <- cor(mod02.ord$scores, site_res02[,cor.preds]))

corrplot(cor(site_res02[,cor.preds]), 
         method = "ellipse",
         type = "lower",
          mar=c(0,0,4,0))

corrplot(t(mod02.cor), is.corr = F,
         method = "ellipse", cl.pos = "n")

# add interaction of oldgrowth and elavatiaon
mod03 <- manyglm(otu.pa ~ elevation_m * oldGrowthIndex + insideHJA + lg_YrsDisturb, 
                    family = "negative.binomial", # family = "binomial"
                    data = env.csv)
plot(mod03)

# anova(mod03, nboot = 50)

dAIC <- AIC(mod01) - AIC(mod02) # will be generally positive  across species if mod2 is better... 
boxplot(dAIC); abline(h =0, col = "red")

# do ordination
mod03.ord <- cord(mod03)

# Site scores, species scores
site_res03 <- data.frame(mod03.ord$scores, env.csv, XY.csv)
sp_res03 <- data.frame(mod03.ord$loadings, species = colnames(otu.pa))
(mod03.cor <- cor(mod03.ord$scores, site_res03[,cor.preds]))

plot_factors(alpha, "pa", "mod03", sp_res03, site_res03)
plot_xy_factor1("pa", "mod03", site_res03)
plot_xy_factor2("pa", "mod03", site_res03)

mod04 <- manyglm(otu.pa ~ 
                   elevation_m * oldGrowthIndex + insideHJA + 
                   lg_YrsDisturb + mean.NDVI + lg_DistRoad + canopyHeight_m, 
                family = "negative.binomial", # family = "binomial"
                data = env.csv)

mod04.ord <- cord(mod04)
site_res04 <- data.frame(mod04.ord$scores, env.csv, XY.csv)
sp_res04 <- data.frame(mod04.ord$loadings, species = colnames(otu.pa))
(mod04.cor <- cor(mod04.ord$scores, site_res04[,cor.preds]))

t_corrplot(mod04.cor, "mod04")

plot_factors(alpha, "pa", "mod04", sp_res04, site_res04)
plot_xy_factor1("pa", "mod04", site_res04)
plot_xy_factor2("pa", "mod04", site_res04)


form05 <- formula(otu.pa ~ elevation_m * oldGrowthIndex + insideHJA + canopyHeight_m * elevation_m + lg_YrsDisturb + mean.NDVI + lg_DistRoad)
mod05 <- manyglm(form05, 
                family = "negative.binomial", # family = "binomial"
                data = env.csv)

mod05.ord <- cord(mod05)
site_res05 <- data.frame(mod05.ord$scores, env.csv, XY.csv)
sp_res05 <- data.frame(mod05.ord$loadings, species = colnames(otu.pa))
(mod05.cor <- cor(mod05.ord$scores, site_res05[,cor.preds]))


plot_xy_factor1("pa", "mod05", site_res05)
ggsave("Hmsc_CD/local/factor1_8predictors.png")
plot_xy_factor2("pa", "mod05", site_res05)
ggsave("Hmsc_CD/local/factor2_8predictors.png")


png("Hmsc_CD/local/corplot_pa.png", width = 300, height = 200, units = "mm", res = 100)
par(mfrow =c(1,6))
t_corrplot(mod0.cor, "mod0")
t_corrplot(mod01.cor, "mod01")
t_corrplot(mod02.cor, "mod02")
t_corrplot(mod03.cor, "mod03")
t_corrplot(mod04.cor, "mod04")
t_corrplot(mod05.cor, "mod05")
dev.off()

dAIC <- list(mod0_mod01 = AIC(mod0)-AIC(mod01),
             mod0_mod02 = AIC(mod0)-AIC(mod02),
             mod0_mod03 = AIC(mod0)-AIC(mod03),
             mod0_mod04 = AIC(mod0)-AIC(mod04),
             mod0_mod05 = AIC(mod0)-AIC(mod05))

par(mfrow=c(1,1))
png("Hmsc_CD/local/AIC_boxplot.png")
boxplot(dAIC); abline(h = 0, col = "red", lty = 2, las = 2)
dev.off()

png("Hmsc_CD/local/cor_boxplot.png")
boxplot(list(cor0 = abs(mod0.cor),
             cor1 = abs(mod01.cor),
             cor2 = abs(mod02.cor),
             cor3 = abs(mod03.cor),
             cor4 = abs(mod04.cor),
             cor5 = abs(mod05.cor)))
dev.off()


getwd()
save.image("Hmsc_CD/working_cd/wk_tmp_cor_ecocop_impage.rdata")
