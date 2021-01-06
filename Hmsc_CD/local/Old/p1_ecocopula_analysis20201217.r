#### ecocopula analysis

library(dplyr)
library(mvabund)
library(ecoCopula)
library(here)
library(glue)
library(vegan)
library(ggplot2)
library(RColorBrewer)

here()

samtoolsfilter <- "F2308" # F2308 filter only
samtoolsqual <- "q48"
minimaprundate <- 20200929
kelpierundate <- 20200927
primer <- "BF3BR2"

outputidxstatstabulatefolder <- glue::glue("outputs_minimap2_{minimaprundate}_{samtoolsfilter}_{samtoolsqual}_kelpie{kelpierundate}_{primer}_vsearch97")

otuenv <- read.csv(here("Kelpie_maps", 
                        outputidxstatstabulatefolder, glue("sample_by_species_table_{samtoolsfilter}_minimap2_{minimaprundate}_kelpie{kelpierundate}_uncorr.csv")))


# M1S1
trap <- "M1"
period <- "S1"
otuenv <- otuenv %>% 
  filter(trap == trap[[1]] & period == period[[1]]) 

XY.csv <- otuenv %>% 
  select(UTM_E, UTM_N) %>% 
  #scale() %>%  # keep orignal for plotting
  as_tibble()

# keep OTUs with >=5 incidences
# original read number abundance
minocc <- 5 # set to high number (e.g. 20) for testing
otu.ab.csv <- otuenv %>% select(contains("__"))
otu.ab.csv <- otu.ab.csv[ , specnumber(otu.ab.csv, MARGIN = 2) >= minocc] 

# log(FSL) correction and scale to quasiprobability
otu.qp.csv <- otu.ab.csv %>% 
  mutate(across(contains("__"), 
                ~ .x /(otuenv$COISpike_sum*otuenv$lysis_ratio))) %>% 
  mutate(across(contains("__"), ~ log(.x + 0.001))) %>% 
  mutate(across(contains("__"), ~ scales::rescale(.x))) # {scales}
max(otu.qp.csv) == 1 # should be TRUE

otu.qp.csv[1:10, 1:5]


# convert to presence/absence data
otu.pa.csv <- otu.ab.csv
otu.pa.csv[otu.pa.csv > 0] <- 1
min(colSums(otu.pa.csv)) == minocc # should be TRUE

# env covariates
otuenv %>% 
  select(!contains("__"), -UTM_E, -UTM_N, -starts_with("nor")) %>% 
  names(.)
#  [1] "SiteName"           "trap"               "period"            
#  [4] "lysis_ratio"        "COISpike_sum"       "clearcut"          
#  [7] "insideHJA"          "oldGrowthIndex"     "elevation_m"       
# [10] "canopyHeight_m"     "minT_annual"        "maxT_annual"       
# [13] "precipitation_mm"   "distToRoad_m"       "distToStream_m"    
# [16] "YrsSinceDist"       "B1_20180717"        "B2_20180717"       
# [19] "B3_20180717"        "B4_20180717"        "B5_20180717"       
# [22] "B6_20180717"        "B7_20180717"        "B10_20180717"      
# [25] "B11_20180717"       "NDVI_20180717"      "EVI_20180717"      
# [28] "B_20180717"         "G_20180717"         "W_20180717"        
# [31] "B1_20180726"        "B2_20180726"        "B3_20180726"       
# [34] "B4_20180726"        "B5_20180726"        "B6_20180726"       
# [37] "B7_20180726"        "B10_20180726"       "B11_20180726"      
# [40] "NDVI_20180726"      "EVI_20180726"       "B_20180726"        
# [43] "G_20180726"         "W_20180726"         "B1_20180802"       
# [46] "B2_20180802"        "B3_20180802"        "B4_20180802"       
# [49] "B5_20180802"        "B6_20180802"        "B7_20180802"       
# [52] "B10_20180802"       "B11_20180802"       "NDVI_20180802"     
# [55] "EVI_20180802"       "B_20180802"         "G_20180802"        
# [58] "W_20180802"         "B1_20180818"        "B2_20180818"       
# [61] "B3_20180818"        "B4_20180818"        "B5_20180818"       
# [64] "B6_20180818"        "B7_20180818"        "B10_20180818"      
# [67] "B11_20180818"       "NDVI_20180818"      "EVI_20180818"      
# [70] "B_20180818"         "G_20180818"         "W_20180818"        
# [73] "mean.NDVI"          "mean.EVI"           "mean.bright"       
# [76] "mean.green"         "mean.wet"           "mean.NDVI.scale"   
# [79] "mean.EVI.scale"     "mean.green.scale"   "mean.bright.scale" 
# [82] "mean.wet.scale"     "l_Cover_2m_max"     "l_Cover_2m_max_all"
# [85] "l_Cover_2m_4m"      "l_Cover_2m_4m_all"  "l_Cover_4m_16m"    
# [88] "l_p25"              "l_p25_all"          "l_p95"             
# [91] "l_p95_all"          "l_rumple"          

# how many otu in total?
sum(grepl("__", colnames(otuenv)))

# select subset of covariates chosen by Christian
# remove OTUs, XY, and normalised NDVI and EVI
# average, optionally log, select, and scale env covariates
env.csv <- otuenv %>% 
  select(!contains("__"), -UTM_E, -UTM_N, -starts_with("nor")) %>%
  mutate(B1_mean = rowMeans(across(starts_with("B1_")))) %>% 
  mutate(B4_mean = rowMeans(across(starts_with("B4_")))) %>%  
  mutate(lg_DistStream = log(distToStream_m + 0.001)) %>% 
  mutate(lg_DistRoad = log(distToRoad_m + 0.001)) %>% 
  mutate(lg_YrsDisturb = log(YrsSinceDist + 0.001)) %>% 
  mutate(lg_cover2m_max = log(l_Cover_2m_max + 0.001)) %>% 
  mutate(lg_cover2m_4m = log(l_Cover_2m_4m + 0.001)) %>%
  mutate(lg_cover4m_16m = log(l_Cover_4m_16m + 0.001)) %>%   
  select(SiteName, trap, period, lysis_ratio, COISpike_sum, clearcut, insideHJA, oldGrowthIndex, elevation_m, canopyHeight_m, precipitation_mm, mean.NDVI, mean.EVI, mean.green, mean.wet, l_p25, l_rumple, B1_mean, B4_mean, lg_DistStream, lg_DistRoad, lg_YrsDisturb, lg_cover2m_max, lg_cover2m_4m, lg_cover4m_16m) %>% 
  mutate(across(where(is.numeric) & !lysis_ratio & !COISpike_sum, scale)) %>%  # scale here
  mutate(offset = (1/COISpike_sum * lysis_ratio)) %>% 
  relocate(offset, .after = COISpike_sum)

str(env.csv)

## offset term
# the smaller the lysis_ratio, the larger the original sample. 
# the larger the COISpike_sum, the larger the original sample
# to create an offset term, I thus use 1/COISpike_sum*lysis_ratio

source("Hmsc_CD/local/ecoCopula_plot_fn.r")

## do the models
# mvabund models, mod0 and mod1
abund <- "pa"
glue("otu.{abund}.csv")

# make a pa matrix as mvabund object
otu.pa <- mvabund(otu.pa.csv)
otu.qp <- mvabund(otu.qp.csv)

is.mvabund(otu.pa)
is.mvabund(otu.qp)
str(otu.pa)

# with no env covariates
mod0 <- manyglm(otu.pa ~ 1 + offset(log(offset)), 
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
                 color = elevation_m, 
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

mod05 <- manyglm(otu.pa ~ elevation_m + offset(log(offset)), 
                family = "negative.binomial", # family = "binomial"
                data = env.csv) # , show.time = "all" # ??

mod05.ord <- cord(mod05)
site_res05 <- data.frame(mod05.ord$scores, env.csv, XY.csv)
sp_res05 <- data.frame(mod05.ord$loadings, species = colnames(otu.pa))
(mod05.cor <- cor(mod05.ord$scores, site_res05[,cor.preds]))


# with 3 predictors, + offset
mod1 <- manyglm(otu.pa ~ elevation_m + insideHJA + lg_YrsDisturb + offset(log(offset)), 
                    family = "negative.binomial", # family = "binomial"
                    data = env.csv) # , show.time = "all" # ??
plot(mod1)

AIC(mod1)


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
mod1.ord <- cord(mod1)

# Site scores, species scores
site_res1 <- data.frame(mod1.ord$scores, env.csv, XY.csv)
sp_res1 <- data.frame(mod1.ord$loadings, species = colnames(otu.pa))

ggplot() + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = Factor1 * alpha, 
                   yend = Factor2 * alpha), 
               data = sp_res1, 
               size = .1) +
  geom_point(aes(x = Factor1, y = Factor2, 
                 color = elevation_m, 
                 size = exp(lg_YrsDisturb), 
                 shape = as.factor(insideHJA)),
             data = site_res1) + 
  scale_color_gradientn(colours = brewer.pal(n = 10, name = "RdYlBu")) +
  theme_classic() + # or "PuOr" or "RdYlBu"
  xlab("Factor 1") +
  ylab("Factor 2")

plot_factors(4, abund, "mod1", sp_res1, site_res1)

# get correlations with selected predictors
# num.pred <- sapply(site_res1, is.numeric)
# cat(paste(colnames(site_res1[,num.pred]), collapse = '", "'))
preds <- c("Factor1", "Factor2", "offset", "oldGrowthIndex", "elevation_m", "canopyHeight_m", "precipitation_mm", "mean.NDVI", "mean.EVI", "mean.green", "mean.wet", "l_p25", "l_rumple", "B1_mean", "B4_mean", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m")

cor.preds <- c("oldGrowthIndex", "elevation_m", "canopyHeight_m", "precipitation_mm", "mean.NDVI", "mean.EVI", "mean.green", "mean.wet", "l_p25", "l_rumple", "B1_mean", "B4_mean", "lg_DistStream", "lg_DistRoad", "lg_YrsDisturb", "lg_cover2m_max", "lg_cover2m_4m", "lg_cover4m_16m")


cor(site_res1[,preds])
(mod1.cor <- cor(mod1.ord$scores, site_res1[,cor.preds]))

corrplot(cor(site_res1[,preds]), 
         method = "ellipse",
         type = "lower",
          mar=c(0,0,4,0))

corrplot(t(mod1.cor), is.corr = F,
         method = "ellipse", cl.pos = "n")

plot_corrplot(abund, "mod1", site_res1)


# add interaction of oldgrowth and elavatiaon
mod2 <- manyglm(otu.pa ~ elevation_m * oldGrowthIndex + insideHJA + lg_YrsDisturb + offset(log(offset)), 
                    family = "negative.binomial", # family = "binomial"
                    data = env.csv)
plot(mod2)

anova(mod2, nboot = 50)

dAIC <- AIC(mod1) - AIC(mod2) # will be generally positive  across species if mod2 is better... 
boxplot(dAIC); abline(h =0, col = "red")

# ordination 
mod2.ord <- cord(mod2)
site_res2 <- data.frame(mod2.ord$scores, env.csv, XY.csv)
sp_res2 <- data.frame(mod2.ord$loadings, species = colnames(otu.pa))

# make plots
plot_factors(4, abund, "mod2", sp_res2, site_res2)
plot_xy_factor1(abund = abund, model = "mod2", siteData = site_res2)

# correlations
(mod2.cor <- cor(mod2.ord$scores, site_res2[,cor.preds]))

corrplot(t(mod2.cor), is.corr = F,
         method = "ellipse", cl.pos = "n")

mod2.cor

mod3 <- manyglm(otu.pa ~ elevation_m * oldGrowthIndex + insideHJA + lg_YrsDisturb + mean.NDVI + lg_DistRoad + 
                  canopyHeight_m + offset(log(offset)), 
                family = "negative.binomial", # family = "binomial"
                data = env.csv)

mod3.ord <- cord(mod3)
site_res3 <- data.frame(mod3.ord$scores, env.csv, XY.csv)
sp_res3 <- data.frame(mod3.ord$loadings, species = colnames(otu.pa))
(mod3.cor <- cor(mod3.ord$scores, site_res3[,cor.preds]))

form4 <- formula(otu.pa ~ elevation_m * oldGrowthIndex + insideHJA + canopyHeight_m * elevation_m + lg_YrsDisturb + mean.NDVI + lg_DistRoad + offset(log(offset)))
mod4 <- manyglm(form4, 
                family = "negative.binomial", # family = "binomial"
                data = env.csv)

mod4.ord <- cord(mod4)
site_res4 <- data.frame(mod4.ord$scores, env.csv, XY.csv)
sp_res4 <- data.frame(mod4.ord$loadings, species = colnames(otu.pa))
(mod4.cor <- cor(mod4.ord$scores, site_res4[,cor.preds]))

png("Hmsc_CD/local/corplot.png", width = 300, height = 200, units = "mm", res = 100)
par(mfrow =c(1,5))
t_corrplot(mod0.cor, "mod0")
t_corrplot(mod1.cor, "mod1")
t_corrplot(mod2.cor, "mod2")
t_corrplot(mod3.cor, "mod3")
t_corrplot(mod4.cor, "mod4")
dev.off()

dAIC <- list(mod0_mod1 = AIC(mod0)-AIC(mod1),
             mod0_mod05 = AIC(mod0)-AIC(mod05),
             mod1_mod2 = AIC(mod1)-AIC(mod2),
             mod1_mod3 = AIC(mod1)-AIC(mod3),
             mod1_mod4 = AIC(mod1)-AIC(mod4))


png("Hmsc_CD/local/AIC_boxplot.png")
boxplot(dAIC); abline(h = 0, col = "red", lty = 2, las = 2)
dev.off()

png("Hmsc_CD/local/cor_boxplot.png")
boxplot(list(cor0 = abs(mod0.cor),
             cor05 = abs(mod05.cor),
             cor1 = abs(mod1.cor),
             cor2 = abs(mod2.cor),
             cor3 = abs(mod3.cor),
             cor4 = abs(mod4.cor)))
dev.off()
