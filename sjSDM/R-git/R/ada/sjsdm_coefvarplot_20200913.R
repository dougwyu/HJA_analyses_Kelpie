# packages
library(sjSDM)
library(here)
library(tidyverse)
library(fs)
library(glue)
library(RColorBrewer)

rundate <- 20200907 # sjsdm_cv run date
envvar <- "gismslidar" # gismslidar, mslidar, gis, ms, lidar
abund <- "qp" # "qp" # pa is 0/1 data, qp is quasiprob data

minocc <- 5 # minimum occupancy (incidence) per OTU, value from dataprep.Rmd

resultsfolder <- glue("results_{rundate}_{minocc}minocc_{envvar}_{abund}_loocv")
resultsfolder
datafolder <- glue("data_{rundate}_{minocc}minocc_{envvar}")
datafolder

# load results
result <- readRDS(here("results", "crossvalidation", resultsfolder, glue("sjsdm_result_HJA_{rundate}.RDS")))
summary.p <- readRDS(here("results", "crossvalidation", resultsfolder, glue("sjsdm_summary.p_HJA_{rundate}.RDS")))

# source function
source(here("R", "ada", "coef.figure_20200912.R"))

# set variables
minsize <- 1
maxsize <- 200000L # default 200000 to include all large OTUs
taxon <- "Coleoptera"

# function(summary.p, result, minsize, maxsize=200000, taxon="all")
p1 <- coef.figure(summary.p, result, minsize, maxsize, taxon)

p1 

pdf(file = here("results", "crossvalidation", resultsfolder, 
      glue("coef_figure_min{minsize}_max{maxsize}_{taxon}_{rundate}.pdf")),
    width = 20, 
    height = 11
    )
print(p1)
dev.off()




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
