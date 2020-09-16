# packages
# library(sjSDM)
library(tidyverse)
library(fs)
library(glue)
library(RColorBrewer)
library(here)

rundate <- 20200915 # sjsdm_cv run date
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
taxon <- "Dolichovespula"

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

