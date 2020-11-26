#### Data preparation #### 

## Christian D
## 12/11/2020

# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

## Only local: 
# setwd("J:/UEA/Oregon/Oregon_winscp")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"  - 
# if you send the .sub from here. whereever you send the sub from is wd
# with folders Oregon, etc... 

# setwd("~/Oregon") # tilde expansion relative to HOME system env variable
## on ADA, HOME is "/gpfs/home/hsp20azu"
setwd("~/Oregon_winscp")  # ./ relative to working directory 
dir()

# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 
# data is in data directory
	
SXY = read.csv(file.path('data', 'incidence_lidar_mulspec_5_sample_by_species_table_F2308_minimap2_20200929_kelpie20200927.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')

# head(SXY)

# quasi -- indicator of species abundance
# SXY = read.csv(file.path("data", 'quasiP_lidar_mulspec_5_sample_by_species_table_F2308_minimap2_20200929_kelpie20200927.csv'), header=T, sep=',', stringsAsFactors = F, na.strings='NA')


## Initial data filtering ###

# train data only S1 M1
# test data S1 M2
table(SXY$trap, useNA = "always")
table(SXY$session, useNA = "always")
addmargins(table(SXY[,c("trap", "session")]))

SXY.train <- subset(SXY, trap == "M1" & session == "S1")
SXY.test <- subset(SXY, trap == "M2" & session == "S1")

## Filter for minumum incidences
head(SXY[, 56:60])

# Species richness per site
# hist(rowSums(SXY[,56:580]))

## prevalence per species
# hist(colMeans(SXY[,56:580]))

# summary of species incidences
summary(colSums(SXY[,56:580]))
# already minimum of 5 incidences

summary(colMeans(SXY[,56:580]))

## Only include unscaled variables and let Hmsc scale internally on each run
Xcols <- colnames(SXY[,c(5:13,16:55)])
sort(Xcols)

X_cols_use <- Xcols[!grepl("scale", Xcols)]
X_cols_scale <- Xcols[grepl("scale", Xcols)]
X_cols_use
X_cols_scale

rm(Xcols, X_cols_scale)

S <- SXY[,c(1:4,14:15)]

S.train <- SXY.train[,c(1:4,14:15)]
X.train <- SXY.train[,X_cols_use]
Y.train <- SXY.train[,56:580]

S.test <- SXY.test[,c(1:4,14:15)]
X.test <- SXY.test[,X_cols_use]
Y.test <- SXY.test[,56:580]


## Taxonomy form colnames ###
# R2431.10__Insecta_Coleoptera_Scraptiidae_Anaspis_olympiae_BOLD_ACC3109_size.60 
spp <- data.frame(col = colnames(SXY[,56:580]))
#spp$spID <- sprintf("sp%03d", 1:nrow(spp))
spp$BOLD <- sub(".*_BOLD_([[:alnum:]]*)_.*", "\\1", spp$col)
spp$class <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\1", spp$col)
spp$class <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\1", spp$col)
spp$order <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\2", spp$col)
spp$family <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\3", spp$col)
spp$genus <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\4", spp$col)
spp$epiphet <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\5", spp$col)

# Change "NA" to NA
for(c in seq_along(spp)[-1]) spp[,c] <- sub("NA", NA, spp[,c])

spp$sciname <- ifelse(!is.na(spp$genus) & !is.na(spp$epiphet), paste(spp$genus, spp$epiphet), NA) 
head(spp)


# check for duplicates
sum(is.na(spp$sciname))
sum(!is.na(spp$sciname))
length(levels(spp$sciname))
sum(duplicated(spp$sciname[!is.na(spp$sciname)]))
ind <- which(duplicated(spp$sciname[!is.na(spp$sciname)]))
dups <- spp$sciname[!is.na(spp$sciname)][ind]
spp[spp$sciname %in% dups,]

# two BOLD are assigned to same species...  fix later?
dupCols <- spp[spp$sciname %in% dups,"col"]
dupCols

## duplicated species in these columns
# "R1240.70__Insecta_Diptera_Syrphidae_Blera_scitula_BOLD_ABY7981_size.998"
# "R6954.9__Insecta_Diptera_Syrphidae_Blera_scitula_BOLD_AAI8832_size.2230"
# "R3500.3__Insecta_Hymenoptera_Apidae_Bombus_sitkensis_BOLD_ACU8557_size.1191"
# "R1891.25__Insecta_Hymenoptera_Apidae_Bombus_sitkensis_BOLD_AAI4757_size.45"

sum(is.na(spp$order))

## add dummy taxa
sum(is.na(spp$family))
#spp$family[is.na(spp$family)] <- paste0(spp$order[is.na(spp$family)], "_dF")
spp$family[is.na(spp$family)] <- sprintf("fam%03d", 1:sum((is.na(spp$family))))

# Add dummy genus
sum(is.na(spp$genus))
#spp$genus[is.na(spp$genus)] <- paste0(spp$family[is.na(spp$genus)], "_dG")
spp$genus[is.na(spp$genus)] <- sprintf("gen%03d", 1:sum((is.na(spp$genus))))

head(spp, 20)

# convert to factors for ape
spp <- spp[order(spp$class, spp$order, spp$family, spp$genus),]
tax.cols <- c("class", "order", "family", "genus", "epiphet", "BOLD", "sciname", "col")
for(i in tax.cols) spp[,i] <- factor(spp[,i])
str(spp)

rownames(spp) <- NULL

P <- ape::as.phylo(~class/order/family/genus/col, data = spp, collapse = F)
str(P)

# plot(P, show.tip.label = F, edge.width = 0.1)
# plot(P, type = "fan", show.tip.label = F, tip.color = as.numeric(spp$order))
# plot(P, type = "radial", show.tip.label = F)

P$edge.length = rep(1, length(P$edge)) # make all lengths eqaul between tax levels
# plot(P)
ape::is.rooted(P)

dist <- ape::cophenetic.phylo(P)
hist(dist) # distances are just number of taxonomic levels apart.... ie it doesn't matter whether families, genera, etc are in correct order.

# ape::write.tree(spp.tree, "P.tre")


## Filter for prevalence < 0.05???
## Why is other data set less species?
## ver low number of sites for training and testing, compares to species

## Traits?

## dummy genus and family??

## predictor choices - see exel doc for description
# cbind(X_cols_use)
# [1,] "elevation"         
# [2,] "canopy.ht"         
# [3,] "min.T"             
# [4,] "max.T"             
# [5,] "precipitation"     
# [6,] "metre.road"        
# [7,] "metre.stream"      
# [8,] "yrs.disturb.min"   
# [9,] "hja"    # hja inside (no logging really, nearer to primary forest) or outside experimental forest (logging)
# make a domindant land  cover variable around each point. % of structure... 

# [10,] "lysis.ratio"   #  ignore    lysis buffer ratio. different samples are different sizes - control for this?
# [11,] "spike.mean"        # 
# [12,] "l_Cover_2m_4m"     
# [13,] "l_Cover_2m_4m_all" 
# [14,] "l_Cover_2m_max"    
# [15,] "l_Cover_2m_max_all"
# [16,] "l_Cover_4m_16m"    # understorey
# [17,] "l_p25"             
# [18,] "l_p25_all"         
# [19,] "l_p95"             
# [20,] "l_p95_all"         
# [21,] "l_rumple"          
# [22,] "mean.NDVI"         
# [23,] "mean.EVI"          
# [24,] "mean.bright"       
# [25,] "mean.green"        
# [26,] "mean.wet" 



# What is not always easy is to decide what goes to S and what to X.
# As a general rule, include in S those variables that you think should be modelled as random effect,
# and in X those that you think should be modelled as fixed effects.
# Don't worry if you are not sure which one is the "right choice", we will discuss this with you.


# check that community data are numeric and have finite numbers. If the script
# writes "Y looks OK", you are ok.
is.numeric(as.matrix(Y.train)) || is.logical(as.matrix(Y.train)) && is.finite(sum(Y.train, na.rm=TRUE))

# Check that the study design data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(S.train))) {
  print("S has NA values - not allowed for")
} else {
  print("S looks ok")	}

# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(X.train))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}


# READING IN TP: traits (T) and/or phylogenetic information in table format (P)

  # Read in the species names as rownames, not as a column of the matrix
# The script below checks if the species names in TP are identical and in the same order as in Y
all(P$tip.label %in% colnames(Y.train))

save(SXY, S,S.train, X.train, Y.train, S.test, Y.test, X.test, P, file = file.path("data", "allData.Rdata"))
