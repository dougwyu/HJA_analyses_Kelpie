#### Data preparation #### 

## Christian D
## 12/11/2020

# source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

## Only local: 
setwd("J:/UEA/Oregon/Oproject/oregon_ada/results/RRR_test")
dir()

fp <- "otuenv_M1S1_minimap2_20200929_kelpie20200927.csv"

# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 
# data is in data directory

# train data
SXY.train <- read.csv(file.path("data", fp))

cols <- colnames(SXY.train)
sppCols <- cols[grepl("__", cols)]

# which are the species columns?
range(which(cols %in% sppCols))

Y.train <- SXY.train[, sppCols]
#summary(Y.train[,1:10])

# convert to presence absnece
Y.train[Y.train > 0] <- 1
#table(unlist(Y.train), useNA = "always")

# filter species
# sp richness (per site)
#hist(rowSums(Y.train))
#summary(rowSums(Y.train))

# no of records per species
#hist(colSums(Y.train))
#summary(colSums(Y.train))

# prevalence (proportion of sites occupied - per species)
#hist(colSums(Y.train)/nrow(Y.train))
#summary(colSums(Y.train)/nrow(Y.train))

#summary(colSums(Y.train[colSums(Y.train)>5])/nrow(Y.train))

## Make a small subset for testing RRR, with mainly common, but some rare too

# first subset some sites - preference for more species
rich <- rowSums(Y.train); hist(rich)

set.seed(999)
rKeep <- sample(1:nrow(Y.train), size = 50, prob = rich^2)
Y.train <- Y.train[rKeep,]

# Filter at > 5 species occurrences
Y.train <- Y.train[,colSums(Y.train)>5]

## subset species, keeping some rare, mainly common
prev <- colSums(Y.train)/nrow(Y.train)
hist(prev)

Y.train <- Y.train[,sample(seq_along(Y.train), size = 15, prob = prev^2)] # samples biased to higher prevalence
prev <- colSums(Y.train)/nrow(Y.train); hist(prev)

## Study design data
S.train <- SXY.train[rKeep, 1:6]
# Make unique id
S.train$site_trap_period <- paste(S.train$SiteName, S.train$trap, S.train$period, sep = "_")
#head(S.train)


X.all <- SXY.train[rKeep, 7:101]
#colnames(X.all)

## Summarise and log predictors(see S0_explore_data.r)
# Initial selection made looking at correlations and pair plots

# mean of bands across all dates
X.all$B1_mean <- apply(X.all[, cols[grepl("B1_", cols)]], 1, mean)
X.all$B4_mean <- apply(X.all[, cols[grepl("B4_", cols)]], 1, mean)

X.all$lg_DistRoad <- log10(X.all$distToRoad_m)
X.all$lg_YrsDisturb <- log10(X.all$YrsSinceDist)

X.all$lg_cover2m_max <- log10(X.all$l_Cover_2m_max)
X.all$lg_cover2m_4m <- log10(X.all$l_Cover_2m_4m)
X.all$lg_cover4m_16m <- log10(X.all$l_Cover_4m_16m+0.01)

X.all$YrsDist_gt00 <- X.all$YrsSinceDist > 100

# subset some predictors
preds <- c("clearcut","insideHJA","oldGrowthIndex","elevation_m","canopyHeight_m", "mean.NDVI","mean.EVI","mean.green","mean.wet","l_p25","l_rumple","B1_mean","B4_mean","lg_DistRoad","lg_YrsDisturb","lg_cover2m_max","lg_cover2m_4m")

X.train <- X.all[, preds]

# Make factors
str(X.train)
X.train$clearcut <- factor(X.train$clearcut, levels = c("no", "yes"))
X.train$insideHJA <- factor(X.train$insideHJA, levels = c("no", "yes"))

## Taxonomy from colnames ###
# R2431.10__Insecta_Coleoptera_Scraptiidae_Anaspis_olympiae_BOLD_ACC3109_size.60 
spp <- data.frame(col = colnames(Y.train))
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
# BOLD is unique species identifiers

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
# plot(P, type = "radial", show.tip.label = F)

P$edge.length = rep(1, length(P$edge)) # make all lengths eqaul between tax levels
ape::is.rooted(P)

dist <- ape::cophenetic.phylo(P)
# hist(dist) # distances are just number of taxonomic levels apart.... 
# ie it doesn't matter whether families, genera, etc are in correct order.

# ape::write.tree(spp.tree, "P.tre")

## very low number of sites for training and testing, compares to species
## Traits?

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

# clean up 
rm(dist, spp, X.all, c, cols, dups, i, ind, preds, sppCols, tax.cols, fp, prev, rich, rKeep)

# save
save(SXY.train, S.train, X.train, Y.train, P, file = file.path("data", "allData.Rdata"))
# S.test, Y.test, X.test,