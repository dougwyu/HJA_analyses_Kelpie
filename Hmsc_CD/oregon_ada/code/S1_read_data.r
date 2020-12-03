#### Data preparation #### 

## Christian D
## reworked 27/11/2020

## Only local: 
# setwd("J:/UEA/Oregon/Oproject/oregon_ada")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"  ## whereever you send the sub from is wd
## on ADA, HOME is "/gpfs/home/hsp20azu"
# setwd("~/Oregon_ada")

# Data files, separated by trap and season
files <- c("otuenv_M1S1_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M2S1_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M1S2_minimap2_20200929_kelpie20200927.csv",
           "otuenv_M2S2_minimap2_20200929_kelpie20200927.csv")

# train data
SXY.train <- read.csv(file.path("data", files[1]))

cols <- colnames(SXY.train)
sppCols <- cols[grepl("__", cols)]
Y.train <- SXY.train[, sppCols]

# convert to presence absnece
Y.train[Y.train > 0] <- 1

# Filter at > 5 species occurrences
Y.train <- Y.train[,colSums(Y.train)>5]

## Study design data
S.train <- SXY.train[, 1:6]

# Make unique id
S.train$site_trap_period <- paste(S.train$SiteName, S.train$trap, S.train$period, sep = "_")

X.all <- SXY.train[, 7:101]

## Summarise and log predictors(see S0_explore_data.r)
# Initial selection made looking at correlations and pair plots (see data S0_explore_data.r in local)

# mean of bands across all dates
X.all$B1_mean <- apply(X.all[, cols[grepl("B1_", cols)]], 1, mean)
X.all$B4_mean <- apply(X.all[, cols[grepl("B4_", cols)]], 1, mean)

X.all$lg_DistRoad <- log10(X.all$distToRoad_m)
X.all$lg_YrsDisturb <- log10(X.all$YrsSinceDist)

X.all$lg_cover2m_max <- log10(X.all$l_Cover_2m_max)
X.all$lg_cover2m_4m <- log10(X.all$l_Cover_2m_4m)
X.all$lg_cover4m_16m <- log10(X.all$l_Cover_4m_16m+0.01)

preds <- c("clearcut","insideHJA","oldGrowthIndex","elevation_m","canopyHeight_m", "precipitation_mm","distToStream_m","mean.NDVI","mean.EVI","mean.green","mean.wet","l_p25","l_rumple","B1_mean","B4_mean","lg_DistRoad","lg_YrsDisturb","lg_cover2m_max","lg_cover2m_4m","lg_cover4m_16m")

X.train <- X.all[, preds]

# Make factors
# str(X.train)
X.train$clearcut <- factor(X.train$clearcut, levels = c("no", "yes"))
X.train$insideHJA <- factor(X.train$insideHJA, levels = c("no", "yes"))

## Taxonomy from colnames ###
# R2431.10__Insecta_Coleoptera_Scraptiidae_Anaspis_olympiae_BOLD_ACC3109_size.60 
spp <- data.frame(col = colnames(Y.train))
#spp$spID <- sprintf("sp%03d", 1:nrow(spp))
spp$BOLD <- sub(".*_BOLD_([[:alnum:]]*)_.*", "\\1", spp$col)
spp$class <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\1", spp$col)
spp$order <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\2", spp$col)
spp$family <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\3", spp$col)
spp$genus <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\4", spp$col)
spp$epiphet <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\5", spp$col)

# Change "NA" to NA
for(c in seq_along(spp)[-1]) spp[,c] <- sub("NA", NA, spp[,c])

spp$sciname <- ifelse(!is.na(spp$genus) & !is.na(spp$epiphet), paste(spp$genus, spp$epiphet), NA) 
# head(spp)

# check for duplicates
sum(is.na(spp$sciname)); sum(!is.na(spp$sciname))
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

rownames(spp) <- NULL

P <- ape::as.phylo(~class/order/family/genus/col, data = spp, collapse = F)

P$edge.length = rep(1, length(P$edge)) # make all lengths eqaul between tax levels
ape::is.rooted(P)

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

# Read in the species names as rownames, not as a column of the matrix
# The script below checks if the species names in TP are identical and in the same order as in Y
all(P$tip.label %in% colnames(Y.train))

# clean up 
rm(spp, X.all, c, cols, dups, files, i, ind, preds, sppCols, tax.cols)

# Data goes straight to model set up. 
