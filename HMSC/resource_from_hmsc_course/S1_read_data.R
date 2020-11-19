setwd("h572/hmsc_course/li_yuanheng/data/")
# OR ON A MAC:
# setwd("/Volumes/group/h572/hmsc_course/li_yuanheng/data/")

# MODIFY THIS SCRIPT SO THAT IT RUNS WITH YOUR OWN DATA

# You need to provide an SXY file.
# The files TP and P are optional, so indicate with TRUE/FALSE if they are included or not
is.TP = TRUE
is.P = FALSE


# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 

SXY = read.csv('SXY.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA') 
	
# Modify the next three lines to split your SXY file to components that relate to

# S: study design, including units of study and their possible coordinates (named as Route_x and Route_y to indicate that they relate to the units of Route)

# X: covariates to be used as predictors

# Y: species data

# If you don't have variables that define the study design, indicate this by S=NULL

# If you don't have covariate data, indicate this by X=NULL

S=SXY[,c(1:3,13:15)]

X=SXY[,c(4:12,16:55)]

Y=SXY[,56:573]



# What is not always easy is to decide what goes to S and what to X.

# As a general rule, include in S those variables that you think should be modelled as random effect,

# and in X those that you think should be modelled as fixed effects.

# Don't worry if you are not sure which one is the "right choice", we will discuss this with you.



# Check that the data looks as it should!

View(S)

View(X)

View(Y)

# check that community data are numeric and have finite numbers. If the script

# writes "Y looks OK", you are ok.

if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {

    print("Y looks OK")

} else {

	print("Y should be numeric and have finite values")	}

# Check that the study design data do not have missing values (they are allowed for Y but not S, X, P or Tr)

if (any(is.na(S))) {

  print("S has NA values - not allowed for")

} else {

  print("S looks ok")	}

# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)

if (any(is.na(X))) {

  print("X has NA values - not allowed for")

} else {

  print("X looks ok")	}




# READING IN TP: traits (T) and/or phylogenetic information in table format (P)
if(is.TP){
  # Read in the species names as rownames, not as a column of the matrix
  TP = read.csv("TP.csv", stringsAsFactors=TRUE,row.names = 1)
  # The script below checks if the species names in TP are identical and in the same order as in Y
  # If the script prints "species names in TP and SXY match", you are ok.
  # If it says that they do not match, you need to modify the files so that they match 
  if(all(rownames(TP)==colnames(Y))) {
    print("species names in TP and SXY match")
  } else{
    print("species names in TP and SXY do not match")
  }
  # Modify the next two lines to split your TP file to components that relate to
  # Tr: species traits (note that T is a reserved word in R and that's why we use Tr)
  # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species
  # If you don't have trait data, indicate this by Tr=NULL. 
  # If TP does not have phylogenetic data (because you don't have such data at all, or because
  # it is given in tree-format, like is the case in this example), indicate this with P=NULL 
  Tr = NULL
  P = NULL # taxonomic levels
  # Check that the data looks as it should!
  View(Tr)
  View(P)
  # Check that the Tr data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(Tr))) {
    print("Tr has NA values - not allowed for")
  } else {
    print("Tr looks ok")	}
  # Check that the phylogenetic/taxonomic data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(P))) {
    print("P has NA values - not allowed for")
  } else {
    print("P looks ok")	}
}

# Creating a proxy tree from taxonomic data on subfamily, tribe, genus and species 
# as found in the P file.

# Creating a proxy tree from taxonomic data on order, family, genus and species in the TP file.
# Adding a dummy genus and a dummy species to each family/genus, respectively, where the 
# taxonomy is not resolved to the genus and/or species level.
Tax = data.frame(TP)
Tax = apply(Tax, 2, as.character)
Tax[which(is.na(Tax[,2])),2] = paste(Tax[which(is.na(Tax[,2])),1], "dF", sep="_")
Tax[which(is.na(Tax[,3])),3] = paste(Tax[which(is.na(Tax[,3])),2], "dG", sep="_")
Tax[which(is.na(Tax[,4])),4] = paste(Tax[which(is.na(Tax[,4])),3], "dS", sep="_")
Tax = data.frame(Tax)
Tax$Yspecies<-as.factor(names(Y))

library(ape)
P = as.phylo(~order/family/genus/Yspecies, data = Tax, collapse = FALSE)
P$edge.length = rep(1, length(P$edge))
plot(P, cex=0.3)

Y = as.matrix(Y)
save(S,X,Y,P,Tax, file="allData.R")

