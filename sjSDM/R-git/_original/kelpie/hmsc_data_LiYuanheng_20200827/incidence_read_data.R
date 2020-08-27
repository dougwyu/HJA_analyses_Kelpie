# MODIFY THIS SCRIPT SO THAT IT RUNS WITH YOUR OWN DATA



# You need to provide an SXY file.

# The files TP and P are optional, so indicate with TRUE/FALSE if they are included or not

is.TP = F
is.P = F
	


# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 

SXY = read.csv('incidence_lidar_mulspec_5_sample_by_species_corr_table_F2308_minimap2_20200221_kelpie20200214.csv', header=T, sep=',', stringsAsFactors = F, na.strings='NA') 
	
# Modify the next three lines to split your SXY file to components that relate to

# S: study design, including units of study and their possible coordinates (named as Route_x and Route_y to indicate that they relate to the units of Route)

# X: covariates to be used as predictors

# Y: species data

# If you don't have variables that define the study design, indicate this by S=NULL

# If you don't have covariate data, indicate this by X=NULL

S=SXY[,c(1:3,13:15)]
S$Route_x = S$UTM_E
S$Route_y = S$UTM_N
S$UTM_N = NULL
S$UTM_E = NULL

X=SXY[,c(4:12,16:53)]

Y=SXY[,54:573]



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

