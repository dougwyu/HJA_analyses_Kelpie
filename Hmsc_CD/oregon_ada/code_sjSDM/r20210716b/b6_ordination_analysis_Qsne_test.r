#### ecocopula analysis

## copied initially from oregon_ada/code/P1_ecocopula_analysis_v3.r

options(echo=TRUE) # if you want see commands in output file

# Qsne 2020-01-24 is now installed on ADA as an environment module
# It can be accessed by executing the following command
# # module add qsne/2020-01-24
# 
# For accessing in R 3.6.2  you have to run the following command inside the R
# 
# module add R/3.6.2
# R

source('code_sjSDM/qsne.r')
qsne

system("/gpfs/software/ada/qsne/2020-01-24/bin/qsne")




# 
# 
# 
# ## Get the Fisher iris data from R
# 
# get data, strip labels
X <- iris
X <- as.matrix( X[, !(colnames(X) == 'Species')] )

# put back the errors
X[35, ] <- c( 4.9, 3.1, 1.5, 0.1 )
X[38, ] <- c( 4.9, 3.1, 1.5, 0.1 )

write.table(X, file = "code_sjSDM/iris.txt", sep = "\t", row.names = F, col.names = F)

# from the command line:
## but need to initialise the output files for full output... 
system("/gpfs/software/ada/qsne/2020-01-24/bin/qsne -d2 -p30 -o code_sjSDM/qsne_out.txt -O code_sjSDM/qsne_obj.txt -S code_sjSDM/qsne_ent.txt code_sjSDM/iris.txt")


# pX <- tempfile("X")
# write.table(X, pX, sep="\t", row.names=F, col.names=F)
# pX


# 
# # testing
# 
# ## qsne function
# # qsne <- function(X, dims=2, perplexity=0, hess_rank=10, kld_ival=10, pca=F, max_iter=1000, tol=1e-6,verbose=F) ## taken this out: Y_init=NA,
# 

perplex <- "20:30"

rr <- qsne(X, dims=2, perplexity = perplex, hess_rank=10, pca=F, max_iter=1000, kld_ival=2, tol=1e-5)
print(rr$res_Y)
getwd()
save(rr, file = "code_sjSDM/qsne_test.rdata")