## For version sjSDM 0.1.8
options(echo=TRUE) # if you want see commands in output file
Sys.setenv(RETICULATE_PYTHON="/gpfs/scratch/hsp20azu/sjSDM_env/bin/python")
myPaths <- .libPaths()
myPaths <- c("/gpfs/scratch/hsp20azu/newrlib",myPaths[2],myPaths[1])
.libPaths(myPaths)
library("sjSDM")
packageVersion("sjSDM")
# [1] ‘0.1.8’
