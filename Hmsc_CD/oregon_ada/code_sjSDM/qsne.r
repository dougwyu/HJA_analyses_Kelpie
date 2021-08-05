#!/usr/bin/Rscript

qsne <- function(X, dims=2, perplexity=0, hess_rank=10, kld_ival=10, max_iter=1000, tol=1e-6, nCores = NULL) {
  # Y_init=NA,
  
  if(is.null(nCores)) nCores <- parallel::detectCores()-1
  print(nCores)
  
  # temp paths
  pY <- tempfile("Y") # results
  pO <- tempfile("obj")
  pE <- tempfile("ent")
  pX <- tempfile("X")
  pYi <- tempfile("Yi")
  write.table(X, pX, sep="\t", row.names=F, col.names=F)
  #write.table(Y_init, pYi, sep="\t", row.names=F, col.names=F)
  
  #cmd <- sprintf("/gpfs/software/ada/qsne/2020-01-24/bin/qsne -d %d -p %s -m %d -K %d -o %s -O %s -S %s %s", dims, perplexity, hess_rank, kld_ival, pY, pO, pE, pX) # pYi
  cmd <- sprintf("/gpfs/software/ada/qsne/2020-01-24/bin/qsne -d %d -p %s -m %d -K %d -k %d -T %d -o %s %s", dims, perplexity, hess_rank, kld_ival, max_iter, nCores, pY, pX) # pYi pE  pO
  # run
  system(cmd)
  #return(list(res_Y=read.delim(pY, header=F), res_o=read.delim(pO, header=F), res_E=read.delim(pE, header=F)))
  return(list(res_Y=read.delim(pY, header=F)))
}

# # testing
# shell.exec(pX)
# shell.exec(pYi)
# in_X=read.delim(pX, header=F)
# perplexity <- 20
# perplexity <- "20:40"
# sprintf("qsne -p %s", perplexity) # needs to be %s to accept either numeric or character

# X <- read.delim('./test/test1.in', header=F)
# Yi <- read.delim('./test/test1.init', header=F)
# rr <- qsne(X, 2, 3, hess_rank=10, pca=F, max_iter=1000, kld_ival=2, Y_init=Yi, tol=1e-5)
# print(rr$res_Y)