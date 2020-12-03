

## Load data for Oregon HMsc analysis ####

source("https://raw.githubusercontent.com/Cdevenish/R-Material/master/Functions/w.xls.r")

## Data from github.. for moment
dataF <- "https://raw.githubusercontent.com/dougwyu/HJA_analyses_Kelpie/master/sjSDM/R-git/data/kelpie_data/for_adagpu/data_20201024_5minocc_gismslidar"

XY <- read.csv(file.path(dataF, "XY.csv"))
pa <- read.csv(file.path(dataF, "otu.pa.csv"))
qp <- read.csv(file.path(dataF, "otu.qp.csv"))
env <- read.csv(file.path(dataF, "scale.env.csv"))

head(XY)
pa[1:10, 1:5]
qp[1:10, 1:5]
head(env)

spp <- data.frame(col = colnames(pa))
spp$spID <- sprintf("sp%03d", 1:nrow(spp))
spp$class <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\1", spp$col)
spp$order <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\2", spp$col)
spp$family <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\3", spp$col)
spp$genus <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\4", spp$col)
spp$epiphet <- sub(".*__([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*)_([[:alpha:]]*).*", "\\5", spp$col)
head(spp)

spp <- spp[order(spp$class, spp$order, spp$family, spp$genus),]
tax.cols <- c("class", "order", "family", "genus", "epiphet", "spID")
for(i in tax.cols) spp[,i] <- factor(spp[,i])
str(spp)

head(spp, 15)
sum(is.na(spp$genus))

spp.tree <- ape::as.phylo(~class/order/family/genus/spID, data = spp, collapse = F)

plot(spp.tree)
spp.tree$edge.length = rep(1, length(spp.tree$edge)) # make all lengths eqaul between tax levels
plot(spp.tree)
ape::is.rooted(spp.tree)

dist <- ape::cophenetic.phylo(spp.tree)
dist # distances are just number of taxonomic levels apart.... ie it doesn't matter whether families, genera, etc are 
## in correct order.

sum(dist)

# w.xls(dist, row = T)
# ape::write.tree(spp.tree, "P.tre")



