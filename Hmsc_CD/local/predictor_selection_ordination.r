##### Ordinations

here::i_am("local/predictor_selection_ordination.r")
wd <- here::here()
setwd(wd)

## Only local:
getwd()
# setwd("J:/UEA/Oregon/Oproject")
dir()

library(vegan)

load(file = file.path("ada/data", "allData.Rdata")) # S,X,Y,P & Tax (where P is based on Tax)

raretaxa <- which(colSums(Y.train > 0) < 10)
length(raretaxa)

Y.train <- Y.train[, -raretaxa]
Y.train[1:5, 1:5]

try = 100
trymax = 5000
parallel = 3

sum(no.shared(Y.train))

nmds <- vegan::metaMDS(Y.train, distance = "bray", k = 2, 
                       try = try, trymax = trymax, parallel = parallel, binary = TRUE)

nmds

save(nmds, file = "results/predSelection/nmds.rdata")
load("results/predSelection/nmds.rdata")

str(X.train)
envFit <- vegan::envfit(nmds, X.train)
envFit
# 
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# oldGrowthIndex   -0.66983 -0.74252 0.0411  0.179    
# elevation_m      -0.97893 -0.20422 0.2048  0.001 ***
#   canopyHeight_m    0.08737 -0.99618 0.0501  0.093 .  
# precipitation_mm -0.69913 -0.71499 0.1783  0.001 ***
#   distToStream_m    0.16777  0.98583 0.0347  0.213    
# mean.NDVI        -0.36583 -0.93068 0.1805  0.001 ***
#   mean.EVI          0.23659 -0.97161 0.0950  0.018 *  
#   mean.green        0.98925  0.14622 0.0086  0.705    
# mean.wet         -0.30160 -0.95344 0.1919  0.001 ***
#   l_p25            -0.40298 -0.91521 0.0983  0.013 *  
#   l_rumple         -0.12950 -0.99158 0.0631  0.071 .  
# B1_mean           0.80066  0.59912 0.3160  0.001 ***
#   B4_mean           0.46272  0.88651 0.2160  0.001 ***
#   lg_DistRoad      -0.96074  0.27746 0.1792  0.001 ***
#   lg_YrsDisturb    -0.14343 -0.98966 0.1834  0.001 ***
#   lg_cover2m_max   -0.12954 -0.99157 0.1822  0.001 ***
#   lg_cover2m_4m     0.41353  0.91049 0.1421  0.003 ** 
#   lg_cover4m_16m    0.83093 -0.55637 0.0041  0.835    
# ---

#               r2      Pr(>r)    
# clearcut     0.0237  0.145    
# insideHJA    0.1552  0.001 ***
# YrsDist_gt00 0.0311  0.072 .  

pdf("results/predSelection/nmds_plot.pdf")
plot(nmds)
plot(envFit, p.max = 0.05)
# help see which labels are overlapped
#rownames(envFit$vectors$arrows)[order(envFit$vectors$arrows[,1])]
dev.off()

ordisurf(nmds, X.train$elevation_m) # default displays sites - open circles
ordisurf(nmds, X.train$distToStream_m, main = "Dist to stream")

pdf("results/predSelection/ordisurf%02d.pdf")
par(mfrow = c(3,2))
for(x in colnames(X.train)[!sapply(X.train, is.factor)]) ordisurf(nmds, X.train[,x], main = x)
dev.off()


## Compare to rda
rda <- vegan::rda(X = Y.train, Y = X.train)
plot(rda)
