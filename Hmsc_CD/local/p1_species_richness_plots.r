### Species richness plots comparing high and low prevalence species

wd <- here::here()
wd
setwd(wd)
dir()

# (follows quantiles of prevalence - as in models of 20201209 - see S2_define_models.r)

# get data
source("Hmsc_CD/local/L1_read_data.r")

# otu.pa.csv is Y.train.pa; otu.qp.csv is Y.train.qp

## Make data sets with  highest species prevalence
sp.prev <- colSums(otu.pa.csv)/nrow(otu.pa.csv)
hist(sp.prev)

# cut into quantiles
quantile(sp.prev)

# no of species in each quantile
table(cut(sp.prev, quantile(sp.prev), include.lowest = T, right = F), useNA= "always") # closed on left

# make groups
qs <- c("[0% 25%)", "[25% 50%)", "[50% 75%)", "[75% 100%)") # labels
qGrps <- cut(sp.prev, quantile(sp.prev), include.lowest = T, right = F, labels = qs)
table(qGrps, useNA = "always")

# Species richness across sites for all species.
hist(rowSums(otu.pa.csv))

## histogram of each quantile
png("Hmsc_CD/local/plots/spRich_histograms.png", width = 200, height = 200, units= "mm", res = 200)
par(mfrow = c(2,2))
for(i in qs) hist(rowSums(otu.pa.csv[, qGrps == i]), main = paste("quantile",i), xlab = "Sp richness")
dev.off()

# scatter plot of lowest prevalence species against all other quantiles
png("Hmsc_CD/local/plots/spRich_correlations.png", width = 200, height = 200, units= "mm", res = 200)
par(mfrow = c(2,2))
for(i in qs) plot(rowSums(otu.pa.csv[, qGrps == i]),rowSums(otu.pa.csv[, qGrps == qs[1]]), 
                  ylab = "sp rich quantile 1", xlab = paste("sp rich quantile",i), pch = 16)
dev.off()

# correlation test
for(i in qs) print(cor.test(rowSums(otu.pa.csv[, qGrps == i]),rowSums(otu.pa.csv[, qGrps == qs[1]]), method = "spearman"))


# Spearman's rank correlation rho
# 
# data:  rowSums(otu.pa.csv[, qGrps == i]) and rowSums(otu.pa.csv[, qGrps == qs[1]])
# S = 0, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
#   1 
# 
# 
# 	Spearman's rank correlation rho
# 
# data:  rowSums(otu.pa.csv[, qGrps == i]) and rowSums(otu.pa.csv[, qGrps == qs[1]])
# S = 88287, p-value = 0.03713
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2225762 
# 
# 
# Spearman's rank correlation rho
# 
# data:  rowSums(otu.pa.csv[, qGrps == i]) and rowSums(otu.pa.csv[, qGrps == qs[1]])
# S = 103923, p-value = 0.4316
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.08489717 
# 
# 
# 	Spearman's rank correlation rho
# 
# data:  rowSums(otu.pa.csv[, qGrps == i]) and rowSums(otu.pa.csv[, qGrps == qs[1]])
# S = 94911, p-value = 0.1262
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.1642518
