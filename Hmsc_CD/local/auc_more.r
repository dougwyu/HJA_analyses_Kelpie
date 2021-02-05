## more AUC


# sample data
set.seed(101)
n.spp <- 5
n.occ <- 100
np <- floor(runif(n.occ, 20,60)) # No of presences
na <- n.occ - np # No of absences

# Matrix of species presence / absence
pa <- mapply(function(x,y) rep(c(1,0), c(x, y)), np, na)

# Mean probabilities of presence (p.p) and absence (p.a)
p.p <- c(0.5, 0.8,0.4, 0.7,0.4)
p.a <- c(0.5, 0.6,0.6, 0.4,0.8)
cbind(p.p, p.a)

val <- mapply(function(w,x,y,z) c(rnorm(w, x, 0.1),
                       rnorm(y, z, 0.1)), np, p.p, na, p.a)
       


df1 <- data.frame(species = rep(letters[1:5], each = n.occ),
                  pa = as.vector(pa),
                  value = as.vector(val))

head(df1)

# 5 species, with different patterns of modelled values
library(ggplot2)
ggplot(df1, aes(x = as.factor(pa), y = value))+
geom_boxplot()+
facet_wrap(~species)+
xlab("Presence - Absence")


library(dismo)
## All species
eList <- sapply(split(df1, df1$species), function(x) evaluate(p = x$value[x$pa == 1],a = x$value[x$pa == 0]))
# AUC values
sapply(eList, function(x) x@auc)

par(mfrow = c(2,5))
sapply(split(df1, df1$species), function(x) boxplot(x$value ~ x$pa))
sapply(eList, plot, "ROC")


## AUC is wilcoxon
auctest <- function(e) {
  w <- wilcox.test(e@presence, e@absence)
  pauc <- w$p.value
  auc <- as.vector(w$statistic) / (e@na * e@np)
  return(auc)
}

sapply(eList, auctest)

# PROC - flips the values if direction not specified
sapply(split(df1, df1$species), function(x) pROC::roc(data = x, response = "pa", predictor = value)$auc)


png("Hmsc_CD/local/plots/auc_roc_plot.png", width = 200, height = 200, units = "mm", res = 100)
par(mfrow = c(4,5))
sapply(split(df1, df1$species), function(x) boxplot(x$value ~ x$pa, xlab = "", ylab = "predicted prob"))
sapply(eList, plot, "ROC")
sapply(split(df1, df1$species), function(x) plot(pROC::roc(data = x, response = "pa", predictor = value)))
mtext("pROC direction: 'auto'", 3, 2.5, adj = 1)
sapply(split(df1, df1$species), function(x) plot(pROC::roc(data = x, response = "pa", predictor = value, direction = "<")))
mtext("pROC Direction: '<'", 3, 2.5, adj = 1)
dev.off()
