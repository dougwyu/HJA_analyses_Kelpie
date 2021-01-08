
## Extract and display results by model

getwd()

library(Hmsc)
wd <- here::here()
wd
setwd(wd)

source("Hmsc_CD/local/fn_getAUC.r")

## load results
# get results folders 
resF <- list.files("Hmsc_CD/oregon_ada/results", pattern = "res\\d*_\\d{2}$", include.dirs = TRUE, full.names = T)
# resF <- c(resF, list.files("Hmsc_CD/oregon_ada/results", pattern = "RRR_test2", include.dirs = TRUE, full.names = T))

resF

# subset results
rF <- resF[!resF %in% c("Hmsc_CD/oregon_ada/results/res20210107_01",
                          "Hmsc_CD/oregon_ada/results/res20201217_01",
                        "Hmsc_CD/oregon_ada/results/res20201127_01",
                        "Hmsc_CD/oregon_ada/results/res20201204_01", # qp 
                        "Hmsc_CD/oregon_ada/results/res20201209_01", # qp
                        "Hmsc_CD/oregon_ada/results/res20201216_02" # error
                        )]

## Loop through results folders, extract model definition, convergence results and evaluation
modRes <- lapply(rF, getAUC, rMod = TRUE)

all.df <- do.call(rbind, modRes)

head(all.df)
rm(modRes, resF, rF)

all.df <- all.df[order(all.df$AUC_pred, decreasing = TRUE),]
head(all.df)
tail(all.df)

getwd()
# save(all.df, all.me, file = "local/modelResults.rdata")


### TO HERE ON THIS VERSION 2....  moved above block to rmd for table.


## Relationship of AUC with prevalence
library(tidyr)
library(ggplot2)

me.l <- pivot_longer(all.me, cols = colnames(all.me)[grepl("_expl|_pred", colnames(all.me))], 
                            names_to = c("metric", "type"), names_sep ="_")
# ends_with("_expl", "_pred")
head(me.l)
tail(me.l)
unique(me.l$type)

unique(me.l$rf)

res <- c("res20201202_01", "res20201204_01", "res20201209_01", "RRR_test2")

me.l %>%
  filter(!grepl("RRR|adj_evi|q", name) & rf %in% res & type == "pred") %>%
  ggplot(aes(x = prev, y = value))+
  geom_point()+
  facet_grid(rows = vars(metric), cols = vars(name), scales = "free_y")

getwd()
ggsave("local/prevPlot.png")  

ggplot(subset(me.l, rf %in% res & type == "pred"), aes(x = prev, y = value))+
  geom_point()+
  facet_wrap(~ metric + name)


m <- models[[2]]
name <- names(models)[2]
m
name

# Get order of species prevalence - no of sites present at
sp.prev <- colSums(m$Y)/nrow(m$Y) # pa data

str(MFCV, max.level = 1)
mfcv <- MFCV[[paste(name, thin, samples, sep = "_")]]

plot(mfcv$AUC ~ sp.prev)
## lower limit of AUC little bit explained by prevalence.... ???

plot(mfcv$TjurR2 ~ sp.prev)
plot(mfcv$RMSE ~ sp.prev)

## Top predictor coefficients for highest AUC,... highest prevalence???

## Species associations 
# residual associations (taking into account covariates)

# sort species names by prevalence
which.max(sp.prev)
which.min(sp.prev)

sort(sp.prev)[c(1, length(sp.prev))]
prev.sort <- sort(sp.prev)

omegaCor <- computeAssociations(m)
str(omegaCor, max.level = 2)

# support level to show positive/negative interactions
sLevel <- 0.95
sLevel <- 0.99

toPlot <- ((omegaCor[[1]]$support > sLevel)
           + (omegaCor[[1]]$support < (1-sLevel)) > 0) * omegaCor[[1]]$mean

# str(toPlot)
# rearrange toPlot in order of species prevalence, rarest to most common
dimnames(toPlot)

head(dimnames(toPlot)[[1]])
prev.names[1]

toPlot_prev <- toPlot[match(names(prev.sort), dimnames(toPlot)[[1]]), match(names(prev.sort), dimnames(toPlot)[[2]])]
toPlot_prev[1:2, 1:2]
tmp <- (length(sp.prev)-1):length(sp.prev)
toPlot_prev[tmp, tmp]

labels <-attr(toPlot_prev, "dimnames")[[1]]
tx.order <- sub(".*__[[:alpha:]]*_([[:alpha:]]*)_.*", "\\1", labels)
tx.order <- as.factor(tx.order)
levels(tx.order)


# assign new labels to cor matrix
attr(toPlot_prev, "dimnames") <- list(substr(tx.order, 1,3), substr(tx.order, 1,3))

# No labels
#pdf(file.path(fp,"spCorPlot.pdf"))
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "n", mar = c(1,1,1,1),
                   tl.cex = 0.4)
# dev.off()

# tax order labels
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(1,1,1,1), tl.col = rainbow(length(unique(tx.order)))[tx.order], 
                   tl.cex = 0.4)

# Show quantiles of rarity/commonness with labels
quantile(prev.sort)
qLabs <- as.numeric(cut(prev.sort, quantile(prev.sort), include.lowest = T, labels = 1:4))
table(qLabs, useNA = "always")

print("\U25A0")
attr(toPlot_prev, "dimnames") <- list(rep("\U25A0", dim(toPlot_prev)[1]), rep("\U25A0", dim(toPlot_prev)[1]))

# show rarity quantiles
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(1,1,1,1), tl.col = c("grey", "black", "grey", "black")[qLabs], 
                   tl.cex = 0.45)


pdf(file.path(fp,"spCorPlot.pdf"))
corrplot::corrplot(toPlot_prev, method = "color", col = c("red", "white", "blue"),
                   cl.pos = "r", tl.pos = "lt", mar = c(2,1,2,1), tl.col = c("grey", "black", "grey", "black")[qLabs], 
                   tl.cex = 0.45)
dev.off()

## Use latent variables to show species associations

etaPost = getPostEstimate(m, "Eta") # sites, circles, latent factor
lambdaPost = getPostEstimate(m, "Lambda") # species, triangles, latent factors

# how many latent factors are there? decided in model..
dim(etaPost$mean) # 88, 5 factors.

# as with ordination, first factors will have most information... 

m$XFormula
Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = "lg_YrsDisturb", main = name, spNames = NA) 
# blue lowest values, red, highest values.  colorRampPalette(c("blue","white","red"))

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = , main = name, spNames = NA) 


Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes
str(m$XData)

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
             colVar = "oldGrowthIndex", main = name, spNames = NA) # numeric, low to high(red)



# 
Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,3),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes

Hmsc::biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(3, 4),
             colVar = "insideHJA", main = name, spNames = NA,
             colors = c("blue", "red")) # levels' no, yes

source(file.path(wd, "local","biplot_hmsc.r"))
# modified biplot function to allow species colors... eg by taxonomy or trait.. and removed par setting asp... see
# if species patterns can be seen better....

biPlot_hmsc(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,2),
            colVar = "insideHJA", main = name, spNames = NA,
            colors = "black", spCols = rainbow(length(unique(tx.order)))[tx.order], asp = 1) # levels' no, yes

legend("top", ncol = 2, bty = "n", cex = 0.7, legend = unique(tx.order), pch = 17, col = rainbow(length(unique(tx.order))))
