data(spider)
head(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

#To fit a log-linear model assuming counts are poisson:
glm.spid <- manyglm(spiddat~X, family="poisson")
glm.spid

summary(glm.spid, resamp="residual")

## ALternative data specification as list
dataN <- c(list(spiddat = mvabund(spider$abund)), as.list(data.frame(spider$x)))
dataN 

glm1 <- manyglm(spiddat~soil.dry+bare.sand+fallen.leaves+moss+herb.layer+reflection, data=dataN, family="poisson")
glm1


glm1
glm.spid

summary(glm1, resamp = "residual")


