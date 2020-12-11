
here::i_am("local/ecocopula_analysis.r")
wd <- here::here()
wd
setwd(wd)

## Ecocopula analysis and plots #####

library(mvabund)
library(ecoCopula)
library(ggplot2)

## get data
load(file.path(wd, "oregon_ada/data", "allData.Rdata")) # SXY.train, S.train, X.train, Y.train, P

# predictors
str(X.train)
pred <- data.frame(lapply(X.train, as.numeric)) # scale.. here.. 
head(pred)

corrplot::corrplot(cor(pred), method = "ellipse", type = "lower")


Y.train[1:3, 1:5]
str(Y.train)
head(S.train)

# make pa matrix
Y.pa <- mvabund(Y.train)
is.mvabund(Y.pa)
str(Y.pa)

# mv glm with no covariates... 
mvBin <- mvabund::manyglm(Y.pa ~ 1, family = "binomial")
plot(mvBin) # chk residuals
##### summary(mvBin) ## causing problems.. crashing R studio.. takes a while to run... 

# do ordination
Y.ord <- cord(mvBin)

# Site scores: join ordination axes to data frame of predictors, coordinates
site_res <- data.frame(Y.ord$scores, pred, S.train[, c("UTM_E", "UTM_N")])
head(site_res)

# Species scores
sp_res <- data.frame(Y.ord$loadings, species = colnames(Y.train))
str(sp_res)

corrplot::corrplot(cor(site_res), method = "ellipse", type = "lower")



plot(site_res[, c("Factor1", "Factor2")], pch = 16, col = "darkblue")
points(sp_res[, c("Factor1", "Factor2")], pch = 17)

sF <- 3

ggplot()+
  geom_point(aes(x = Factor1, y = Factor2, 
                  color = elevation_m,
                 shape = as.factor(insideHJA)),
             data = site_res) + 
  geom_point(aes(x = Factor1*sF, y = Factor2*sF),
                 data = sp_res)+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu"))



alpha <- 4
p0 <- ggplot() + 
  geom_segment(aes(x = 0, y = 0, 
                  xend = Factor1 * alpha * 0.95, 
                   yend = Factor2 * alpha * 0.95), 
               data = sp_res, 
               size = .1) +
  geom_point(aes(x = Factor1, y = Factor2, 
                 color = elevation_m, 
                 size = lg_YrsDisturb, 
                 shape = as.factor(insideHJA)
  ), 
  data = site_res) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")) +
  theme_classic() + # or "PuOr" or "RdYlBu"
  labs(title = "HJA ecoCopula ordination pa Y train")
p0

