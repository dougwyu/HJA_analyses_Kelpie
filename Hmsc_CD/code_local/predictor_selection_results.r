## Predictor selection  #####

## load univariate and all but one models and evaluate predictor contributions

## Only local: 
setwd("J:/UEA/Oregon/Oregon_winscp")
dir()

# models from S0_pSelection_univariate.r and S0_pSelection_allBut.r

thin <- 10
samples <- 250

# Univariate
load(paste0("results/predSelection/prRes_univariate_thin", thin, "_samp", samples, ".rdata"))
# user  system elapsed 
# 2.582   0.382 197.759  - same time as job email log 
# 1095326 Name=Hmsc_pred_uni Ended, Run time 00:03:21
# 50 samples, about 1.5 mins
# Job_id=1142366 Name=Hmsc_pred_uni_250 Ended, Run time 00:16:03
# with CV, samples 250, thin 10,  Job_id=1144277 Name=Hmsc_pred_uni_10 Ended, Run time 00:46:20

thin <- 5
samples <- 50
# All But one
load(paste0("results/predSelection/prRes_allB_thin", thin, "_samp", samples, ".rdata"))
# user    system   elapsed 
# 2.974     0.444 26780.803  - same as below..
# id=1096121 Name=Hmsc_pred_uni Ended, Run time 07:26:25

rm(thin, samples)


uni <- unlist(prRes, recursive = F)
# str(uni, max.level = 1)
names(uni)
MF.uni <- uni[grepl("MF", names(uni))]
beta.uni <- uni[grepl("beta", names(uni))]

allB <- unlist(prRes_allB, recursive = F)
MF.allB <- allB[grepl("MF", names(allB))]
beta.allB <- allB[grepl("beta", names(allB))]


### Beta convergence
# sapply(beta.uni, mean)
# sapply(beta.allB, mean)
sapply(beta.uni, function(x) sprintf("%.3f \U00B1 %.3f",mean(x), sd(x)))
sapply(beta.allB, function(x) sprintf("%.3f \U00B1 %.3f",mean(x), sd(x)))

# boxplot(beta.uni, las = 2)
par(mfrow = c(2,1), mar = c(10,2,1,1))
vioplot::vioplot(beta.uni, las = 2, names = sub("beta_(.*)_uni", "\\1", names(beta.uni)), main = "univariate")
vioplot::vioplot(beta.allB, las = 2,names = sub("beta_(.*)_all_but", "\\1", names(beta.allB)), main = "all but one")
par(op)


## Predictor evaluation

mfRes.uni <- sapply(MF.uni, function(x) sapply(x, mean))
mfRes.allB <- sapply(MF.allB, function(x) sapply(x, mean))

mfRes.uni
mfRes.allB

# ## Bar plot of all values, sorted by AUC
# op <- par(mar = c(9,5,4,1))
# nm <- sub("MF_(.*)_uni", "\\1", colnames(mfRes.uni[,order(mfRes.uni["AUC",], decreasing = TRUE)]))
# barplot(mfRes.uni[,order(mfRes.uni["AUC",], decreasing = TRUE)], 
#         beside = T, las = 2, names.arg = nm, legend = T, args.legend = list(x = "top", horiz = T, inset = -0.2))
# par(op)

# sort by AUC
# sort(mfRes[2,], decreasing = TRUE)

# sort by RMSE
# nm.rmse <- sub("MF_(.*)_uni", "\\1", colnames(mfRes[,order(mfRes["RMSE",])]))
# plot(sort(mfRes[1,]), pch = 16)
# text(x = 1:ncol(mfRes), y = sort(mfRes[1,]), labels = nm.rmse, srt = 90, cex = 0.7, pos =3)
# sort(mfRes[1,])


df.uni <- data.frame(mfRes.uni)
df.uni$metric <- rownames(df.uni)
rownames(df.uni) <- NULL
df.uni

df.allB <- data.frame(mfRes.allB)
df.allB$metric <- rownames(df.allB)
rownames(df.allB) <- NULL
df.allB


all.res <- rbind(tidyr::pivot_longer(df.uni, cols = !metric, 
                                     names_pattern = "(MF|MFCV)_(.*)_(.*$)", 
                                     names_to = c("eval", "predictor", "type")),
                 tidyr::pivot_longer(df.allB, cols = !metric, 
                                     names_pattern = "(MF|MFCV)_(.*)_(all_but$)", 
                                     names_to = c("eval", "predictor", "type")))

rm(df.uni, df.allB)

head(all.res)
tail(all.res)

library(ggplot2)
library(dplyr)

# order AUC descending for CV univariate
pred.uni.CV <- all.res %>%
  filter(metric == "AUC" & type == "uni" & eval == "MFCV") %>%
  arrange(desc(value))
  

all.res %>%
  filter(metric == "AUC" & type == "uni") %>%
  mutate(pred = factor(predictor, levels = pred.uni.CV$predictor)) %>%
  #arrange(type, desc(value)) %>%
  ggplot(aes(x = pred, y = value, fill=eval))+
  geom_bar(stat="identity", position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_discrete(name = "Univariate", labels = c("Explanatory", "Predictive (CV)"))+
  ylab("AUC")+
  xlab("Predictor")

ggsave("results/predSelection/univariate_AUC.png")


# order AUC for CV  - jacknife
pred.allB.CV <- all.res %>%
  filter(metric == "AUC" & type == "all_but" & eval == "MF") %>%
  arrange(value)


all.res %>%
  filter(metric == "AUC" & type == "all_but") %>%
  mutate(pred = factor(predictor, levels = pred.allB.CV$predictor)) %>%
  #arrange(type, desc(value)) %>%
  ggplot(aes(x = pred, y = value, fill=eval))+
  geom_bar(stat="identity", position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_discrete(name = "All but one", labels = c("Explanatory", "Predictive (CV)"))+
  ylab("AUC")+
  xlab("Predictor")


ggsave("results/predSelection/jacknife_AUC.png")
