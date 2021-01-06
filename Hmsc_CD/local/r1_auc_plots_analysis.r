
sP

### AUC  plots and analysis


## For a given MF or MFCV list

library(Hmsc)
library(dplyr)
library(ggplot2)

wd <- here::here()
wd
setwd(wd)

# get AUC data function
(fn <- list.files(pattern = "^fn_.*\\.r", recursive = TRUE))
source(fn[1])

# get results folders
resF <- list.files("Hmsc_CD/oregon_ada/results", pattern = "res\\d*_\\d{2}$", include.dirs = TRUE, full.names = T)
resF

rf <- resF[length(resF)]
rf

auc.df <- getAUC(rf, rMod = FALSE)
mod.df <- getAUC(rf, rMod = TRUE)

head(auc.df)
head(mod.df)

## how many OTUs per order?
auc.df %>%
  left_join(y = mod.df[, c("modID", "name")]) %>%
  group_by(name, order)%>%
  summarise(noSp = n())%>%
  arrange(name, desc(noSp))

# how many OTUs per family
auc.df %>%
  left_join(y = mod.df[, c("modID", "name")]) %>%
  group_by(name, order, family)%>%
  summarise(noSp = n())%>%
  arrange(name, desc(noSp))


# All orders
auc.df %>%
  left_join(y = mod.df[, c("modID", "name")]) %>%
  #filter(name %in% c("pa_nsp_tp8", "pa_sp_tp8")) %>%
  tidyr::pivot_longer(cols = c("AUC", "TjurR2", "RMSE"), names_to = "metric") %>%
  ggplot(aes(y = value, x = order))+
  geom_boxplot(varwidth = TRUE)+
  facet_wrap(~ name+metric)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Order")+
  ylab("Cross validated evaluation metric")

ggsave("local/eval_metric_All_orders_vif.png")

# Order with > 5 OTUs
auc.df %>%
  left_join(y = mod.df[, c("modID", "name")]) %>%
  # filter(name %in% c("pa_nsp_tp8", "pa_sp_tp8")) %>%
  tidyr::pivot_longer(cols = c("AUC", "TjurR2", "RMSE"), names_to = "metric") %>%
  group_by(name, order, metric)%>%
  filter(n() > 5) %>%
  ungroup() %>%
  ggplot(aes(y = value, x = order))+
  geom_boxplot(varwidth = TRUE)+
  facet_wrap(~ name + metric)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Order")+
  ylab("Cross validated evaluation metric")

getwd()
ggsave("local/eval_metric_gt5OTUs_order_vif.png")
  

