### RESPONSE CURVES

## Gets data from b1_pred_response_curves.r

#### Read data on Ada  #####

## Only testing local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
wd <- here::here()
wd
setwd(file.path(wd, "Hmsc_CD"))
dir()

getwd()

resFolder <-"oregon_ada/code_sjSDM/r20210311a/results"
load(file.path(resFolder, "responseData.rdata"))
# modFull_nsp, modFull_sp, predList_sp, predList_nsp


## each list element is one variable
## inside is matrix, col1 is focal variable on unscaled. cols 2:, are species columns

## Plot all species, and means.. 
vars <- names(predList_nsp)
spp <- colnames(predList_nsp$be10[,2:ncol(predList_nsp$be10)])

ylim <- range(predList_nsp$be10[,2:ncol(predList_nsp$be10)])

# plot(predList_nsp$be10[,1], predList_nsp$be10[,2], type = "l", ylim = ylim, col = "grey80")
# for(sp in spp) points(predList_nsp$be10[,1],predList_nsp$be10[,sp], col = "grey80", lty = 1, type = "l")
# points(predList_nsp$be10[,1], rowMeans(predList_nsp$be10[,2:ncol(predList_nsp$be10)]), col = "black", lwd = 2, type = "l")

matplot(predList_nsp$be10[,1], predList_nsp$be10[,2:ncol(predList_nsp$be10)], type = "l", 
      col = "grey80", lty = 1)
points(predList_nsp$be10[,1], rowMeans(predList_nsp$be10[,2:ncol(predList_nsp$be10)]), col = "black", lwd = 2, type = "l")

### rearrange data ...
str(predList_nsp, max.level = 1)

# and join with species incidence
load(file.path(resFolder, "modelData.rdata"))
rm(otu.pa.csv, otu.qp.csv, otuenv, env.vars, spChoose, k, minocc, noSteps, vars, varsName, abund, device, iter, sampling)

# x <- predList_sp[[1]]

predRes <- lapply(predList_nsp, function(x){ ## CHANGE HERE FOR spatial / non-spatial
# predRes <- lapply(predList_sp, function(x){ ## CHANGE HERE FOR spatial / non-spatial
  
  
 tmp <- data.frame(x) %>%
   tidyr::pivot_longer(cols = contains("__"), 
                       names_to = "OTU", 
                       values_to = "pred") %>%
   tidyr::pivot_longer(col = 1, names_to = "var", values_to = "grad") %>%
   dplyr::left_join(y = spIncid, by = "OTU") %>%
   tidyr::separate(col = OTU, into = c("OTU", "empty", "class", "order", "family","genus", 
                                             "epithet", "BOLD", "BOLDID","size"),
                     remove = FALSE, sep = "_") %>%
     dplyr::select(-c(empty, BOLD, BOLDID, size)) %>%
   arrange(class, order, family, genus, grad) # ascending by grad for plotting - not necessary
   
 # head(tmp)
 # sum(is.na(tmp$incidence))
  
 
  
}
)

nsp.df <- do.call(rbind, predRes)
head(nsp.df)

222*250*42
length(unique(nsp.df$OTU))
unique(nsp.df$class)

tapply(nsp.df$OTU, nsp.df$class, function(x) length(unique(x)))

library(ggplot2)
library(ggforce)

## Get means for plotting on top of ind spp lines
varMeans <- nsp.df %>%
  group_by(var, grad) %>%
  summarise(pred = mean(pred))

head(varMeans)


# i = 1

length(vars) / 9

pdf("local/plots/nsp_response_plot_incidence.pdf")
#pdf("Hmsc_CD/local/plots/sp_response_plot_incidence.pdf")
for(i in 1:4){
  print(nsp.df %>%
          ggplot(aes(x = grad, y = pred, col = incidence))+ # col = order
          geom_line(aes(group = OTU), alpha = 0.3, size = 0.5) + # col = "grey80"
          geom_line(data = varMeans, col = "black", size = 1) +
          scale_color_viridis_c(option = "B", direction = -1) +  # D is viridis
          facet_wrap_paginate(~var, scales = "free_x", nrow= 3, ncol = 4, page = i) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ggtitle("Partial response per species + mean - pa min6 vars2 DNN")+
          ylab("Prediction")+
          xlab(""))
}
dev.off()


pdf("local/plots/nsp_response_plot_order.pdf")
# pdf("Hmsc_CD/local/plots/sp_response_plot_order.pdf")
for(i in 1:4){
  print(nsp.df %>%
          ggplot(aes(x = grad, y = pred, col = order))+ # col = order
          geom_line(aes(group = OTU), alpha = 0.3, size = 0.5) + # col = "grey80"
          geom_line(data = varMeans, col = "black", size = 1) +
          #scale_color_viridis_d(option = "A") + 
          facet_wrap_paginate(~var, scales = "free_x", nrow= 3, ncol = 4, page = i) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ggtitle("Partial response per species + mean - pa min6 vars2 DNN")+
          ylab("Prediction")+
          xlab(""))
}
dev.off()

pdf("local/plots/nsp_response_plot_grey.pdf")
# pdf("Hmsc_CD/local/plots/sp_response_plot_grey.pdf")
for(i in 1:4){
  print(nsp.df %>%
          ggplot(aes(x = grad, y = pred))+ # col = order
          geom_line(aes(group = OTU), alpha = 0.3, size = 0.5, col = "grey80") + # 
          geom_line(data = varMeans, col = "black", size = 1) +
          scale_color_viridis_d(option = "A") + 
          facet_wrap_paginate(~var, scales = "free_x", nrow= 3, ncol = 4, page = i) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ggtitle("Partial response per species + mean - pa min6 vars2 DNN")+
          ylab("Prediction")+
          xlab(""))
}
dev.off()


##### Get range of responses ... as sort of measure of magnitude of response. 
## could do a more sophisticated wiggliness of curve, etc. but probably similar

head(nsp.df)

nsp_resp_range <- nsp.df %>%
  group_by(var, OTU, order, genus, epithet, incidence) %>%
  summarise(mean = mean(pred),
            range = diff(range(pred)),
            IQR = IQR(pred))
  
## order vars by IQR
sortVar <- nsp_resp_range %>%
  group_by(var)%>%
  summarise(sumIQR = sum(IQR)) %>%
  arrange(desc(sumIQR))
  
nsp_resp_range$var <- factor(nsp_resp_range$var, levels = sortVar$var)

head(nsp_resp_range)

pdf("local/plots/nsp_response_plot_IQR_incidence.pdf", width = 10, height = 6)
# pdf("Hmsc_CD/local/plots/sp_response_plot_grey.pdf")
ggplot(nsp_resp_range, aes(x = var, y = range))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col= incidence), alpha = 0.3, size = 0.5)+
  scale_color_viridis_c(option = "B", direction = -1) +  # D is viridis  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Range of partial response, ordered by IQR - pa min6 vars2 DNN")
dev.off()


pdf("local/plots/nsp_response_plot_IQR.pdf", width = 10, height = 6)
# pdf("Hmsc_CD/local/plots/sp_response_plot_grey.pdf")
ggplot(nsp_resp_range, aes(x = var, y = range))+
  geom_boxplot()+
  #geom_jitter(aes(col= incidence), alpha = 0.4, size = 0.75)+
  # scale_color_viridis_c(option = "B", direction = -1) +  # D is viridis  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Range of partial response, ordered by IQR - pa min6 vars2 DNN")

dev.off()




ggplot(nsp_resp_range, aes(x = var, y = range))+
  geom_boxplot()+
  coord_flip()

## how does range (size of response) vary with incidence and mean

head(nsp_resp_range)

ggplot(nsp_resp_range, aes(x = range, y = incidence, col = var))+
  geom_point()
  

ggplot(nsp_resp_range, aes(x = mean, y = incidence, col = var))+
  geom_point()

ggplot(nsp_resp_range, aes(x = range, y = incidence, col = var))+
  geom_point()



##  are some speices always the ones to repsonse strongly to all variables.. or is more unviform?
head(nsp.df)
head(nsp_resp_range)

rankSp <- nsp_resp_range %>%
  group_by(var) %>%
  mutate(rank = rank(range),
         spp = paste(genus, epithet))%>%
  arrange(rank, OTU)
  
rankSp

spF <- rankSp %>%
  group_by(OTU, order, genus, epithet)%>%
  summarise(sumRank = sum(rank)) %>%
  arrange(sumRank)

spF
rankSp$OTU <- factor(rankSp$OTU, levels = spF$OTU)

rankSp


rankSp %>%
  group_by(OTU, order, genus, epithet)%>%
  summarise(sumRank = sum(rank)) %>%
  mutate(spp = paste(genus, epithet)) %>%
  ggplot(aes(y = sort(sumRank), x = 1:length(sumRank), col = order)) + # 
  geom_point()


rankSp %>%
  group_by(OTU, order, genus, epithet)%>%
  summarise(sumRank = sum(rank)) %>%
  arrange(sumRank)
  
pdf("local/plots/nsp_response_range_rank_x_species.pdf", width = 12, height = 6)
rankSp %>%
  group_by(OTU, order, spp)%>%
  ggplot(aes(x = OTU, y = rank, col= order))+
  geom_boxplot()+
  scale_x_discrete(labels = paste(spF$genus, spF$epithet))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))+
  ylab("Ranks of species response range across variables")
  
dev.off()
