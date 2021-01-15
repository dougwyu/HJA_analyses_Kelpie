#  Functions to plot ecoCopula output (cord())

library(dplyr)
library(ggplot2)
library(corrplot)
library(ggrepel)
library(RColorBrewer)
library(tidygraph)
library(ggraph)
library(glue)

site_res_fxn <- function(cordord, env, xy) {
  data.frame(cordord$scores, env, xy)
}

sp_res_fxn <- function(cordord, abund) {
  data.frame(cordord$loadings, 
             species = colnames(get(glue::glue("otu.{abund}.csv")))) %>%
    separate(species, into = c("OTU", "empty", "class", "order", "family",
                               "genus", "epithet", "BOLD", "BOLDID",
                               "size"), remove = FALSE, sep = "_") %>%
    select(-empty)
}

plot_factors <- function(alphanum, abund, model, spData, siteData, 
                         cont_pred1 = elevation_m, cont_pred2 = lg_YrsDisturb, cat_pred = insideHJA) {
  alpha <- alphanum
  cont_pred1 <- ensym(cont_pred1)
  cont_pred2 <- ensym(cont_pred2)
  cat_pred <- ensym(cat_pred)
  
  
  ggplot() + 
    geom_segment(aes(x = 0, y = 0, 
                     xend = Factor1 * alpha * 0.95, 
                     yend = Factor2 * alpha * 0.95), 
                 data = spData, 
                 size = .1) +
    geom_point(aes(x = Factor1, y = Factor2, 
                   color = !!cont_pred1, 
                   size = !!cont_pred2, 
                   shape = !!cat_pred),
               data = siteData) + 
    # geom_text_repel(aes(x = Factor1, y = Factor2, 
    #                     label = env.csv$SiteName), 
    #                 data = siteData, 
    #                 size = 2) +
    # geom_label_repel(aes(x = Factor1 * alpha, y = Factor2 * alpha, label = species), data = sp_res, size = 2.5) +
    scale_color_gradientn(colours = brewer.pal(n = 10, name = "RdYlBu")) +
    theme_classic() + # or "PuOr" or "RdYlBu"
    labs(title = glue("HJA, ecoCopula ordination, {abund}"), 
         subtitle = model) +
    xlab("Factor 1") +
    ylab("Factor 2")
}

plot_xy_factor1 <- function(abund, model, siteData,
                            cont_pred2 = lg_YrsDisturb, cat_pred = insideHJA) {
  
  cont_pred2 <- ensym(cont_pred2)
  cat_pred <- ensym(cat_pred)
  
  ggplot() + 
    geom_point(aes(x = UTM_E, y = UTM_N, 
                   color = Factor1, 
                   size = !!cont_pred2, # YrsSinceDist, l_rumple
                   shape = !!cat_pred),
               data = siteData) +
    # geom_text_repel(aes(x = UTM_E, y = UTM_N, 
    #                     label = env.csv$SiteName), 
    #                 data = siteData, 
    #                 size = 2) +
    # # geom_label_repel(aes(x = Factor1 * alpha, y = Factor2 * alpha, label = species), data = sp_res, size = 2.5) +
    scale_color_gradientn(colours = brewer.pal(n = 10, name = "PuOr")) +
    theme_classic() + # or "PuOr" or "RdYlBu"
    labs(title = glue("HJA, XY plot, Factor1, {abund}"),
         subtitle = model) +
    xlab("UTM_E") +
    ylab("UTM_N")
}

plot_xy_factor2 <- function(abund, model, siteData,
                            cont_pred2 = lg_YrsDisturb, cat_pred = insideHJA) {
  
  cont_pred2 <- ensym(cont_pred2)
  cat_pred <- ensym(cat_pred)
  
  ggplot() + 
    geom_point(aes(x = UTM_E, y = UTM_N, 
                   color = Factor2, 
                   size = !!cont_pred2, # YrsSinceDist, l_rumple
                   shape = !!cat_pred),
               data = site_res) +
    # geom_text_repel(aes(x = UTM_E, y = UTM_N, 
    #                     label = env.csv$SiteName), 
    #                 data = siteData, 
    #                 size = 2) +
    # geom_label_repel(aes(x = Factor1 * alpha, y = Factor2 * alpha, label = species), data = sp_res, size = 2.5) +
    scale_color_gradientn(colours = brewer.pal(n = 10, name = "PuOr")) +
    theme_classic() + # or "PuOr" or "RdYlBu"
    labs(title = glue("HJA, XY plot, Factor2, {abund}"),
         subtitle = model) +
    xlab("UTM_E") +
    ylab("UTM_N")
}

plot_corrplot <- function(abund, model, siteData) {
  site_res_num <- siteData %>% 
    select(-SiteName, -trap, -period, -lysis_ratio, -COISpike_sum,
           -clearcut, -insideHJA)
  
#   # save corrplot
#   pdf(file <- here("Hmsc_CD/local",
#                   glue("corrplot_{trap}{period}_minimap2_{minimaprundate}_kelpie{kelpierundate}_{abund}.pdf")),
#   width = 11.7, 
#   height = 8.3
#   )
  corrplot(cor(site_res_num),
           method = "ellipse",
           type = "lower",
           title = glue("{abund}, {model}"),
           mar=c(0,0,2,0)
  )
#   dev.off()
#   print(file); rm(file)
}


t_corrplot <- function(mod.cor, title){
  
  corrplot(t(mod.cor), 
           is.corr = F,
           method = "ellipse",
           cl.pos = "n",
           cl.lim = c(-1,1),
           title = title,
           oma = c(0,0,0,0),
           #oma = c(2,2,5,1),
           mar = c(0,0,1,0))
}
  