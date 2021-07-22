

## Only local: 
# setwd("J:/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")
setwd("D:/CD/UEA/gitHRepos/HJA_analyses_Kelpie/Hmsc_CD/oregon_ada")


## trial run 
options(echo=TRUE) # if you want see commands in output file


library(dplyr)
library(raster)
library(sf)
library(ggplot2)
library(cowplot)

# utm10N <- 32610

## on ADA
# gis_out <- gis_in <- "data/gis"

## local
# gis_out <- "J:/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "J:/UEA/Oregon/gis/raw_gis_data"
# gis_out <- "D:/CD/UEA/Oregon/gis/processed_gis_data"
# gis_in <- "D:/CD/UEA/Oregon/gis/raw_gis_data"


baseFolder <- "code_sjSDM/r20210716a"

resFolder <- file.path(baseFolder, "results")
plotsFolder <- file.path(baseFolder, "plots")
# if(!dir.exists(plotsFolder)) dir.create(plotsFolder, recursive = TRUE)
abund <- "pa"

#### Figure paper 


load(file.path(resFolder, "plotData.rdata"))
# save(hsd.df, df2, cut.40, hja.utm, xy.utm, spRich, rStack.sum, spp.in, spp, 

# st_bbox(aoi.pred.sf_edit)
xlim <- c(554650, 572250)
ylim <- c(4890780, 4908610)

head(df2)
## Prepare ordination data

## Load ordination results
load(file.path(resFolder, "ord_tsne_res.rdata")) # 

## function to load site scores into raster
# makeR <- function(r, siteScores, NAs) {
#   
#   rSites <- raster(r)
#   rSites[] <- NA
#   rSites[NAs] <- siteScores
#   rSites
#   
# }
# 
# rSites1 <- makeR(r, tsne$Y[,1], NAs)
# rSites2 <- makeR(r, tsne$Y[,2], NAs)

stck <- raster::stack(rSites1, rSites2)
names(stck) <- c("tsne1", "tsne2")

# prepare for ggplot
coords <- xyFromCell(stck, seq_len(ncell(stck)))
tsne.df <- data.frame(values(stck))
tsne.df <- cbind(coords, tsne.df)
head(tsne.df)

unique(df1$sp.richness)


myTheme <- theme(panel.grid.major = element_line(color = '#cccccc' 
                                      ,linetype = 'dashed'
                                      ,size = .3),
                 axis.title = element_blank(),
                 legend.key = element_rect(fill = "grey80"))

myTheme2 <- theme(panel.grid.major = element_line(color = '#cccccc' 
                                                 ,linetype = 'dashed'
                                                 ,size = .3),
                 axis.title = element_blank(),
                 legend.key = element_rect(fill = "grey80"),
                 axis.text = element_blank())


## species richness
p1 <- ggplot() + 
  geom_tile(data = df1, aes(x, y, fill = sp.richness)) +
  scale_fill_fermenter(palette = "YlGnBu", type = "seq", name = "Species richness")+ ## binned scale
  # scale_fill_gradientn(colours = heat.colors(225), name = "Species richness") +
  # scale_fill_viridis_c(name = "Species richness", option = "C", alpha = 0.7)+
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, aes(color = "A"), bg = NA, inherit.aes = FALSE, show.legend = "polygon")+
  geom_sf(data= xy.utm, aes(color = "B"), size = 2.5, inherit.aes = FALSE, show.legend = "point")+
  geom_sf(data = cut.40, alpha = 0.3, bg = NA, aes(color = "C"), show.legend = "polygon")+
  scale_color_manual(values = c("A" = "blue", "B"= "white", "C" = "grey50"),
                     labels = c("HJA border", "Sample sites", "Logged within 40 years"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(linetype = c("solid", "blank", "solid"), 
                                                              shape = c(NA, 16, NA))))+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
  myTheme

#unlink(file.path(plotsFolder, "results_p1.png"))
ggsave(file.path(plotsFolder, "results_p1.png"), p1, width = 300, height = 300, units= "mm")

## ordination
p2 <- ggplot() + 
  geom_tile(data = tsne.df, aes(x, y, fill = tsne1)) +
  scale_fill_viridis_c(name = "T-SNE axis 1", option = "B", alpha = 0.7)+
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "white", size = 2.5, inherit.aes = FALSE)+
  geom_sf(data = cut.40, col = "grey50", bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
  myTheme

#unlink(file.path(plotsFolder, "results_tsne1.png"))
ggsave(file.path(plotsFolder, "results_tsne1.png"), p2, width = 300, height = 300, units= "mm")

## ordination
p2 <- ggplot() + 
  geom_tile(data = tsne.df, aes(x, y, fill = tsne2)) +
  scale_fill_viridis_c(name = "T-SNE axis 2", option = "B", alpha = 0.7)+
  geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
  scale_alpha(range = c(0.25, 0.65), guide = "none")+
  geom_sf(data = hja.utm, alpha = 0, col = "blue", bg = NA, inherit.aes = FALSE)+
  geom_sf(data= xy.utm, col = "white", size = 2.5, inherit.aes = FALSE)+
  geom_sf(data = cut.40, col = "grey50", bg = NA)+
  coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
  myTheme

#unlink(file.path(plotsFolder, "results_tsne2.png"))
ggsave(file.path(plotsFolder, "results_tsne2.png"), p2, width = 300, height = 300, units= "mm")



# tax groups
# "Diptera"       "Coleoptera"    "Hymenoptera"   "Lepidoptera"

# p3 <- ggplot() + 
#   geom_tile(data = df1, aes(x, y, fill = Diptera)) +
#   scale_fill_viridis_c(name = "Diptera richness", option= "A") +
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "darkblue", bg = NA, inherit.aes = FALSE)+
#   #geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme2
# 
# p4 <- ggplot() + 
#   geom_tile(data = df1, aes(x, y, fill = Coleoptera)) +
#   scale_fill_viridis_c(name = "Coleoptera richness", option= "A") +
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "darkblue", bg = NA, inherit.aes = FALSE)+
#   #geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme2
# 
# 
# p5 <- ggplot() + 
#   geom_tile(data = df1, aes(x, y, fill = Lepidoptera)) +
#   scale_fill_viridis_c(name = "Lepidoptera richness", option= "A") +
#   geom_tile(data = hsd.df, aes(x, y, alpha = hsd), fill = "grey20") +
#   scale_alpha(range = c(0.25, 0.65), guide = "none")+
#   geom_sf(data = hja.utm, alpha = 0, col = "darkblue", bg = NA, inherit.aes = FALSE)+
#   #geom_sf(data= xy.utm, col = "darkred", size = 1, inherit.aes = FALSE)+
#   coord_sf(datum = sf::st_crs(32610), ylim = ylim, xlim = xlim)+
#   myTheme2
# 
# top_row <- plot_grid(p1, p2, labels = c("a", "b"), ncol = 2, label_size = 10)
# bottom_row <- plot_grid(p3, p4, p5, ncol = 3, labels = c("c", "d", "e"), label_size = 10)
# 
# finP <- plot_grid(top_row, bottom_row, ncol = 1)
# ggsave(file.path(plotsFolder, "results_figure.png"), finP, width = 300, height = 300, units= "mm")
# 
# 
# bottom_row <- plot_grid(p3, p4, p5, ncol = 3)
# ggsave(file.path(plotsFolder, "results_tax.png"), bottom_row, width = 300, height = 100, units= "mm")
