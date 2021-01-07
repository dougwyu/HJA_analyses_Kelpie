
## Ecocopula plots ######

library(ggplot2)
library(cowplot)

## load ecocopula results
load("Hmsc_CD/local/ecocopula_modList_pilot.rdata") # fm, cor.preds, modList

str(modList, max.level =2)

##### Plots
source("Hmsc_CD/local/ecoCopula_plot_fn.r")## modified plot functions...

mod <- modList[[1]]$mod
sp_res <- modList[[1]]$sp
site_res <- modList[[1]]$site

alpha <- 4*0.95

# chk residuals
plot(mod) 

#plot factors
plot_factors(alpha, "pa", "mod1", sp_res, site_res)

# Plot over coords with factors as colour fill
plot_xy_factor1("pa", "mod1", site_res)
plot_xy_factor2("pa", "mod1", site_res)

## correlation plot
cor.preds <- colnames(env.vars)[sapply(env.vars, is.numeric)]
mod.cor <- cor(mod.ord$scores, site_res[,cor.preds])
t_corrplot(mod.cor, "mod1")


## loop over list and do all plots

# i = 1
alpha <- 4*0.95

for(i in seq_along(modList)){
  

  mod <- modList[[i]]$mod
  mod.ord <- modList[[i]]$ord
  sp_res <- modList[[i]]$sp
  site_res <- modList[[i]]$site
  

  # chk residuals
  # plot(mod) 
  
  modName <- paste0("predictors: ", as.character(fm[[i]])[3])
  
  #plot factors
  p1 <- plot_factors(alpha, "pa", "", sp_res, site_res)
  
  
  ## correlation plot
  mod.cor <- cor(mod.ord$scores, site_res[,cor.preds])
  p2 <- ~t_corrplot(mod.cor, "")
  
  # Plot over coords with factors as colour fill
  p3 <- plot_xy_factor1("pa", "", site_res)
  p4 <- plot_xy_factor2("pa", "", site_res)
  
  # put rows together
  pA <- plot_grid(p2, p1, rel_widths = c(4,5), align = "h", axis = "t")
  pB <- plot_grid(p3, p4)
  
  title <- ggdraw() + 
    draw_label(modName, fontface = 'bold', size = 18, x = 0.3, y = 0.5, hjust = 0)
    
  m1 <- plot_grid(pA, pB, title, nrow = 3, rel_heights = c(4,5,1))
  # m1 <- plot_grid(pA, pB, nrow = 2, rel_heights = c(4,5))
  # m1
  ggsave(paste0("Hmsc_CD/local/plots/mod", i, ".png"), m1, width = 300, height = 250, units = "mm")
}





# ggsave("Hmsc_CD/local/factor1_8predictors.png")
# ggsave("Hmsc_CD/local/factor2_8predictors.png")


# # get correlations with selected predictors
# (mod02.cor <- cor(mod02.ord$scores, site_res02[,cor.preds]))
# corrplot(cor(site_res02[,cor.preds]), 
#          method = "ellipse",
#          type = "lower",
#           mar=c(0,0,4,0))
# 
# corrplot(t(mod02.cor), is.corr = F,
#          method = "ellipse", cl.pos = "n")
# 
# dAIC <- AIC(mod01) - AIC(mod02) # will be generally positive  across species if mod2 is better... 
# boxplot(dAIC); abline(h =0, col = "red")




# png("Hmsc_CD/local/corplot_pa.png", width = 300, height = 200, units = "mm", res = 100)
# par(mfrow =c(1,6))
# t_corrplot(mod0.cor, "mod0")
# t_corrplot(mod01.cor, "mod01")
# t_corrplot(mod02.cor, "mod02")
# t_corrplot(mod03.cor, "mod03")
# t_corrplot(mod04.cor, "mod04")
# t_corrplot(mod05.cor, "mod05")
# dev.off()
# 
# dAIC <- list(mod0_mod01 = AIC(mod0)-AIC(mod01),
#              mod0_mod02 = AIC(mod0)-AIC(mod02),
#              mod0_mod03 = AIC(mod0)-AIC(mod03),
#              mod0_mod04 = AIC(mod0)-AIC(mod04),
#              mod0_mod05 = AIC(mod0)-AIC(mod05))
# 
# par(mfrow=c(1,1))
# png("Hmsc_CD/local/AIC_boxplot.png")
# boxplot(dAIC); abline(h = 0, col = "red", lty = 2, las = 2)
# dev.off()
# 
# png("Hmsc_CD/local/cor_boxplot.png")
# boxplot(list(cor0 = abs(mod0.cor),
#              cor1 = abs(mod01.cor),
#              cor2 = abs(mod02.cor),
#              cor3 = abs(mod03.cor),
#              cor4 = abs(mod04.cor),
#              cor5 = abs(mod05.cor)))
# dev.off()



