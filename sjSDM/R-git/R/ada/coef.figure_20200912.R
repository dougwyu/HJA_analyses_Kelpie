# code from Cai Wang, edited by Doug

coef.figure <- function(summary.p, result, minsize, maxsize=200000, taxon="all") {
  effect <- data.frame(summary.p$coefmat)
  effect$rownames=rownames(effect)
  # effect <- rownames_to_column(effect, var = "rownames")
  effect<-separate(data = effect, 
                   col = rownames, 
                   into = c("species", "coef"), 
                   sep = " "
                   )
  
  ifelse(taxon != "all",
     effect <- filter(effect, grepl(taxon, species)),
     effect <- effect
     )
  
  effect<-effect %>% filter(coef != "(Intercept)")
  effect<-effect %>% dplyr::select(-c(Z.value,Pr...z..))
  effect$coef<-as.factor(effect$coef)
  effect$max<-NA
  #str(effect$coef)
  n=ncol(result$sigma) # number of species
  m=ncol(data.frame(result$beta$env))-1
  for (i in 1:n) {
    max_effects = as.numeric(
      apply(data.frame(effect[((m*i-(m-1)):(m*i)),1]),2, 
            function(e) which.max(abs(e))
            )
      ) 
    max_effects=max_effects + ((i-1)*m)
    effect$max[max_effects]= "max"
  }
  effect$max[is.na(effect$max)] <- " "
  effect$group=NA
  effect$group[grep("Diptera",effect$species)]="Diptera"
  effect$group[grep("Hymenoptera",effect$species)]="Hymenoptera"
  effect$group[grep("Araneae",effect$species)]="Araneae"
  effect$group[grep("Psocodea",effect$species)]="Psocodea"
  effect$group[grep("Hemiptera",effect$species)]="Hemiptera"
  effect$group[grep("Coleoptera",effect$species)]="Coleoptera"
  effect$group[grep("Neuroptera",effect$species)]="Neuroptera"
  effect$group[grep("Lepidoptera",effect$species)]="Lepidoptera"
  effect$group[grep("Lepidoptera",effect$species)]="Lepidoptera"
  effect$group[grep("Raphidioptera",effect$species)]="Raphidioptera"
  effect$group[grep("Orthoptera",effect$species)]="Orthoptera"
  effect$group[grep("Opiliones",effect$species)]="Opiliones"  
  effect$group[grep("Plecoptera",effect$species)]="Plecoptera"  
  effect$coef=as.factor(effect$coef)
  #levels species 
  
  effect <- effect %>% 
    tidyr::separate(species, into = c("OTU", "taxon"), sep = "__", remove = TRUE) %>% 
    tidyr::separate(taxon, into = c("taxon", "size"), sep = "_size.", remove = TRUE) %>% 
    mutate(
      taxon = str_replace(taxon, "_sp", "_genus_sp")
    ) %>% 
    tidyr::separate(taxon, into = c("class", "order", "family", "genus", "genus2", "species", "BOLD", "BOLDID")) %>% 
    unite(OTU, c("OTU", "class", "order", "family"), sep = "_", ) %>% 
    select(-genus, -genus2, -species, -BOLD, -BOLDID) %>%
    mutate(
      size = as.numeric(size)
    ) %>% 
    rename(species = OTU) %>% 
    filter(size>minsize & size < maxsize)

  
  effect$species=as.factor(effect$species)
  effect$species=factor(effect$species, levels=rev(levels(effect$species)))
  #effect$species=factor(effect$species,levels = levels(effect$species)[c(26:31,33:40,15:25,41:47,1:10,11:14,32)])
  #effect$species <- relevel(effect$species, "Wafl_Pele_Arde_Ardecine")
  
  #reorder levels coef 
  # if (length((levels(effect$coef)))==13)  effect$coef=factor(effect$coef, levels = c('HSI2_Pond_area','HSI3_Pond_drying','HSI4_Water_quality','HSI5_Shade','HSI8_Pond_count','HSI10_Macrophytes','Inflow_present','Outflow_present','PC1','PC2','PC3','PC4','PC5')) 
  # 
  # #  else effect$coef=factor(effect$coef,levels = c('OS_area','cropN','cropP','cropK'))
  # effect <- effect %>% filter(effect$max == "max")
  t= -(abs(effect$Estimate[which.max(abs(effect$Estimate))]))
  #t=-0.8
  #figure
  col <- brewer.pal(12, "Paired")
  p1 <- ggplot(effect, 
               aes(x = species, 
                   y = Estimate, 
                   fill = group)) +
    geom_bar(position = position_dodge(0.6), 
             stat="identity", width = 0.5, size=0.5) +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(reverse=F)) +
    xlab("species") + 
    ylab("coef") + 
    labs(fill="Species Index") + 
    coord_flip(expand=F) + 
    geom_hline(aes(yintercept = 0), 
               linetype="dotted", 
               size=0.5) +
    theme_classic()+ 
    facet_wrap(~coef, ncol = 5) +
    geom_text(aes(y = t, label = max), 
              position = position_dodge(0.6),
              size = 1, fontface = "bold") + 
    # geom_errorbar(aes(ymax = Estimate + Std.Err, ymin = Estimate - Std.Err), width = 0.3) +
    scale_y_continuous(limits = c(-1.1,1)) +
    theme(axis.text.x = element_text(size = 3)) +
    theme(axis.text.y = element_text(size = 3))
  return(p1)
}


# effecthist <- effect %>% distinct(species, .keep_all = TRUE) 
# hist(effecthist$size, breaks = 5000, xlim = c(0, 3000))
