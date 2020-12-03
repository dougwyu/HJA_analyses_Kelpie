

## Only local: 
# setwd("J:/UEA/Oregon")
# dir()

## On ADA
## getwd() will be "/gpfs/home/hsp20azu"
# with folders Oregon, etc... 
setwd("~/Oregon_winscp")
dir()
rm(list = ls())

library(Hmsc)


fs <- list.files("models", "models_thin.*", full.names = TRUE)
fs

# i = 3
# m <- models[[2]]

pdf("results/biplots_%02d.pdf")

for (i in seq_along(fs)) {
  
  filename = fs[i]
  samples <- sub(".*samples_([[:digit:]]{1,4}).*", "\\1", filename)
  
  load(filename)
  
  m.names <- names(models)[grepl("_sp", names(models))]
  
  for(nm in m.names){
  #nm = m.names[1]
    
    m <- models[[nm]]
    etaPost = getPostEstimate(m, "Eta")
    lambdaPost = getPostEstimate(m, "Lambda")
    
    ### without asp = 1
    biPlot_hmsc(m, etaPost = etaPost, lambdaPost = lambdaPost,
           colVar = "canopy.ht", xlim = xlims, ylim = ylims,
           main = paste0(nm, "_smpl_", samples), spNames = NA)

    biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost,
                colVar = "yrs.disturb.min", xlim = xlims, ylim = ylims,
                main = paste0(nm, "_smpl_", samples), spNames = NA)
    
    
  }
  
  

    
  
}


dev.off()