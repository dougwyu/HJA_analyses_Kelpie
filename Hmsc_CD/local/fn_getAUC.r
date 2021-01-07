

## For a given MF or MFCV list

# library(Hmsc)
# library(dplyr)
# wd <- here::here()
# wd
# setwd(file.path(wd, "Hmsc_CD"))
# 
# ## load results
# # get results folders 
# 
# resF <- list.files("oregon_ada/results", pattern = "res\\d*_\\d{2}$", include.dirs = TRUE, full.names = T)
# resF
# rf <- resF[4]
# rf


getAUC <- function(rf, rMod = TRUE){

  library(Hmsc)
  library(dplyr)
  
  # get model paths
  ms <- list.files(rf, pattern = "^models_.*\\.[Rr]data$", full.names = TRUE, recursive = TRUE)
  
  # get model evaluation paths
  mf <- list.files(rf, pattern = "^MF_.*\\.[Rr]data$", full.names = TRUE, recursive = TRUE)
  
  # check lists same length
  if(!length(ms) == length(mf)) stop("Check lists")
  
  # get highest thin and samples
  thin_vec <- sub(".*_thin_(\\d*)_.*", "\\1", ms)
  samples_vec <- sub(".*_samples_(\\d*)_.*", "\\1", ms)
  thin <- thin_vec[which.max(as.numeric(thin_vec))] # generally higher thin has highest samples.. for now.
  samples <- samples_vec[which.max(as.numeric(thin_vec))]
  rm(thin_vec, samples_vec)
  
  # load model 
  load(ms[grepl(paste0(".*_thin_", thin, "_samples_", samples, ".*\\.Rdata$"), ms)])
  
  # check for errors
  tryE <- !sapply(models, inherits, what = "try-error")
  
  # Get model summary information
  # print models, predictors, no of Obs, random levels, etc
  #m <- models[[1]]; str(m)
  mod.specs <- lapply(models[tryE], function(m) {
    
    if(is.null(m$ranLevelsUsed)) rL <- NA else rL <- m$ranLevelsUsed
    if(is.null(m$XRRRFormula)) RRR_predictors = NA else RRR_predictors = as.character(m$XRRRFormula)[2]
    
    data.frame(
      predictors = as.character(m$XFormula)[2],
      obs = m$ny,
      nSp = m$ns,
      phylo = !is.null(m$phyloTree),
      rL = rL,
      RRR_predictors = RRR_predictors,
      ncRRR = m$ncRRR
    )
    
  })
  
  mod.df <- do.call(rbind, mod.specs)
  mod.df$name <- names(models)[tryE]
  rownames(mod.df) <- NULL  
  mod.df
  rm(mod.specs)
  
  # add ID of resfolder plus number
  mod.df$modID <- paste(basename(rf), sprintf("m%02d", 1:nrow(mod.df)), sep = "_")
  
  head(mod.df)
  
  
  ## Do convergence 
  load(list.files(rf, "beta.*\\.[rR]data", full.names = TRUE, recursive = TRUE))
  
  # Scale reduction factor - beta parameter, mean +/- sd
  all.psrf <- unlist(lapply(beta, function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y)))))
  mod.df$psrf <- all.psrf[match(paste(mod.df$name, thin, samples, sep = "_"), names(all.psrf))]
  rm(beta, all.psrf)
  
  #  get species names
  sp.names <- lapply(seq_along(models), function(i) data.frame(species = colnames(models[[i]]$YScaled),
                                                               modID = mod.df$modID[i]))
  
  sp.names <- do.call(rbind, sp.names)
  
  sp.names <- sp.names %>%
    tidyr::separate(species, into = c("OTU", "empty", "class", "order", "family",
                             "genus", "epithet", "BOLD", "BOLDID",
                             "size"), remove = FALSE, sep = "_") %>%
    select(-empty)%>%
    as.data.frame()
  
  head(sp.names)
  
  # get model evaluation
  load(mf[grepl(paste0(".*_thin_", thin, "_samples_", samples, ".*\\.Rdata$"), mf)]) # WAIC, MF, MFCV
  
  ## get summary of evaluation metrics and add to mod.df
  mfres.list <- lapply(MF[!is.na(MF)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))
  MFres <- t(do.call(data.frame,mfres.list))
  
  mfcvres.list <- lapply(MFCV[!is.na(MFCV)], function(x) sapply(x, function(y) sprintf("%.3f \U00B1 %.3f",mean(y), sd(y))))
  MFCVres <- t(do.call(data.frame, mfcvres.list))
  
  colnames(MFres) <- paste0(colnames(MFres), "_expl")
  colnames(MFCVres) <- paste0(colnames(MFCVres), "_pred")
  
  mfs <- base::merge.default(MFres,MFCVres, by = "row.names")
  mod.df <- cbind(mod.df,  mfs[match(paste(mod.df$name, thin, samples, sep = "_"), mfs$Row.names),2:7])
  
  mod.df$thin <- as.numeric(thin)
  mod.df$samples <- as.numeric(samples)
  
  
  ## Collect all AUC data with species prevalence
  prevList <- lapply(models[tryE], function(m) {
    mev <- data.frame(prev = colSums(m$Y)/nrow(m$Y))
  })
  
  ind <- is.na(MF) # check for and remove NA in MF list (from try errors)
  
  ME_stats <- mapply(function(x,y,z,k,s) {
    
    mfdf <- do.call(data.frame, x)
    colnames(mfdf) <- paste0(colnames(mfdf), "_expl")
    
    mfcvdf <- do.call(data.frame, y)
    colnames(mfcvdf) <- paste0(colnames(mfcvdf), "_pred")
    
    n <- nrow(mfdf)
    
    cbind(rf = rep(basename(rf), n),
          name = rep(k,n),
          s,
          mfdf,
          mfcvdf, 
          prev = z)},
    MF[!ind], MFCV[!ind], prevList[!ind], names(prevList[!ind]), list(sp.names), SIMPLIFY = FALSE)
  
  MF_all <- do.call(rbind,ME_stats)
  rownames(MF_all) <- NULL
  head(MF_all)
  
  # ## Get all evaluation metrics per species and compile in dagta frame
  # tmp1 <- lapply(MFCV[!is.na(MFCV)], function(x) t(do.call(rbind, x)))
  # 
  # # add model name
  # tmp2 <- lapply(seq_along(MFCV), function(i) {
  #   
  #   data.frame(modID2 = mod.df$modID[i], data.frame(tmp1[[i]]))
  # })
  # 
  # # combine all
  # tmp3 <- do.call(rbind, tmp2)
  # head(tmp3)
  # MF.all <- cbind(sp.names, tmp3); rm(tmp1, tmp2, tmp3)
  # 
  # head(MF.all)
  # 
  # # checck names
  # if(!nrow(MF.all) == sum(MF.all$modID == MF.all$modID2)) stop("mod IDs don't match") # should be true
  # MF.all$modID2 <- NULL
  # 
  if(rMod) return(mod.df) else return(MF.all)

}


