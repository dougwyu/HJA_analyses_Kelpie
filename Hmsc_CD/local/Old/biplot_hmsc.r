biPlot_hmsc <- function (hM, etaPost, lambdaPost, factors = c(1, 2), colVar = NULL, 
          colors = NULL, spCols = "black", spNames = hM$spNames, ...) 
{
  if (!is.null(colVar)) {
    col = hM$XData[, colVar]
    if (!is.factor(col)) {
      if (is.null(colors)) {
        colors = colorRampPalette(c("blue", "white", 
                                    "red"))
      }
      cols = colors(100)
      plotorder = order(hM$XData[, which(colnames(hM$XData) == 
                                           colVar)])
    }
    else {
      if (is.null(colors)) {
        colors = palette("default")
        cols = colors[col]
      }
      cols = colors[col]
      plotorder = 1:hM$ny
    }
  }
  else {
    cols = "grey"
    plotorder = 1:nrow(hM$XData)
  }
  scale1 = abs(c(min(etaPost$mean[, factors[1]]), max(etaPost$mean[, 
                                                                   factors[1]])))
  scale2 = abs(c(min(etaPost$mean[, factors[2]]), max(etaPost$mean[, 
                                                                   factors[2]])))
  scale1 = min(scale1/abs(c(min(lambdaPost$mean[factors[1], 
  ]), max(lambdaPost$mean[factors[1], ]))))
  scale2 = min(scale2/abs(c(min(lambdaPost$mean[factors[2], 
  ]), max(lambdaPost$mean[factors[2], ]))))
  scale <- min(scale1, scale2)
  
  # plot sites
  plot(etaPost$mean[, factors[1]][plotorder], etaPost$mean[,factors[2]][plotorder], 
       pch = 16, col = cols, xlab = paste("Latent variable",factors[1]), 
       ylab = paste("Latent variable", factors[2]), ...) # removed asp=1 here.. 
  
  #plot species
  points(scale * lambdaPost$mean[factors[1], ], scale * lambdaPost$mean[factors[2],], 
         pch = 17, cex = 1, col = spCols) ## added species colors here
  
  
  text(scale * lambdaPost$mean[factors[1], ], scale * lambdaPost$mean[factors[2],], spNames, pos = 1, cex = 1)
}