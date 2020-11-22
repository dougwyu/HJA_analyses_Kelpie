Manyany<-function (fn, yMat, formula, data, family = "negative.binomial", 
          composition = FALSE, block = NULL, get.what = "details", 
          var.power = NA, na.action = "na.exclude", ...) 
{
  tol = 1e-08
  yMat = as.matrix(yMat)
  data = as.data.frame(data)
  yNames = dimnames(yMat)
  n.rows = dim(yMat)[1]
  n.vars = dim(yMat)[2]
  if (is.null(yNames[[1]])) 
    yNames[[1]] = 1:n.rows
  if (length(yNames) == 1) 
    yNames[[2]] = paste("y", 1:n.vars, sep = "")
  call = match.call()
  if (composition == FALSE) 
    formula = formula(paste("y~", formula[3], sep = ""))
  else {
    yVec = as.vector(yMat)
    ref = factor(rep(1:n.rows, n.vars))
    spp = factor(rep(1:n.vars, each = n.rows))
    data = data.frame(ref, spp, data[ref, ])
    formula = formula(paste("y~ref+spp+spp:(", formula[3], 
                            ")", sep = ""))
    n.rows.orig = n.rows
    n.vars.orig = n.vars
    n.rows = length(yVec)
    n.vars = 1
    names(yVec) = paste(yNames[[2]][spp], ".", yNames[[1]][ref], 
                        sep = "")
    yMat = as.matrix(yVec)
    if (class(family) != "family" & length(family) > 
        1) 
      stop("when using composition=TRUE, family argument must have length one.")
    if (is.null(block)) 
      block = ref
    else block = block[ref]
  }
  if (class(family) == "family" || length(family) == 
      1) {
    fam = family
    family = vector("list", n.vars)
    for (i.var in 1:n.vars) family[[i.var]] = fam
  }
  if (length(family) != n.vars) 
    stop("family argument has length more than one but not equal to the number of columns in yMat (!?)")
  if (length(var.power) == 1) 
    var.power = rep(var.power, n.vars)
  fam = family
  if (fn == "clm") {
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if ("link" %in% names(dots)) 
      link <- dots$link
    else link = "logit"
    if (link == "loglog") 
      fam.i = binomial("cloglog")
    else fam.i = binomial(link)
    fam.i$family = "ordinal"
    fam = vector(mode = "list", length = n.vars)
    family = fam
    yMat = data.frame(yMat,stringsAsFactors=TRUE)
  }
  for (i.var in 1:n.vars) {
    if (is.character(family[[i.var]])) {
      if (family[[i.var]] == "negbinomial" || family[[i.var]] == 
          "negative.binomial") {
        fam[[i.var]] = negative.binomial(10^6)
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=logit)") {
        fam[[i.var]] = binomial()
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=cloglog)" || 
               family[[i.var]] == "cloglog") {
        fam[[i.var]] = binomial("cloglog")
        fam[[i.var]]$family = family[[i.var]]
      }
      else {
        fam.fn = get(fam[[i.var]], mode = "function", 
                     envir = parent.frame())
        fam[[i.var]] = fam.fn()
      }
    }
    if (fn == "clm") 
      fam[[i.var]] = fam.i
    if (fam[[i.var]]$family == "binomial") 
      warning("The binomial option of manyany currently assumes you have binary (presence/absence) response")
  }
  manyfit = vector(mode = "list", length = n.vars)
  fits = matrix(NA, n.rows, n.vars)
  etas = matrix(NA, n.rows, n.vars)
  params = manyfit
  logL = rep(NA, n.vars)
  for (i.var in 1:n.vars) {
    data$y = yMat[, i.var]
    manyfit[[i.var]] = do.call(fn, list(formula = formula, 
                                        family = family[[i.var]], data = data, na.action = na.action, 
                                        ...))
    logL[i.var] = logLik(manyfit[[i.var]])
    if (is.na(logL[i.var])) 
      logL[i.var] = -0.5 * deviance(manyfit[[i.var]])
    if (get.what == "details" || get.what == "models") {
      fits[, i.var] = fitted(manyfit[[i.var]])
      etas[, i.var] = switch(fn, lmer = manyfit[[i.var]]@eta, 
                             clm = predict(manyfit[[i.var]], type = "linear.predictor", 
                                           newdata = data[, names(data) != "y"])$eta1[, 
                                                                                      1], predict(manyfit[[i.var]]))
      if (substr(fam[[i.var]]$family, 1, 3) == "bin" || 
          fam[[i.var]]$family == "ordinal") {
        etas[, i.var] = pmax(etas[, i.var], fam[[i.var]]$linkfun(tol)/2)
        etas[, i.var] = pmin(etas[, i.var], fam[[i.var]]$linkfun(1 - 
                                                                   tol)/2)
      }
      if (fam[[i.var]]$link == "log" || fam[[i.var]]$link == 
          "mu^0") 
        etas[, i.var] = pmax(etas[, i.var], log(tol)/2)
      if (i.var == 1) {
        cf = try(coef(manyfit[[i.var]]), silent = TRUE)
        if (inherits(cf, "try-error")) {
          do.coef = FALSE
          coefs = NULL
        }
        else {
          coefs = vector(mode = "list", n.vars)
          coefs[[1]] = cf
          names(coefs[[1]]) = dimnames(cf)[[1]]
          if (composition == FALSE & n.vars > 1) 
            names(coefs) = yNames[[2]]
          do.coef = TRUE
        }
      }
      else {
        if (do.coef == TRUE) 
          coefs[[i.var]] = coef(manyfit[[i.var]])
      }
      if (fam[[i.var]]$family == "poisson") 
        params[[i.var]] = list(q = yMat[, i.var], lambda = fits[, 
                                                                i.var])
      if (substr(fam[[i.var]]$family, 1, 3) == "bin") 
        params[[i.var]] = list(q = yMat[, i.var], prob = fits[, 
                                                              i.var], size = 1)
      if (fam[[i.var]]$family == "Tweedie") 
        params[[i.var]] = list(q = yMat[, i.var], power = var.power[i.var], 
                               mu = fits[, i.var], phi = summary(manyfit[[i.var]])$disp)
      if (fam[[i.var]]$family == "ordinal") 
        params[[i.var]] = list(q = yMat[, i.var], mu = predict(manyfit[[i.var]], 
                                                               type = "cum.prob"), muAll = predict(manyfit[[i.var]], 
                                                                                                   type = "cum.prob", newdata = data[, names(data) != 
                                                                                                                                       "y"])$cprob2)
      if (grepl("egative", fam[[i.var]]$family) || 
          fam[[i.var]]$family == "negbinomial") {
        if (any(names(manyfit[[i.var]]) == "theta")) 
          theta = manyfit[[i.var]]$theta
        else {
          if (any(names(manyfit[[i.var]]) == "phi")) 
            theta = 1/manyfit[[i.var]]$phi
          else theta = 1/(fam[[i.var]]$var(1) - 1)
        }
        params[[i.var]] = list(q = yMat[, i.var], mu = fits[, 
                                                            i.var], size = theta)
      }
      if (fam[[i.var]]$family == "gaussian") {
        s.ft = summary(manyfit[[i.var]])
        if (any(names(s.ft) == "sigma")) 
          sd = s.ft$sigma
        else sd = s.ft$scale
        params[[i.var]] = list(q = yMat[, i.var], mean = fits[, 
                                                              i.var], sd = sd)
      }
    }
  }
  object = list(logL = logL, get.what = get.what)
  if (get.what == "details" || get.what == "models") {
    if (composition == TRUE) {
      fits = matrix(fits, n.rows.orig, n.vars.orig)
      etas = matrix(etas, n.rows.orig, n.vars.orig)
    }
    else {
      if (n.vars > 1) {
        names(logL) = yNames[[2]]
        names(params) = yNames[[2]]
      }
    }
    attributes(logL)$df = attributes(logLik(manyfit[[i.var]]))$df
    attributes(logL)$nobs = n.rows
    class(logL) = "logLik"
    resids = residuals.manyany(list(params = params, family = fam, 
                                    composition = composition, fitted.values = fits, 
                                    get.what = get.what))
    dimnames(resids) = yNames
    dimnames(fits) = yNames
    dimnames(etas) = yNames
    mf = if (any(names(manyfit[[i.var]]) == "data")) 
      manyfit[[i.var]]$data
    else model.frame(manyfit[[i.var]])
    object = list(logL = logL, fitted.values = fits, residuals = resids, 
                  linear.predictor = etas, family = fam, coefficients = coefs, 
                  call = call, params = params, model = mf, terms = terms(manyfit[[i.var]]), 
                  formula = formula, block = block, composition = composition, 
                  get.what = get.what)
  }
  if (get.what == "models") {
    object$fits = manyfit
    names(object$fits) = yNames[[2]]
  }
  class(object) = c("manyany", class(manyfit[[i.var]]))
  return(object)
}