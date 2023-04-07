###--------------------------------------------------------------------------###

# generate a spline with 1 knot at the median for one variable

simSpline <- function(x, deg=3, knots=1, coefSpl=NULL){
  
  if (is.null(coefSpl)) coefSpl <- runif(deg + knots + 1, 0, 1)
  bas <- splines::bs(x = x, knots = median(x), degree = deg, intercept = TRUE)
  return(as.vector(bas %*% coefSpl))

}

###--------------------------------------------------------------------------###

# generate a binary outcome with probit gam (splines for cont and nonlinear
# transformations of discrete variables) for a fixed R^2 in the latent variable

simProbit <- function(params, X, Y){

  z <- qnorm(1-params$pm)
  out <- data.table::data.table(LP=rep(NA, params$N))
  is.mnar <- params$mechanism == "mnar"
  X.temp <- cbind(data.table::copy(X), data.table::copy(Y[,..is.mnar]))

  # apply splines to continuous X (and y if is.mnar)
  contNms <- c(names(X.temp)[grep("^RelCont", names(X.temp))], names(Y)[is.mnar])
  X.temp[, (contNms) := lapply(.SD, simSpline, deg = params$deg, coefSpl = unlist(params$theta)),
        .SDcols = contNms]
  
  # apply some transformations to any discrete X
  catNms <- names(X.temp)[grep("(^RelBin)|(^RelOrd)", names(X.temp))]
  if (is.null(params$trans)){
    transf <- lapply(1:length(catNms), function(x) identity)
  } else transf <- unlist(params$trans)
  X.temp[, (catNms) := mapply({function(f, x) f(x)}, transf, .SD, SIMPLIFY = FALSE), 
        .SDcols = catNms]
  
  # scale variables and assign moderately large coefs, the largest if mnar
  allNms <- c(contNms, catNms)
  coefsRel <- rep(1, length(allNms))
  X.temp[, (allNms) := lapply(.SD, scale), .SDcols = allNms]
  coefsRel <- runif(length(allNms), 0.3, 0.4)
  if (is.mnar) coefsRel[which(allNms == names(Y))] <- max(coefsRel)*1.5
  
  # define residual variance for the linear predictor
  explSig <- t(coefsRel) %*% cov(X.temp[, ..allNms]) %*% coefsRel
  if (params$R2r %in% c(0,1)) stop('R**2 for indicator must not be 0 or 1.')
  resSig <- (explSig / params$R2r) - explSig
  
  # define theoretical sd of lp
  sd.LP <- as.numeric(sqrt(explSig + resSig))
    
  # discretise at z*theoretical SD above mean of linear predictor
  out[, LP := as.matrix(X.temp[, ..allNms]) %*% (coefsRel) + 
        stats::rnorm(params$N, 0, sqrt(resSig))]
  return(out[, as.numeric(LP > (mean(LP) + (z) * (sd.LP)))])
  
}

###--------------------------------------------------------------------------###

# simulate missingness with a probit model or randomly

simR <- function(params, X, Y){
  
  if (params$mechanism == "mcar"){
    preds <- NULL
    out <- rbinom(params$N, size = 1, prob = params$pm)
  } else {
    preds <- names(X)[grepl("^Rel", names(X))]
    out <- simProbit(params, X[,..preds], Y)
  }
  return(list(R=out, predsR=preds))
}

###--------------------------------------------------------------------------###
