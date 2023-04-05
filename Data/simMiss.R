###--------------------------------------------------------------------------###

# generate a spline with 1 knot at the median for one variable

simSpline <- function(x, deg=3, knots=1, coefSpl=NULL){
  
  if (is.null(coefSpl)) coefSpl <- runif(deg + knots + 1, 0, 1)
  bas <- splines::bs(x = x, knots = median(x), degree = deg, intercept = TRUE)
  # plot(x, bas %*% coefSpl)
  return(as.vector(bas %*% coefSpl))

}

###--------------------------------------------------------------------------###

simProbit <- function(params, Xrel, Y){

  z <- qnorm(1-params$pm)
  out <- data.table::data.table(LP=rep(NA, params$N))
  is.mnar <- params$mechanism == "mnar"

  # apply splines to cont vars (and y if is.mnar)
  contNms <- names(Xrel)[grepl("^RelCont", names(Xrel))]
  splns <- cbind(data.table::copy(Xrel[,..contNms]), Y[,..is.mnar]) 
  splns[, (contNms) := lapply(.SD, simSpline, deg = params$deg),
        .SDcols = contNms]
  
  # apply some transformations to any discrete X
  catNms <- names(Xrel)[grepl("(^RelBin)|(^RelOrd)", names(Xrel))] # for safety
  if (is.null(params$trans)){
    transf <- lapply(1:length(catNms), function(x) identity)
  } else transf <- params$trans
  transformed <- data.table::copy(Xrel[,..catNms])
  transformed[, (catNms) := mapply({function(f, x) f(x)}, transf, .SD, SIMPLIFY = FALSE), 
              .SDcols = catNms]
  
  preds.all <- cbind(splns, transformed)
  coefsRel <- rep(1, length(preds.all))
  
  # LP
  out[, LP := as.matrix(preds.all) %*% coefsRel]
  
  # define residual variance for the linear predictor
  s <- t(coefsRel) %*% cov(preds.all) %*% coefsRel
  if (params$R2r %in% c(0,1)) stop('R**2 for indicator must not be 0 or 1.')
  resSig <- (s / params$R2r) - s
  
  # define theoretical sd of lp
  sd.LP <- as.numeric(sqrt(s + resSig))
    
  # set mean as threshold and discretize at z*theoretical SD above it
  out[, LP := LP + MASS::mvrnorm(params$N, 0, resSig)]
  return(out[, as.numeric(LP > (mean(LP)+ (z)*(sd.LP)))])
  
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
