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
# mar or mnar data

simProbit <- function(params, dat){

  z <- qnorm(1-params$pm)
  out <- data.table::data.table(LP=rep(NA, params$N))
  is.mnar <- params$mechanism == "mnar"
  contNms <- c(names(dat)[grep("^RelCont", names(dat))], c("Y")[is.mnar])
  catNms <- names(dat)[grep("(^RelBin)|(^RelOrd)", names(dat))]
  allNms <- c(contNms, catNms)
  dat.rel <- data.table::copy(dat[,..allNms])

  # apply splines to continuous X (and y if is.mnar)
  dat.rel[, (contNms) := lapply(.SD, simSpline, deg = params$deg, coefSpl = unlist(params$theta)),
          .SDcols = contNms]
  
  # apply some transformations to any discrete X
  if (is.null(params$trans)){
    transf <- lapply(1:length(catNms), function(x) identity)
  } else transf <- unlist(params$trans)
  dat.rel[, (catNms) := mapply({function(f, x) f(x)}, transf, .SD, SIMPLIFY = FALSE), 
          .SDcols = catNms]
  
  # scale variables and assign moderately large coefs, the largest if mnar
  dat.rel[, (allNms) := lapply(.SD, scale), .SDcols = allNms]
  coefsRel <- runif(length(allNms), 0.3, 0.4)
  if (is.mnar) coefsRel[which(allNms == "Y")] <- max(coefsRel)*1.5
  
  # define residual variance for the linear predictor
  explSig <- t(coefsRel) %*% cov(dat.rel) %*% coefsRel
  if (params$R2r %in% c(0,1)) stop('R**2 for indicator must not be 0 or 1.')
  resSig <- (explSig / params$R2r) - explSig
  
  # define theoretical sd of lp
  sd.LP <- as.numeric(sqrt(explSig + resSig))
    
  # discretise at z*theoretical SD above mean of linear predictor
  out[, LP := as.matrix(dat.rel) %*% (coefsRel) + 
        stats::rnorm(params$N, 0, sqrt(resSig))]
  return(list(R=out[, as.numeric(LP > (mean(LP) + (z) * (sd.LP)))],
              preds=allNms))
  
}

###--------------------------------------------------------------------------###

# simulate missingness with a probit model or randomly

simR <- function(params, dat){
  
  if (params$mechanism == "mcar"){
    preds <- NULL
    R <- rbinom(params$N, size = 1, prob = params$pm)
  } else {
    out.lst <- simProbit(params, dat)
    R <- out.lst$R
    preds <- out.lst$preds
  }
  return(list(R=R, preds=preds))
}

###--------------------------------------------------------------------------###
