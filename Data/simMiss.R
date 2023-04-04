###--------------------------------------------------------------------------###

# generate a cubic spline with 1 knot at the median

simSpline <- function(dat, knots=1, deg=4, beta=NULL){
  
  if (is.null(dim(dat))) stop("Data has to be 2D.")
  spl <- apply(dat, 2, function(x) {
    if (is.null(beta)) beta <- sample(0:4, deg + knots + 1, replace=T)
    basis <- splines::bs(x = x, knots = median(x), degree = deg, intercept = TRUE)
    basis %*% beta
    }
  )
  
  return(spl)
   # this could also be performed on columns of datatable! in the main func
}

###--------------------------------------------------------------------------###

simProbit <- function(params, Xrel, Y){

  z <- qnorm(1-params$pm)
  out <- data.table(LP=rep(NA, params$N))
  is.mnar <- params$mechanism == "mnar"

  # apply splines to cont vars and y if mnar
  # to do: y is data.table, make it empty if condition not satisfied
  contNms <- names(Xrel)[grepl("^RelCont", names(Xrel))]
  splns <- simSpline(cbind(
    Xrel[,..contNms],
    Y[is.mnar]),
    deg = params$deg)
  
  # apply some transformations to any discrete X (also ordinal with ncat > 5)
  catNms <- names(Xrel)[grepl("(^RelBin)|(^RelOrd)", names(Xrel))]
  if (is.null(params[["trans"]])){
    transf <- lapply(1:length(catNms), function(x) identity) 
  } else transf <- params[["trans"]]
  transformed <- data.table::copy(Xrel[,..catNms]) # not copied by reference!
  dim(transformed[, catNms := Map({function(f, x) f(x)}, transf, .SD), .SDcols = catNms])
  #dim(transformed[, ..catNms])
  
  preds.all <- cbind(splns, transformed)
  # Define a slope vector (if y.pred make it larger)
  beta <- runif(ncol(preds.all), 1, 2)
  if (predsY) { 
    beta <- c(beta, runif(ncol(Y.preds), max(beta)*10, max(beta)*20))
    preds.all <- cbind(preds.all, Y.preds)
  }
  
  # LP
  out[, LP := as.matrix(preds.all) %*% beta]
  
  # define residual variance for the linear predictor
  sig <- t(beta) %*% cov(preds.all) %*% beta
  if (params$R2r %in% c(0,1)) stop('R**2 for indicator must not be 0 or 1.')
  resVar <- (sig / params$R2r) - sig
  
  # define theoretical sd of lp
  sd_LP <- as.numeric(sqrt(sig + resVar))
    
  # set mean as threshold and discretize at z*theoretical SD above it
  thresh <- mean(out)
  LPs <- out + MASS::mvrnorm(params$N, 0, resVar)
  return(LPs > (thresh + z*sd_LP))
  
}

###--------------------------------------------------------------------------###

# simulate missingness with a probit model or randomly

simR <- function(params, X, Y, predsX){
  
  if (params$mechanism == "mcar"){
    preds <- NULL
    out <- rbinom(params$N, size = 1, prob = params$pm)
  } else {
    preds <- names(X)[grepl("^Rel", names(X))]
    out <- simProbit(params, X[,preds], Y)
  }
  return(list(R=out, predsR=preds))
}

###--------------------------------------------------------------------------###
