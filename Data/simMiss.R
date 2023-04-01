###--------------------------------------------------------------------------###

# generate a cubic spline with 1 knot at the median

simSpline <- function(dat, knots=1, deg=4, beta=NULL){
  
  spl <- apply(dat, 2, function(x) {
    if (is.null(beta)) beta <- sample(0:4, deg + knots + 1, replace=T)
    basis <- splines::bs(x = x, knots = median(x), degree = deg, intercept = TRUE)
    basis %*% beta
    }
  )
  
  return(as.data.frame(spl))
  
}

###--------------------------------------------------------------------------###

simProbit <- function(params, X, Y, trans=F, spl=F){
  
  # now only relevant variables are added so no need to subset
  # just apply transformations / splines 
  # no need to loop over nY since only Y
  z <- qnorm(1-params$pm)
  resVar <- sd_LP <- rep(0, nrow=params$nY)
  dats <- matrix(NA, nrow=params$N, ncol=params$nY)
  
  for (i in 1:params$nY){
    
    predsY <- i * (params$mechanism == "mnar")
    ycat <- (params$typeY)[i] == "cat"
    xcontpr <- params$infoX[2]/params$nY
    xcat.ind <- 1:(params$infoX[1]/params$nY)
    
    # create dummy variables for non-cont vars
    X.dummy <- creatDumm(X[,predsX[[i]][xcat.ind]])
    if (predsY && ycat) {
      Y.preds <- creatDumm(Y[,i,drop=F])
    } else if (predsY && trans) { # could be moved to the lower loop and next line is here instead
      Y.preds <- params$trans[[1]](Y[,predsY,drop=F])
    } else Y.preds <- Y[,predsY,drop=F]
    
    # apply functions with non linear transformations to cont X
    if (trans){
      X.cont <- do.call(cbind, lapply(2:(xcontpr+1), function(c) 
        params$trans[[c]](X[,((predsX[[i]])[-xcat.ind]),drop=F][,c,drop=F])))
    } else X.cont <- X[,((predsX[[i]])[-xcat.ind]),drop=F]
    
    if (spl) X.cont <- simSpline(X.cont, deg=params$deg)
    
    dat.all <- cbind(X.dummy, X.cont)
    
    # Define a slope vector
    beta <- runif(ncol(dat.all), 1, 2)
    if (predsY) { 
      beta <- c(beta, runif(ncol(Y.preds), max(beta)*10, max(beta)*20))
      dat.all <- cbind(dat.all, Y.preds)
    }
    
    # LP
    dats[,i] <- as.matrix(dat.all) %*% beta
    
    # define residual variance for the linear predictor
    sig <- t(beta) %*% cov(dat.all) %*% beta
    if (params$R2r == 0) stop('R**2 must be different from 0 for division.')
    resVar[i] <- (sig / params$R2r) - sig
    
    # define theoretical sd of lp
    sd_LP[i] <- as.numeric(sqrt(sig + resVar[i]))
    
  }
  
  # generate R for each Y independently
  thresholds <- colMeans(dats)
  LPs <- dats + mvrnorm(params$N, 
                        rep(0, params$nY), 
                        diag(x=resVar, nrow=params$nY))
  return(LPs > matrix(thresholds + z*sd_LP, 
                      ncol=params$nY, 
                      nrow=params$N,
                      byrow=T))
  
}

###--------------------------------------------------------------------------###

# simulate missingness in a (non-)linear fashion with a probit-like model

simR <- function(params, X, Y, predsX){
  
  if (params$mechanism == "mcar"){
    preds <- NULL
    dat <- rbinom(params$N, size = 1, prob = params$pm)
    colnames(dat) <- c("Rrand")
  } else {
  preds <- grepl("^rel", colnames(X))
  dat <- cbind(data.frame(simProbit(params, X[,preds], Y)), 
               data.frame(simProbit(params, X[,preds], Y, spl=T)))
  colnames(dat) <- c("Rlin", "Rnonlin")
  }
  return(list(dat=dat,preds=preds))
}

###--------------------------------------------------------------------------###