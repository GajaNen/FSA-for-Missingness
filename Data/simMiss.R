###--------------------------------------------------------------------------###

simR <- function(params, X, Y, predsX){
  
  coln <- paste0("R", 1:params$nY)
  
  if (params$mechanism == "mcar"){
    dat <- data.frame(replicate(params$nY, rbinom(params$N, 1, params$pm)))
    colnames(dat) <- coln
    # sim all: probit, gam, rf
    return(list(dat=dat,preds=NULL))
  }
  
  dat <- data.frame(ampUni(params, X, Y, predsX))
  colnames(dat) <- coln
  return(list(dat=dat, preds=predsX))

}

###--------------------------------------------------------------------------###

simProbit <- function(params, X, Y, predsX){
  
  mult <- qnorm(1-params$pm)
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
      Y.preds <- creatDumm(Y[,i, drop=F])
    } else Y.preds <- params$trans[[max(xcontpr)+1]](Y[,predsY,drop=F])
    
    # apply functions with non linear transformations to cont X
    transX <- do.call(cbind, lapply(1:xcontpr, function(c) 
      params$trans[[c]](X[,((predsX[[i]])[-xcat.ind]),drop=F][,c,drop=F])))
    
    dat.all <- cbind(X.dummy, transX)
    
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
    resVar[i] <- (sig / params$R2r) - sig
    
    # define theoretical sd of lp
    sd_LP[i] <- as.numeric(sqrt(sig + resVar[i]))
    
  }
  
  # generate R for each Y independently
  thresholds <- colMeans(dats)
  LPs <- dats + mvrnorm(params$N, 
                        rep(0, params$nY), 
                        diag(x=resVar, nrow=params$nY))
  return(LPs > matrix(thresholds + mult*sd_LP, 
                      ncol=params$nY, 
                      nrow=params$N,
                      byrow=T))
  
}

###--------------------------------------------------------------------------###
