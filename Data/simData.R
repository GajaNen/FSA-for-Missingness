###--------------------------------------------------------------------------###

# generate ordinal variables from quantiles of logistic distribution

genOrd <- function(logisZ, lP, N){
  
  cprop <- t(apply(lP, 1, cumsum))
  if (any(apply(cprop, 1, max) != rep(1, nrow(cprop)))){
    stop(paste0("Baseline probabilities of ordinal variables sum 
                  up to:", apply(cprop, 1, max), "instead to 1 for each."))
  }
  quant <- stats::qlogis(cprop)
  ordDat <- data.table::data.table(NULL)
  maxs <- max.col(cprop, ties.method = "first")

  for (i in 1:nrow(quant)){
    LogisZi <- logisZ[, i]
    matlp <- matrix(rep(quant[i,], N),
                    ncol = ncol(quant),
                    byrow = TRUE)
    locateGrp <- (LogisZi > cbind(-Inf, matlp))
    cats <- apply(locateGrp, 1, sum)
    if (maxs[i] < 5){
      ordDat[, paste0("Ord", i, "Cat", 2:maxs[i]) := 
               model.matrix(~factor(matrix(cats)))[,-1]]
    } else ordDat[, paste0("Ord",i) := cats]

  }
  
  return(ordDat)
}

###--------------------------------------------------------------------------###

# generate relevant variables of mixed types with a given correlation structure
# (only non-parametric correlation structure retained)

# generete mvrnorm
# norm
# apply each q appropriate number of times
# include Y

getSymMat <- function(rhos, dims){
  
  m <- matrix(NA,nrow=dims,ncol=dims)
  m[lower.tri(m)] <- rhos
  diag(m) <- 1
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

genCorMix <- function(params, Nvar, prfx){
  
  #min(eigen(S, symmetric = TRUE)$values) <= 0
  if (params$corrPred && (prfx=="Rel")) {
    S <- getSymMat(unlist(params$cormat), Nvar)
  } else S <- diag(1, nrow=Nvar)
  paramMarg <- c(unlist(params[[paste0("normParam", prfx)]], recursive = F),
                 unlist(params[[paste0("gamParam", prfx)]], recursive = F),
                 unlist(params[[paste0("binParam", prfx)]], recursive = F),
                 lapply(1:(Nvar / 3),
                        function(x) c(location = 0, scale = 1)))
  funcs <- c(rep(c(qnorm), Nvar / 6),
             rep(c(qgamma), Nvar / 6),
             rep(c(qbinom), Nvar / 3),
             rep(c(qlogis), Nvar / 3))
  dts <- data.table::as.data.table(pnorm(mvnfast::rmvn(1000,rep(0, Nvar),S)))
  
  #freqs <- c(rep(Nvar/6,2),rep(Nvar/3,2))

  dts[, names(dts) :=  mapply(function(f, x) do.call(f, c(list(unname(.SD)),x)), 
                              funcs, paramMarg),
      .SDcols = names(dts)]
  # and Y

  
  ords <- genOrd(logisZ = relt[, (Nvar - Nvar / 3 + 1): Nvar], 
                 lP = params[[paste0("popProbs", prfx)]][[1]],
                 N = params$N)
  relt <- data.table::as.data.table(relt[, 1:(Nvar - Nvar / 3)])
  relt[, names(ords) := ords]
  nm <- c(unlist(lapply(c("Cont", "Bin"), paste0, 1:(Nvar / 3))),
          names(ords)[grepl("^Ord", names(ords))])
  data.table::setnames(relt, paste0(prfx, nm))
  return(relt)
}

###--------------------------------------------------------------------------###

# simulate correlated or independent relevant and irrelevant predictors

simX <- function(params){
  
  Nrel <- params$pr * params$Ntotal
  Nirrl <- params$Ntotal - Nrel
  if ((Nrel %% 6) || (Nirrl %% 6)) stop("Number of the relevant and
                                        irrelevant features each must
                                        be a multiple of 6. Adjust Ntotal
                                        or pr.")
  allv <- genCorMix(params = params, Nvar = Nrel, prfx = "Rel")
  irrelt <- genCorMix(params = params, Nvar = Nirrl, prfx = "Irrel")
  allv[,names(irrelt) := irrelt]
  return(allv)
  
}

###--------------------------------------------------------------------------###

# simulate Y for a given R^2 Y~X

simY <- function(params, X){
  
  rey <- cor(X, X)
  cors <- runif(params$pr * params$Ntotal, 0.3, 0.6)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1] # theoretical bounds for cor?
  cors <- a*cors
  # for given correlations simulate Y
  Xsc <- scale(X)  
  x <- 1:params$N
  e <- residuals(lm(x ~ Xsc)) 
  X.dual <- with(svd(Xsc), (params$N-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - cors %*% cov(X.dual) %*% cors) / var(e))
  y <- data.table::as.data.table(X.dual %*% cors + sqrt(sigma2)*e)
  setnames(y, "Y")
  return(y)
}

###--------------------------------------------------------------------------###
