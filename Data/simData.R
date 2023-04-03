###--------------------------------------------------------------------------###

# TO DO: add potentiall dummy coding already here 
# and add names of type of var
# e.g. RelCont, RelOrd, RelBin, if dummy RelOrd1Dum1
# so then in simR I can use ^relcont for splines and ^relord ^relbin for trans
# also for ordinal and binary a transformation could be applied, e.g. sin, exp

genOrd <- function(logisZ, quant, N, nvar){
  
  ordDat <- matrix(NA, nrow = N, ncol = nvar)
  
  for (i in 1:nvar){
    iLogis <- logisZ[,i]
    matlp <- matrix(rep(quant[i,], N),
                    ncol = ncol(quant),
                    byrow = TRUE
    )
    
    locateGrp <- (iLogis > cbind(-Inf, matlp))
    ordDat[,i] <- apply(locateGrp, 1, sum)
  }
  
  return(as.data.frame(ordDat))
}

#logisZ <- relt[, (Nrel - Nrel / 3 + 1): Nrel, drop = F]
#N <- params$N
#nvar <- (Nrel / 3)

###--------------------------------------------------------------------------###

genIndep <- function(nP, gP, bP, lP, Nvar, Nobs){
  
  cprop <- t(apply(lP, 1, cumsum))
  if (any(apply(cprop, 1, max) != rep(1, length(bP)))){
    stop(paste0("Baseline probabilities of ordinal variables sum 
                  up to:", apply(cprop, 1, max), "instead to 1 for each var."))
  }
  
  vars <- cbind(MASS::mvrnorm(n = Nobs,
                              mu = sapply(nP,"[[",1),
                              Sigma = diag(sapply(nP,"[[",2)**2, nrow = (Nvar / 6))),
                sapply(gP, 
                       function(x) stats::rgamma(Nobs, shape = x[1], rate = x[2])),
                sapply(bP, 
                       function(x) stats::rbinom(Nobs, size = x[1], prob = x[2])),
                genOrd(replicate(Nvar / 3, stats::qlogis(stats::runif(Nobs, 0, 1),
                                                         location = 0,
                                                         scale = 1)),
                       quant = t(apply(cprop, 1, stats::qlogis)),
                       N = Nobs,
                       nvar = (Nvar / 3)
                )
  )
  
  return(vars)
  
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
  if (params$corrPred){
    paramMarg <- c(unlist(params$normParamRel, recursive = F),
                   unlist(params$gamParamRel, recursive = F),
                   unlist(params$binParamRel, recursive = F),
                   lapply(1:(Nrel / 3),
                          function(x) c(location = 0, scale = 1)))
    copula <- copula::normalCopula(param = unlist(params$cormat), 
                                   dim = Nrel,
                                   dispstr = "un")
    joint <- copula::mvdc(copula = copula, 
                          margins = unlist(params$marginals), 
                          paramMargins = lapply(paramMarg, function(x) as.list(x)))
    relt <- as.data.frame(copula::rMvdc(1000, joint))
    cprop <- t(apply(params$popProbsRel[[1]], 1, cumsum)) 
    if (any(apply(cprop, 1, max) != rep(1, nrow(cprop)))){
      stop(paste0("Baseline probabilities of relevant variables sum 
                  up to:", apply(cprop, 1, max), "instead to 1 for each var."))
    }
    quant <- t(apply(cprop, 1, stats::qlogis))
    relt[, (Nrel - Nrel / 3 + 1): Nrel] <- 
      genOrd(logisZ = relt[, (Nrel - Nrel / 3 + 1): Nrel, drop = F], 
             quant = quant,
             N = params$N,
             nvar = (Nrel / 3))
  } else {
    relt <- genIndep(nP = unlist(params$normParamRel, recursive = F),
                     gP = unlist(params$gamParamRel, recursive = F),
                     bP = unlist(params$binParamRel, recursive = F),
                     lP = params$popProbsRel[[1]],
                     Nobs = params$N,
                     Nvar = Nrel)
  }
  colnames(relt) <- paste0("rel", 1:Nrel)
  # generate irrelevant
  irrelt <- genIndep(nP = unlist(params$normParamIrrel, recursive = F),
                     gP = unlist(params$gamParamIrrel, recursive = F),
                     bP = unlist(params$binParamIrrel, recursive = F),
                     lP = params$popProbsIrrel[[1]],
                     Nobs = params$N,
                     Nvar = Nirrl)
  colnames(irrelt) <- paste0("irrel", (Nrel+1):params$Ntotal)

  return(data.frame(cbind(relt, irrelt)))
  
}

###--------------------------------------------------------------------------###

# simulate Y for a given R^2 Y~X
# simulate Y for correlations expected given the population correlation

simY <- function(params, Xrel){
  
  rey <- cor(Xrel)
  # make this the population matrix -- this way it will always be SPD
  # and the correlations will be theoretical, e.g. the R2 will be theoretical
  # and so the actual R2 witht thi data will change with the deviations from
  # the correlations
  cors <- sample(1:10, (params$pr * params$Ntotal), replace = T)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1]
  # make sure that all cors all reasonably high (e.g. above 0.2, so
  # that it's not that one does not predict - important for mnar)
  cors <- a*cors
  # for given correlations simulate Y
  X <- scale(Xrel) 
  x <- 1:params$N
  e <- residuals(lm(x ~ X)) 
  X.dual <- with(svd(X), (params$N-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - cors %*% cov(X.dual) %*% cors) / var(e))
  return(X.dual %*% cors + sqrt(sigma2)*e)
}

###--------------------------------------------------------------------------###