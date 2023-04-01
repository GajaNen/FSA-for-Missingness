###--------------------------------------------------------------------------###

# simulate correlated or independent relevant and irrelevant predictors

simX <- function(params){
  
  Nrel <- params$pr * params$Ntotal
  if ((Nrel %% 6) || (Nrest %% 6)) stop("Number of relevant and
                                        the irrelevant features each must
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
    baseprobs <- genProbs(params$popProbsRel[[1]])
    cprop <- t(apply(baseprobs, 1, cumsum)) 
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
  colnames(relt) <- paste0("rel", 1:ncol(relt))
  # generate irrelevant
  irrelt <- genIndep(nP = unlist(params$normParamIrrel, recursive = F),
                     gP = unlist(params$gamParamIrrel, recursive = F),
                     bP = unlist(params$binParamIrrel, recursive = F),
                     lP = params$popProbsIrrel[[1]],
                     Nobs = params$N,
                     Nvar = params$Ntotal - Nrel)
  colnames(irrelt) <- paste0("irrel", 1:ncol(irrelt))

  return(data.frame(cbind(relt, irrelt)))
  
}

###--------------------------------------------------------------------------###

# simulate Y for a given R^2 Y~X
# simulate Y for correlations expected given the population correlation

simY <- function(params, X){
  
  rey <- cor(X[,1:(params$percRel * params$Ntotal)], 
             X[,1:(params$percRel * params$Ntotal)])
  # make this the population matrix -- this way it will always be SPD
  # and the correlations will be theoretical, e.g. the R2 will be theoretical
  # and so the actual R2 witht thi data will change with the deviations from
  # the correlations
  cors <- sample(1:10, (params$percRel * params$Ntotal), replace = T)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1]
  cors <- a*cors
  # for given correlations simulate Y
  X <- scale(X[,1:(params$percRel * params$Ntotal)]) 
  x <- 1:params$N
  e <- residuals(lm(x ~ X)) 
  X.dual <- with(svd(X), (params$N-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - cors %*% cov(X.dual) %*% cors) / var(e))
  return(X.dual %*% cors + sqrt(sigma2)*e)
}

###--------------------------------------------------------------------------###