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
      cats <- model.matrix(~factor(matrix(cats)))[,-1] # make dummies 
      coln <- paste0("Ord", i, "Cat", 2:maxs[i]) # appropriate naming
    } else coln <- paste0("Ord",i)
    ordDat[, (coln) := cats]
  }
  
  return(ordDat)
}

###--------------------------------------------------------------------------###

varNamSet <- function(nams, nums, nmdDT, pattr=".*Ord"){
  
  return(c(unlist(lapply(nams, paste0, nums)),
    names(nmdDT)[grepl(pattr, names(nmdDT))]))
  
}

###--------------------------------------------------------------------------###

# generate relevant variables of mixed types with a given correlation structure
# (only non-parametric correlation structure retained)

genCorMix <- function(params, prfx="Rel"){
  
  Nrel<-params$Ntotal*params$pr
  paramMarg <- c(unlist(params$normParamRel, recursive = F),
                 unlist(params$gamParamRel, recursive = F),
                 unlist(params$binParamRel, recursive = F),
                 lapply(1:(Nrel / 3),
                        function(x) c(location = 0, scale = 1)))
  rhos <- if (!is.na(params$cormat)) unlist(params$cormat) else rep(0, Nrel*(Nrel-1)/2)
  copula <- copula::normalCopula(param = rhos, dim = Nrel, dispstr = "un")
  joint <- copula::mvdc(copula = copula, 
                        margins = unlist(params$marginals), 
                        paramMargins = lapply(paramMarg, function(x) as.list(x)))
  relt <- copula::rMvdc(1000, joint)
  ords <- genOrd(logisZ = relt[, (Nrel - Nrel / 3 + 1): Nrel], 
                 lP = params$popProbsRel[[1]],
                 N = params$N)
  relt <- data.table::as.data.table(relt[, 1:(Nrel - Nrel / 3)])
  # is there a way to transform these two lines to changing column names at the 
  # same time as assingning by reference
  #set(relt, j= names(ords), value=ords)
  relt[, names(ords) := ords]
  nm <- varNamSet(c("Cont", "Bin"), 
                  1:(Nrel / 3), 
                  relt)
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
  allv <- genIndMix(nP = unlist(params$normParamRel, recursive = F),
                     gP = unlist(params$gamParamRel, recursive = F),
                     bP = unlist(params$binParamRel, recursive = F),
                     lP = params$popProbsRel[[1]],
                     Nobs = params$N,
                     Nvar = Nrel,
                     prfx = "Rel")
  irrelt <- genIndMix(nP = unlist(params$normParamIrrel, recursive = F),
                     gP = unlist(params$gamParamIrrel, recursive = F),
                     bP = unlist(params$binParamIrrel, recursive = F),
                     lP = params$popProbsIrrel[[1]],
                     Nobs = params$N,
                     Nvar = Nirrl,
                     prfx = "Irrel")
  
  allv[,names(irrelt) := irrelt]
  return(allv)
  
}

###--------------------------------------------------------------------------###

# simulate Y for a given R^2 Y~X

simY <- function(params, X){
  
  rey <- cor(X)
  cors <- runif(params$pr * params$Ntotal, 0.3, 0.6)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1] # theoretical bounds for cor?
  cors <- a*cors
  # for given correlations simulate Y
  X <- scale(Xrel)  # is this by reference
  x <- 1:params$N
  e <- residuals(lm(x ~ X)) 
  X.dual <- with(svd(X), (params$N-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - cors %*% cov(X.dual) %*% cors) / var(e))
  return(data.table::as.data.table(X.dual %*% cors + sqrt(sigma2)*e, value.name="Y"))
}

###--------------------------------------------------------------------------###
