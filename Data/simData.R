###--------------------------------------------------------------------------###

# so then in simR I can use ^relcont for splines and ^relord ^relbin for trans
# also for ordinal and binary a transformation could be applied, e.g. sin, exp
# bin bin
# cat cat

genOrd <- function(logisZ, lP, N){
  
  cprop <- t(apply(lP, 1, cumsum))
  if (any(apply(cprop, 1, max) != rep(1, nrow(lP)))){
    stop(paste0("Baseline probabilities of ordinal variables sum 
                  up to:", apply(cprop, 1, max), "instead to 1 for each."))
  }
  quant <- stats::qlogis(cprop)
  maxs <- max.col(quant, ties.method = "first")
  ncols_pv <- (maxs >= 5) + (maxs - 1) * (maxs < 5)
  ordDat <- data.table::data.table(matrix(NA, nrow = N, ncol = sum(ncols_pv)))
  coln <- vector(mode = "list", length = nrow(quant))
  strIx <- 1
  
  for (i in 1:nrow(quant)){
    LogisZi <- logisZ[, i]
    matlp <- matrix(rep(quant[i,], N),
                    ncol = ncol(quant),
                    byrow = TRUE)
    locateGrp <- (LogisZi > cbind(-Inf, matlp))
    cats <- apply(locateGrp, 1, sum)
    if (ncols_pv[i] > 1){
      finalOut <- model.matrix(~factor(matrix(cats)))[,-1]
      coln[[i]] <- paste0("Ord",i,"Cat",1:ncols_pv[i]+1)
    } else {
      finalOut <- cats
      coln[[i]] <- paste0("Ord",i)
    }
    endIx <- strIx + ncols_pv[i] - 1
    ordDat[,paste0("V",strIx:endIx) := finalOut]
    strIx <- endIx + 1
  }
  
  setnames(ordDat,unlist(coln))
  return(ordDat)
}


#unlist(lapply(1:nrow(quant), function(x) 
#paste0("Ord",x,"Cat",setdiff(1:ncols_pv[x]+(ncols_pv[x] > 1), 1))))
# but then "Cat remains
#logisZ <- relt[, (Nrel - Nrel / 3 + 1): Nrel, drop = F]
#N <- params$N
#nvar <- (Nrel / 3)

###--------------------------------------------------------------------------###

genIndMix <- function(nP, gP, bP, lP, Nvar, Nobs, prfx=""){
  
  mixs <- data.table::as.data.table(
                cbind(MASS::mvrnorm(n = Nobs,
                              mu = sapply(nP,"[[",1),
                              Sigma = diag(sapply(nP,"[[",2)**2, nrow = (Nvar / 6))),
                sapply(gP, 
                       function(x) stats::rgamma(Nobs, shape = x[1], rate = x[2])),
                sapply(bP, 
                       function(x) stats::rbinom(Nobs, size = x[1], prob = x[2])),
                genOrd(replicate(Nvar / 3, stats::qlogis(stats::runif(Nobs, 0, 1),
                                                         location = 0,
                                                         scale = 1)),
                       lP = lP,
                       N = Nobs)
                )
    )
  nm <- varNamSet(nams=c("Cont", "Bin"), 
                  nums=1:(Nvar / 3), 
                  nmdDT=mixs)
  data.table::setnames(mixs, paste0(prfx, nm))
  return(mixs)
  
}

###--------------------------------------------------------------------------###

varNamSet <- function(nams, nums, nmdDT, pattr=".*Ord"){
  
  return(c(unlist(lapply(nams, paste0, nums)),
    names(nmdDT)[grepl(pattr, names(nmdDT))]))
  
}

###--------------------------------------------------------------------------###

# generate correlated relevant variables of mixed types

genCorMix <- function(params, prfx="Rel"){
  
  Nrel<-params$Ntotal*params$pr
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
  relt <- copula::rMvdc(1000, joint)
  ords <- genOrd(logisZ = relt[, (Nrel - Nrel / 3 + 1): Nrel], 
                 lP = params$popProbsRel[[1]],
                 N = params$N)
  relt <- as.data.table(relt[, 1:(Nrel - Nrel / 3)])
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
  if (params$corrPred){
    relt <- genCorMix(params = params, prfx="Rel")
  } else {
    relt <- genIndMix(nP = unlist(params$normParamRel, recursive = F),
                     gP = unlist(params$gamParamRel, recursive = F),
                     bP = unlist(params$binParamRel, recursive = F),
                     lP = params$popProbsRel[[1]],
                     Nobs = params$N,
                     Nvar = Nrel,
                     prfx = "Rel")
  }
  # generate irrelevant
  irrelt <- genIndMix(nP = unlist(params$normParamIrrel, recursive = F),
                     gP = unlist(params$gamParamIrrel, recursive = F),
                     bP = unlist(params$binParamIrrel, recursive = F),
                     lP = params$popProbsIrrel[[1]],
                     Nobs = params$N,
                     Nvar = Nirrl,
                     prfx = "Irrel")
  
  # join by reference
  return(merge.data.table(relt, irrelt, all = T))
  
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
