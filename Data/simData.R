###--------------------------------------------------------------------------###

# generate ordinal variables from quantiles of logistic distribution

genOrd <- function(logisZ, lP, N){
  
  cprop <- cumsum(lP)
  if (max(cprop) != 1){
    stop("Baseline probabilities of all ordinal variables must sum to 1.")
  }
  quant <- stats::qlogis(cprop)
  matlp <- matrix(rep(quant, N),
                  ncol = length(quant),
                  byrow = TRUE)
  grp <- (logisZ > cbind(-Inf, matlp))
  cats <- apply(grp, 1, sum)
  return(cats)
}

###--------------------------------------------------------------------------###

# generate relevant variables of mixed types with a given correlation structure
# (only non-parametric correlation structure retained)

simCorMix <- function(params, Nvar, prfx, addY=NULL){
  
  namY <- c("Y")[addY]
  paramMarg <- c(unlist(params[[paste0("normParam", prfx)]], recursive = F),
                 unlist(params[[paste0("gamParam", prfx)]], recursive = F),
                 unlist(params[[paste0("binParam", prfx)]], recursive = F),
                 lapply(1:(Nvar / 3),
                        function(x) c(location = 0, scale = 1)),
                 ifelse(addY==T, list(params$paramY),0))
  if (prfx=="Rel") {
    rhos <- unlist(params$corMatRel)
  } else rhos <- unlist(params$corMatIrrel)
  if (!is.null(addY)) Nsim <- Nvar + 1 else Nsim <- Nvar
  copula <- copula::normalCopula(param = rhos, dim = Nsim, dispstr = "un")
  joint <- copula::mvdc(copula = copula, 
                        margins = c(rep("norm", Nvar / 6),
                                    rep("gamma", Nvar / 6),
                                    rep("binom", Nvar /3),
                                    rep("logis", Nvar / 3),
                                    params$margY[addY]), 
                        paramMargins = lapply(paramMarg, function(x) as.list(x)))
  relt <- data.table::as.data.table(copula::rMvdc(1000, joint))
  nms <- unlist(lapply(list("Cont", "Bin", "Ord"), 
                       function(x) paste0(prfx, x, 1:(Nvar / 3))))
  setnames(relt, c(nms,namY))
  ords <- names(relt)[grep("^.*Ord", names(relt))]
  probs <- unlist(apply(params[[paste0("popProbs", prfx)]][[1]],1,list),recursive = F)
  relt[, (ords) := mapply(genOrd, .SD, probs, MoreArgs = list(N=params$N), SIMPLIFY = F), 
       .SDcols = ords]
  return(relt)
}

###--------------------------------------------------------------------------###

# simulate correlated or independent relevant and irrelevant predictors

simDat <- function(params){
  
  Nrel <- params$pr * params$Ntotal
  Nirrl <- params$Ntotal - Nrel
  if ((Nrel %% 6) || (Nirrl %% 6)) stop("Number of the relevant and
                                        irrelevant features each must
                                        be a multiple of 6. Adjust Ntotal
                                        or pr.")
  allv <- simCorMix(params = params, Nvar = Nrel, prfx = "Rel", addY=TRUE)
  irrelt <- simCorMix(params = params, Nvar = Nirrl, prfx = "Irrel")
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
