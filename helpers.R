###--------------------------------------------------------------------------###

creatDumm <- function(dat) data.frame(apply(dat, 2,function(x) {model.matrix(~factor(x))[,-1]}))

###--------------------------------------------------------------------------###

extrReg <- function(rstr, symb="\\d+") as.numeric(unique(unlist(regmatches(rstr, regexpr(symb, rstr))))) 

###--------------------------------------------------------------------------###

lambdaSE <- function(x) 
{
  index <- 1:length(x$cve)
  bestIndex <- which.min(x$cve)
  perf <- x$cve[bestIndex] + x$cve[bestIndex]
  candidates <- index[x$cve <= perf]
  best <- min(candidates)
  bestMetr <- x$cve[best]
  return(list(best, bestMetr))
}

###--------------------------------------------------------------------------###

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
  if (apply(cprop, 1, max) != rep(1, nrow(bP))){
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