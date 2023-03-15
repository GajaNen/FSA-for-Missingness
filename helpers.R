###--------------------------------------------------------------------------###

creatDumm <- function(dat) data.frame(apply(dat, 2,function(x) {model.matrix(~factor(x))[,-1]}))

###--------------------------------------------------------------------------###

extrReg <- function(rstr, symb="\\d+") as.numeric(unique(unlist(regmatches(rstr, regexpr(symb, rstr))))) 

###--------------------------------------------------------------------------###

lambdaSE <- function (x) 
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