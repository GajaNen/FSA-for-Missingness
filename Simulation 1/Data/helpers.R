###--------------------------------------------------------------------------###

lambdaSE <- function(x) 
{
  # lambdas ordered from largest to smallest
  # lambda: regularization parameter (larger means larger penalty)
  index <- 1:length(x$cve)
  bestIndex <- which.min(x$cve) # take the largest lambda (first min) with the lowest ME
  perf <- x$cve[bestIndex] + x$cvse[bestIndex]
  candidates <- index[x$cve <= perf]
  best <- min(candidates) # take the smallest index (largest lambda)
  bestMetr <- x$cve[best]
  return(list(bestIndex=best, bestME=bestMetr))
}

###--------------------------------------------------------------------------###
