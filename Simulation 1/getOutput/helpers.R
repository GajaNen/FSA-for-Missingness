###--------------------------------------------------------------------------###

# modified function from caret package for choosing subset size in RFE
# specifically to implement one SE rule for tuning subset size

# set of function specifically for random forest
rfFuncs$selectSize <- function (x, metric, maximize) 
{
  index <- 1:nrow(x)
  if (!maximize) {
    bestIndex <- which.min(x[, metric])
    perf <- x[bestIndex, metric] + (x[bestIndex, paste(metric, 
                                                       "SD", sep = "")])/sqrt(1)
    candidates <- index[x[, metric] <= perf]
  }
  else {
    bestIndex <- which.max(x[, metric])
    perf <- x[bestIndex, metric] - (x[bestIndex, paste(metric, 
                                                       "SD", sep = "")])/sqrt(1)
    candidates <- index[x[, metric] >= perf]
  }
  best <- min(candidates)
  x[best, "Variables"]
  
}

###--------------------------------------------------------------------------###

# set of function for any functions (in my case used for svm)
caretFuncs$selectSize <- function (x, metric, maximize) 
{
  index <- 1:nrow(x)
  if (!maximize) {
    bestIndex <- which.min(x[, metric])
    perf <- x[bestIndex, metric] + (x[bestIndex, paste(metric, 
                                                       "SD", sep = "")])/sqrt(1)
    candidates <- index[x[, metric] <= perf]
  }
  else {
    bestIndex <- which.max(x[, metric])
    perf <- x[bestIndex, metric] - (x[bestIndex, paste(metric, 
                                                       "SD", sep = "")])/sqrt(1)
    candidates <- index[x[, metric] >= perf]
  }
  best <- min(candidates)
  x[best, "Variables"]
}

###--------------------------------------------------------------------------###
# one SE rule for sparseSVM's CV results
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
