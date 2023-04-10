###--------------------------------------------------------------------------###

# modified function from caret package for choosing subset size in RFE

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
