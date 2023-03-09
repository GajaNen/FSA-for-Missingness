###--------------------------------------------------------------------------###

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
  min(x[best, "Variables"])
}

###--------------------------------------------------------------------------###

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
  min(x[best, "Variables"])
}

###--------------------------------------------------------------------------###