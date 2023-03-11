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
  
  # check if index does what it should do (return indices of the elements of x
  # with this criteria)
  # if it returns indices, then the smallest index will always refer to the smallest
  # HP or smallest SS
  # and then select of the column variables the one with such index
  # could also rewrite it in a diff way (by indexing directly the column Variables
  # candidates <- x[x[, metric] >= perf, "Variables"]
  # min(candidates)
  # check if it differs or rather if the current way always works
  # if the code for train does not work do the same
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