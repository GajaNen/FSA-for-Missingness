###--------------------------------------------------------------------------###

sampleDT <- function(dt, N){
  
  dts <- data.table::setDT(lapply(seq_along(dt), function(x) numeric(N)))
  data.table::setnames(dts, names(dt))
  idx <- base::sample(1:nrow(dt),N,replace = T)
  dts[, names(dts) := dt[idx,]]
  return(dts)
}

###--------------------------------------------------------------------------###