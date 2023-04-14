###--------------------------------------------------------------------------###

sampleDT <- function(dt, N){
  
  dts <- data.table::setDT(lapply(seq_along(dt), function(x) numeric(N)))
  data.table::setnames(dts, names(dt))
  idx <- base::sample(1:nrow(dt),N,replace = T)
  dts[idx, names(dts) := dt[idx,]]
  return(dts)
}

# should it be stratified on outcome?

###--------------------------------------------------------------------------###