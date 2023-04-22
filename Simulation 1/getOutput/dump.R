# genProbs <- function(nVar, ncat=rep(5, nVar)){
#   
#   for (i in 1:nVar){
#     
#     sm <- 0
#     eq <- rep(1/ncat[i], ncat[i])
#     cm <- cumsum(eq)
#     lapply(cm, function(x) 1-x)
#     
#   }
# }