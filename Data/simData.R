###--------------------------------------------------------------------------###

simX <- function(params){
  
  # generate relevant independent X
  Nrel <- params$percRel * params$Ntotal
  Nrest <- params$Ntotal - Nrel
  if ((Nrel %% 3) || (Nrest %% 3)) stop(cat("Number of relevant and
                                        the rest of features must
                                        be a multiple of 3."))
  relt <- cbind(sapply(params$normParam, 
                      function(x) rnorm(params$N, x[1], x[2]),
                      simplify = T),
               sapply(params$gammaParam, 
                      function(x) rgamma(params$N, shape=x[1], rate=x[2]),
                      simplify = T),
               sapply(params$probBinRel, 
                      function(x) rbinom(params$N, 1, x),
                      simplify = T),
               sapply(1:(Nrel / 3), 
                      function(x) sample(0:(params$ncatsRel[x]-1), 
                                         params$N, 
                                         replace = TRUE, 
                                         prob = params$probOrdRel[[x]]),
                      simplify = T))
  colnames(relt) <- paste0("relevan", 1:ncol(relt))
  # generate redundant features
  if (params$typeRest == "IR"){
    Nred <- Nrel
    redt <- relt # add names
    # add noise 
    colnames(redt) <- paste0("redund", 1:ncol(redt))
  } else {
    Nred <- 0
    redt <- vector(length = params$N)
  }
  # generate irrelevant
  irrelt <- cbind(sapply(params$normParam, 
                        function(x) rnorm(params$N, x[1], x[2]),
                        simplify = T),
                 sapply(params$gammaParam, 
                        function(x) rgamma(params$N, shape=x[1], rate=x[2]),
                        simplify = T),
                 sapply(params$probBinIrrel, 
                        function(x) rbinom(params$N, 1, x),
                        simplify = T),
                 sapply(1:(Nrest - Nred), 
                        function(x) sample(0:(params$ncatsIrrel[x]-1), 
                                           params$N, 
                                           replace = TRUE, 
                                           prob = params$probOrdIrrel[[x]]),
                        simplify = T))
  colnames(irrelt) <- paste0("irrelvan", 1:ncol(irrelt))
  
  # is sample the correct way? i.e. are probs the cumulative probs
  
  return(data.frame(cbind(relt, irrelt, redt)))
  
}

###--------------------------------------------------------------------------###

simY <- function(params, X){
  
  rey <- cor(X[,1:(params$percRel * params$Ntotal)], 
             X[,1:(params$percRel * params$Ntotal)])
  cors <- sample(1:10, (params$percRel * params$Ntotal), replace = T)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1]
  cors <- a*cors
  # for given correlations simulate Y
  X <- scale(X[,1:(params$percRel * params$Ntotal)]) 
  x <- 1:params$N
  e <- residuals(lm(x ~ X)) 
  X.dual <- with(svd(X), (params$N-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - cors %*% cov(X.dual) %*% cors) / var(e))
  return(X.dual %*% cors + sqrt(sigma2)*e)
}

###--------------------------------------------------------------------------###

simRho <- function(params){
  
  # sample from LKJ distribution
  rey <- rethinking::rlkjcorr(1, params$nX, params$eta)
  diag(rey) <- 1
  
  # move cat y to the beginning: requirement for simmulticorr package
  vars_move <- paste0("Y", which(params$typeY == "cat"), "D")
  preds <- rep(list(rep(NA, sum(params$infoX)/params$nY)), params$nY)
  
  for (j in 1:params$nY){
    
    corX <- rep(0, params$nX)
    names(corX) <- paste0("X", 1:params$nX, "D")
    
    nCat <- params$infoX[1]/params$nY
    nCont <- params$infoX[2]/params$nY
    end_non <- params$infoX[1]+params$infoXnon[1]
    
    # each Y has different X predictors to get different patterns
    preds[[j]] <- paste0("X", 
                         c((1:(nCat)) + (nCat)*(j-1), 
                           ((end_non+1):(end_non+nCont)) + (nCont)*(j-1)), 
                         "D")
    corX[preds[[j]]] <- runif(length(preds[[j]]), 4, 6)
    # each y_i uncorrelated with every y_j
    cors <- c(corX, rep(0, j-1))
    n <- params$nX+j-1
    # generate corr between X and Y for a given R2 (defined in terms of deviance)
    a <- sqrt(params$R2y/(t(cors) %*% solve(rey[1:n,1:n]) %*% cors))[1]
    cors <- a*cors
    rey <- cbind(rbind(rey, cors), c(t(cors), 1))
    
    colnames(rey)[params$nX+j] <- paste0("Y", j, "D")
    
    
    
  }
  
  rownames(rey) <- colnames(rey)
  rey <- reorder_mat(rey, c(vars_move, setdiff(colnames(rey), vars_move)))
  
  return(list(rey=rey, preds=preds))
  
}


###--------------------------------------------------------------------------###
