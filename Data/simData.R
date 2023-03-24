###--------------------------------------------------------------------------###

simX <- function(params){
  
  # generate relevant independent X
  relt <- cbind(sapply(params$normParam, 
                      function(x) rnorm(params$N, x[1], x[2]),
                      simplify = T),
               sapply(params$gammaParam, 
                      function(x) rgamma(params$N, shape=x[1], rate=x[2]),
                      simplify = T),
               sapply(params$probBinRel, 
                      function(x) rbinom(params$N, 1, x),
                      simplify = T),
               sapply(1:params$Nrel[3], 
                      function(x) sample(0:(params$ncatsRel[x]-1), 
                                         params$N, 
                                         replace = TRUE, 
                                         prob = params$probOrdRel[[x]]),
                      simplify = T))
  colnames(relt) <- paste0("relevan", 1:ncol(relt))
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
                 sapply(1:params$Nirrel[3], 
                        function(x) sample(0:(params$ncatsIrrel[x]-1), 
                                           params$N, 
                                           replace = TRUE, 
                                           prob = params$probOrdIrrel[[x]]),
                        simplify = T))
  colnames(irrelt) <- paste0("irrelvan", 1:ncol(irrelt))
  # generate redundant
  redt <- relt # add names
  colnames(redt) <- paste0("redund", 1:ncol(redt))
  return(cbind(relt, irrelt, redt))
  
}

###--------------------------------------------------------------------------###

simY <- function(params, X){
  
  rey <- cor(X[,1:sum(params$Nrel)], X[,1:sum(params$Nrel)])
  cors <- sample(1:5, sum(params$Nrel), replace = T)
  a <- sqrt(params$R2y/(t(cors) %*% solve(rey) %*% cors))[1]
  cors <- a*cors
  # for given correlations simulate Y (cont or transformed to bin/cat)
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


simM <- function(params){
  
  # marginal distributions for cont variables
  M <- t(apply(params$contParam, 1, 
               function(x) round(calc_theory(Dist=x[3], 
                                             params=as.numeric(x[1:2])), 
                                 8)))
  
  return(M)
}

###--------------------------------------------------------------------------###

simDat <- function(params, seed){
  
  M <- simM(params)
  res <- simRho(params)
  
  B <- SimMultiCorrData::rcorrvar(n = params$N, k_cont = params$ncont, 
                                  k_cat = params$ncat, method = "Polynomial", 
                                  means = M[,1], vars =  (M[,2])^2,
                                  skews = M[,3], skurts = M[,4],
                                  fifths = M[,5], sixths = M[,6], 
                                  marginal = params$marginal, rho = res$rey, 
                                  errorloop = FALSE, seed=seed)
  

  # separate X and Y and name their columns
  X <- cbind(B$ordinal_variables[,-c(1:params$ycat)] - 1, 
             B$continuous_variables[,1:params$xcont])
  colnames(X) <- paste0("X", 1:params$nX, "D")
  Y <- cbind(B$ordinal_variables[,1:params$ycat] - 1, 
             B$continuous_variables[,-c(1:params$xcont)])
  colnames(Y) <- paste0("Y", 1:params$nY, "D")
  
  return(list(X=X, Y=Y, preds=res$preds, rey=res$rey))
  
  
}

###--------------------------------------------------------------------------###
