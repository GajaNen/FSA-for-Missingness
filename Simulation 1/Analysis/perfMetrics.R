getPF <- function(conds, nrpt){
  
  path <- conds$dir
  subsets <- data.table::setDT(lapply(seq_along(c(conds$subsets, c(0, 0.1, 0.3, 0.5))),
                                      function(x) numeric(nrpt)))
  data.table::setnames(subsets, c(names(conds$subsets),paste0("FCBC_", c(0, 0.1, 0.3, 0.5))))
  #all <- c(names(conds$rankers),names(conds$subsets))
  all <- c(names(subsets)) # just for now when no code for rankings yet
  JI_sim <- array(NA, dim=c(length(all), length(all), nrpt))
  JI_stab <- data.table::setDT(lapply(seq_along(subsets), function(x) numeric(nrpt*(nrpt-1)/2)))
  data.table::setnames(JI_stab, names(subsets))
  colnames(JI_sim) <- rownames(JI_sim) <- names(subsets)
  if (conds$mechanism=="mcar"){
    pos <- NULL
    neg <- unlist( mapply(function(x,y) paste0(x, rep(c("Cont", "Bin", "Ord"),each=y/3),1:(y/3)),
                         c("Rel", "Irrel"), c(conds$Ntotal*conds$pr, conds$Ntotal*(1-conds$pr))),
                  use.names = F)
  } else  {
    pos <- paste0("Rel",rep(c("Cont", "Bin", "Ord"), each=conds$Ntotal*conds$pr/3),1:(conds$Ntotal*conds$pr/3))
    neg <- paste0("Irrel",rep(c("Cont", "Bin", "Ord"), each=conds$Ntotal*(1-conds$pr)/3),1:(conds$Ntotal*(1-conds$pr)/3))
    alpha <- min(0.5, length(pos)/length(neg))
  }
  for (i in 1:nrpt){
    res <- readRDS(file.path(path, paste0("mech_", conds$mechanism,
                                           "_pm_", conds$pm,
                                           "_corrPred_", conds$corrPred,
                                           "_pr_", conds$pr,
                                           "_rep_", i,
                                           ".RDS")))
    # check 
    if (any(names(res$res) != c("rankers", "subsets", "accuracies"))){
      stop(paste("Check results! There have been some warnings or errors in rep", i))
    }
    if (i == 1){ # we don't know in advance how many will have accuracy metrics
      acc <- data.table::setDT(lapply(1:length(res$res$accuracies), function(x) rep(0, nrpt)))
      data.table::setnames(acc, names(res$res$accuracies))
    }
    acc[i, names(acc) := res$res$accuracies]
    metrics <- lapply(res$res$subsets,
                      function(x){
                        TPR <- as.numeric(length(intersect(x, pos)) / length(pos))
                        FPR <- as.numeric(length(intersect(x, neg)) / length(neg))
                        return(c(TPR=TPR, FPR=FPR))
                        }
                      )
    if (conds$mechanism == "mcar"){
      # calculate PF which does not include pos: false positive rate = ME
      subsets[i, names(subsets) := lapply(metrics, function(x) 100*(1-x["FPR"]))]
    } else {
      # success rate TPR - FPR
      subsets[i, names(subsets) := lapply(metrics, function(x) 100*(x["TPR"] - alpha*x["FPR"]))]
    }
    # for rankers apply thresholds
    #res$res$rankers
    # power to select
    # similarity -- pairwise subsets Jaccard's index
    s <- res$res$subsets
    n <- length(s)
    similarity <- lapply(1:(n-1), 
                         function(x) sapply((x+1):n,
                                            function(y) {
                                              x <-
                                                as.numeric(length(intersect(s[x][[1]], s[y][[1]]))/
                                              length(union(s[x][[1]], s[y][[1]])))
                                              ifelse(is.nan(x), 100, 100*x)
                                              }
                                            ))
    
    m <- diag(1, n, n)
    m[lower.tri(m)] <- unlist(similarity)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    JI_sim[,,i] <- m
    diag(JI_sim[,,i]) <- 100
  }
  # stability -- Jaccard's index of all pairs of reps
  # return mean accuracy and all measures
  acc_means <- colMeans(acc)*100
  JI_sim_mean <- apply(JI_sim, c(1,2), mean)
  subsets_mean <- colMeans(subsets)
  return(list(accPred=acc_means,similarity=JI_sim_mean,PFM=subsets_mean))
}


applyThres <- function(rank, type){
  
  # apply thresholds 
}
