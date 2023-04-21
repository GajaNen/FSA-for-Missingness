###--------------------------------------------------------------------------###

getPF <- function(conds, nrpt, nms){
  
  path <- conds$dir
  errs <- 0
  
  if (conds$mechanism=="mcar"){
    pos <- NULL
    neg <- nms
  } else  {
    pos <- nms[grep("^Rel",nms)]
    neg <- nms[grep("^Irrel",nms)]
    alpha <- min(0.5, length(pos)/length(neg))
  }
  for (i in 1:nrpt){
    res <- readRDS(file.path(path, paste0("mech_", conds$mechanism,
                                           "_pm_", conds$pm,
                                           "_corrPred_", conds$corrPred,
                                           "_pr_", conds$pr,
                                           "_rep_", i,
                                           ".RDS")))
    if (!identical(names(res$res),c("rankers", "subsets", "accuracies"))){
      errs <- errs+1
    } else {
      rnkrs_thr <- applyThres(res$res$rankers, nms)
      rnkrs_sub <- rnkrs_thr$res
      all_substs <- c(res$res$subsets, rnkrs_sub)
      if (i == 1){ # we don't know in advance how many will have accuracy metrics
        acc <- data.table::setDT(lapply(1:length(res$res$accuracies), function(x) rep(0, nrpt)))
        data.table::setnames(acc, names(res$res$accuracies))
        subsets <- data.table::setDT(lapply(seq_along(all_substs),
                                            function(x) numeric(nrpt)))
        data.table::setnames(subsets, c(names(all_substs)))
        JI_sim <- array(NA, dim=c(length(all_substs), length(all_substs), nrpt))
        JI_stab <- data.table::setDT(lapply(seq_along(all_substs), function(x) 0))
        data.table::setnames(JI_stab, names(all_substs))
        colnames(JI_sim) <- rownames(JI_sim) <- names(all_substs)
        pwr <- rep(0, length(all_substs))
        alls <- lapply(1:nrpt,
                       function(x) vector(mode="list", length = length(all_substs)))
      }
      alls[[i]] <- all_substs
      acc[i, names(acc) := res$res$accuracies]
      metrics <- lapply(all_substs,
                        function(x){
                          if (class(x)=="try-error"){
                            x <- NULL
                          }
                          else if(class(x)=="data.frame"){
                            x <- rownames(x)
                          }
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
        pwr <- pwr + (unlist(lapply(metrics, "[",1)) == 1)
      }
      # similarity -- pairwise subsets Jaccard's index
      s <- all_substs
      n <- length(s)
      similarity <- lapply(1:(n-1),function(x) sapply((x+1):n, function(y) Jaccard(s[x][[1]], s[y][[1]])))
      
      m <- diag(100, n, n)
      m[lower.tri(m)] <- unlist(similarity)
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      JI_sim[,,i] <- m
      }
    }
  # stability -- Jaccard's index of all pairs of reps
  JI_stab[, names(JI_stab) :=
            lapply(names(JI_stab), 
                   function(x) mean(unlist(sapply(1:(nrpt-1),
                                      function(y)
                                        sapply((y+1):nrpt,
                                      function(n) Jaccard(alls[[y]][[x]],
                                                  alls[[n]][[x]]))
                                      )))
                   )
          ]
  # return mean accuracy and all measures
  acc_means <- colMeans(acc)*100
  acc_sd <- acc[, lapply(.SD, sd)]
  JI_sim_mean <- apply(JI_sim, c(1,2), mean)
  subsets_mean <- colMeans(subsets)
  return(list(accMean=acc_means,accSD=acc_sd,
              similarity=JI_sim_mean,
              PFM=subsets_mean,power=pwr/nrpt,
              stability=JI_stab, rnksM=rnkrs_thr$means,
              rnksSD=rnkrs_thr$sds))
}

###--------------------------------------------------------------------------###

Jaccard <- function(x,y){
  ji <- as.numeric(length(intersect(x, y)) / length(union(x, y)))
  return(ifelse(is.nan(ji), 100, 100*ji))
}

###--------------------------------------------------------------------------###

applyThres <- function(rk, nms){
  
  N <- nrow(rk)
  rs <- list()
  t.LR <- c(0.1, 0.05, 0.01, 0.001)
  ## general:
  rs[paste0("LR_", t.LR)] <- lapply(t.LR, function(x) nms[rk[,LR < x]])
  rs[["Boruta_conf"]] <- nms[rk[, Boruta == "Confirmed"]]
  rs[["Boruta_tent"]] <- nms[rk[, Boruta == "Tentative"]]
  # log2(n)
  t.log2n <- ceiling(log2(N))
  rs[paste0(c("ReliefF", "RF"), "_log2n")] <-
    lapply(rk[, c("ReliefF", "RF")], function(x) nms[order(x)[1:t.log2n]])
  # Fisher's dicriminant ratio
  
  # ReliefF
  #rk$ReliefF
  #seq(0, 1/sqrt(0.05*10), 0.2)
  t.relief <- c(0, 1/sqrt(0.05*10), mean(rk$ReliefF)-2*var(rk$ReliefF))
  rs[paste0("RelieF_", c("0", "alphaM", "meanVar"))] <-
    lapply(t.relief, function(x) nms[rk[, ReliefF > x]])
  # or take the largest gap as cut-off
  ord.relief <- rk$ReliefF[order(rk$ReliefF, decreasing = T)]
  ord.nms.rel <- nms[order(rk$ReliefF, decreasing = T)]
  # min & negative diffs -- in case of 0 diff, this would be selected
  rs[["ReliefF_gap"]] <- ord.nms.rel[0:which.max(abs(diff(ord.relief)))]
  t.perc <- c(0.05, 0.1, 0.2, 0.5)
  ord.nms.rf <- nms[order(rk$RF, decreasing = T)]
  rs[paste0("RF_", t.perc)] <- lapply(t.perc, function(x) ord.nms.rf[1:(N*x)])
  rs[paste0("ReliefF_", t.perc)] <- lapply(t.perc, function(x) ord.nms.rel[1:(N*x)])
  mean_RR <- rk[, c(mean(ReliefF), mean(RF))]
  sd_RR <- rk[, c(sd(ReliefF), sd(RF))]
  names(mean_RR) <- names(sd_RR) <- c("ReliefF", "RF")
  return(list(res=rs,means=mean_RR,sds=sd_RR))
}


###--------------------------------------------------------------------------###