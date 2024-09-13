# dependencies: data.table, FCBF, rlecuyer

###--------------------------------------------------------------------------###

# fit all algorithms at once

# @dat: must not contain incomplete variable Y, must contain all variables (rel&irrel)
# and missingness indicator named target

#return a list of
#feature rankings (DT), feature subsets by algo (list) and accuracies (list)
fitAlgo <- function(params, dat){
  
  substs <- setNames(lapply(params$fcbcThres, #fast correlation-based filter
                                      function(x) 
                                        try(FCBF::fcbf(FCBF::discretize_exprs(t(scale(dat[,!"target"]))),
                                                       t(dat[,target]),
                                                       minimum_su=x), 
                                            silent = T)),  
                               paste0("FCBC_", params$fcbcThres)) # apply different thresholds
  
  return(list(subsets=substs))
}

###--------------------------------------------------------------------------###

# execute one repetition of all conditions

#@fixed: a list of fixed parameters, which don't change over conditions
#@varied: a list of varied parameters, which change over conditions
#@for more information about parameters which have to be contained check setup.R
#@rpt: number of repetition
#@nfac: number of varied (primary) factors

simRep <- function(fixed, varied, rpt="test", nfac=4){
  
  rlecuyer::.lec.SetPackageSeed(rep(fixed$seed, 6))
  if (!rpt %in% rlecuyer::.lec.GetStreams()) rlecuyer::.lec.CreateStream(1:fixed$streams)
  rlecuyer::.lec.CurrentStream(rpt)
  
  prev <- rep("none", nfac)
  for (i in 1:nrow(varied)){
    warns <- errs <- list()
    conds <- c(fixed, varied[i,])
    changes <- varied[i, 1:nfac] != prev 
    names(changes) <- colnames(varied)[1:nfac]
    s <- .Random.seed
    if (changes["pr"] || changes["corrPred"]) XY <- simDat(conds)
    mssng <- simR(conds, XY)
    XY[, target := factor(mssng$R, labels = c("c","m"))]
    res_fcbf <- tryCatch(withCallingHandlers(
        expr = fitAlgo(conds, XY[,!"Y"]), 
        warning = function(w) {
          warns <<- c(warns, list(w))
          invokeRestart("muffleWarning")
        }), 
        error = function(e) {
          errs <<- c(errs, list(e))
        }
      )
    
    # open the previous results and save them
    res <- readRDS(file.path(conds$prevdir,paste0("mech_", conds$mechanism,
                                                   "_pm_", conds$pm,
                                                   "_corrPred_", conds$corrPred,
                                                   "_pr_", conds$pr,
                                                   "_rep_", rpt,
                                                   ".RDS")))
    # then check for fcbf and replace them
    FCBF.in.nms <- grepl("FCBC_*", names(res$res$subsets))
    if (any(FCBF.in.nms) && conds$simN == 1){
      
      for (name in names(res$res$subsets)[FCBF.in.nms]){
        
        res$res$subsets[[name]] <- res_fcbf$subsets[[name]]
        
      }
      
    }
    
    # save res object in the new dir
    saveRDS(list(coef_rerun=mssng$coefs,
                 res_rerun=res,
                 RNGstate_rerun=s,
                 changes_rerun=changes),
            file.path(conds$dir,paste0("mech_", conds$mechanism,
                                       "_pm_", conds$pm,
                                       "_corrPred_", conds$corrPred,
                                       "_pr_", conds$pr,
                                       "_rep_", rpt,
                                       ".RDS")))
    
    prev <- varied[i, 1:nfac]
  }
}

###--------------------------------------------------------------------------###