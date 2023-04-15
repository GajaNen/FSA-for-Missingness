###--------------------------------------------------------------------------###

# fit all algorithms at once

#@dat: must not contain incomplete variable Y, must contain all variables (rel&irrel)
# and missingness indicator named target

#return a list of
#feature rankings (DT), feature subsets by algo (list) and accuracies (list)
fitAlgo <- function(params, dat){
  
  ranks <- data.table::setDT(lapply(seq_along(params$rankers), rep, NA, params$Ntotal))
  data.table::setnames(ranks, names(params$rankers))
  substs <- vector(mode = "list", length = length(params$subsets))
  names(substs) <- names(params$subsets)
  acc <- list()
  all <- c(params$rankers, params$subsets) # parameters for each FSA
  folds <- caret::createFolds(dat[,target], k = params$kOut, returnTrain = T)
  # lasso SVM
  ids <- data.table::setDT(list(CVids=rep(0, params$N))) # get holdout indices instead of train
  for (x in 1:params$kOut) ids[setdiff(1:params$N, folds[[x]]), CVids := x] 
  res <- sparseSVM::cv.sparseSVM(X=as.matrix(dat[,!"target"]), 
                                 y=dat[,target], alpha = 1, gamma = 0.1, nlambda=100,
                                 lambda.min = 0.01, screen = "ASR", max.iter = 1000, eps = 1e-5,
                                 fold.id = ids, nfolds = params$kOut)
  lambda_res <- lambdaSE(res) # one SE rule (largest lambda with 1sdME within smallest)
  substs[["l1SVM"]] <- 
    names(res$fit$weights[-1,lambda_res$bestIndex])[res$fit$weights[-1,lambda_res$bestIndex] != 0]
  acc[["l1SVM"]] <- 1 - lambda_res$bestME
  trControl <- caret::trainControl(method = "cv", number = params$kInn,
                                   selectionFunction = "oneSE")
  #simulated annealing
  for (name in intersect(names(all), (c("rfSA", "linSvmSA", "rbfSvmSA")))){
    # specify controls, appropriate set of functions and method for the given name
    ctrlSA <- caret::safsControl(functions = unlist(ifelse(grepl("^rf", name),
                                                           list(rfSA), 
                                                           list(caretSA))),
                                 method = "cv", number = params$kOut, 
                                 index = folds, improve = 10)
    resSA <- caret::safs(x=dat[,!"target"], y=dat[,target],
                         iters = 10,
                         safsControl = ctrlSA,
                         method = all[[name]],
                         tuneGrid = NULL, # too long otherwise
                         trControl = trControl,
                         differences = F)
    substs[[name]] <- resSA$optVariables
    acc[[name]] <- resSA$averages[resSA$optIter, "Accuracy"]
  }
  binv <- grep(".*Bin",names(dat)) # if using SVM-RFE move RFE chunk before conversion to factors
  dat[, (binv) := lapply(.SD, factor, labels = c("z","o")), .SDcols = binv]
  #recursive feature elimination
  for (name in intersect(names(all), c("rfRFE", "linSvmRFE", "rbfSvmRFE"))){
    # specify controls, appropriate set of functions and method for the given name
    ctrlRFE <- caret::rfeControl(functions = unlist(ifelse(grepl("^rf", name),
                                                           list(rfFuncs), 
                                                           list(caretFuncs))), 
                                 method = "cv", number = params$kOut, 
                                 index = folds)
    resRFE <- caret::rfe(target~., data = dat,
                         sizes = params$sizes,
                         rfeControl = ctrlRFE,
                         trControl = trControl,
                         tuneGrid = params$tuneGrids[[all[[name]]]],
                         method = all[[name]])
    substs[[name]] <- caret::predictors(resRFE)
    acc[[name]] <- resRFE$results[resRFE$results$Variables==resRFE$bestSubset,"Accuracy"]
  }
  # regressions
  for (name in intersect(names(all), c("lasso", "EN", "LR"))){
    
    trControl <- caret::trainControl(method = "cv", number = params$kOut, 
                                     classProbs = T, index = folds,
                                     selectionFunction = "oneSE")
    model <- caret::train(target~., data = dat, family = "binomial",
                          method = all[[name]], trControl = trControl, 
                          tuneGrid = params$tuneGrids[[name]], 
                          metric = "Accuracy")
    if (name == "LR") {
      ranks[, (name) := coef(summary(model))[-1,"Pr(>|z|)"]]
    } else substs[[name]] <- 
      names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    acc[[name]] <- model$results[rownames(model$bestTune), "Accuracy"]
  }
  #hybrid: gini (filter) + whale opt (wrapper)
  filt_eval <- FSinR::filterEvaluator('giniIndex')
  hyb_search <- FSinR::whaleOptimization(10, 50)
  res <- hyb_search(dat, "target", featureSetEval = filt_eval)
  substs[["hyb"]] <- dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]
  # Boruta & ReliefF
  ranks[, Boruta := (Boruta::Boruta(target~., data = dat))$finalDecision]
  res <- ranger::ranger(target~., data=dat, importance = "impurity",splitrule = "gini",
                        oob.error = T)
  acc[["RF"]] <- 1 - res$prediction.error
  ranks[, RF := res$variable.importance]
  ranks[, ReliefF := FSelector::relief(target~., data=dat,sample.size = 10)[,1]]
  substs <- c(substs, setNames(lapply(params$fcbcThres, #fast correlation-based filter
                                      function(x) Biocomb::select.fast.filter(
                                        dat,#last col must be the target
                                        disc.method = "MDL",
                                        threshold = x,
                                        attrs.nominal = binv)),
                               paste0("FCBC_", params$fcbcThres))) # apply different thresholds
  return(list(rankers=ranks, subsets=substs, accuracies=acc))
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
  if (!rpt %in% .lec.GetStreams()) rlecuyer::.lec.CreateStream(1:fixed$streams)
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
    res <- tryCatch(withCallingHandlers(
      expr = fitAlgo(conds, XY[,!"Y"]), 
      warning = function(w) {
        warns <<- c(warns, list(w))
        invokeRestart("muffleWarning")
      }), 
      error = function(e) {
        errs <<- c(errs, list(e))
      }
    )
    if (!dir.exists(conds$dir)) dir.create(conds$dir)
    saveRDS(list(coef=mssng$coefs,
                 res=res,
                 RNGstate=s,
                 changes=changes),
            file.path(conds$dir,paste0("mech_", conds$mechanism,
                                       "_pm_", conds$pm,
                                       "_corrPred_", conds$corrPred,
                                       "_pr_", conds$pr,
                                       "_rep_", rpt,
                                       ".RDS")))
    prev <- varied[i, 1:nfac]
    # cat(paste0("Condition ", i, "out of ", nrow(varied),
    #             "in repetition ", rpt, "finished!"))
  }
}

###--------------------------------------------------------------------------###
