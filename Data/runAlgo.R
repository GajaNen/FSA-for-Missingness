###--------------------------------------------------------------------------###

fitAlgo <- function(params, X, target, rpt){
  
  ranks <- data.table::data.table(array(
    data = NA, dim = c(ncol(X), length(params$rankers))))
  colnames(ranks) <- names(params$rankers)
  substs <- vector(mode = "list", length = length(params$subsets))
  names(substs) <- names(params$subsets)
  acc <- list()
  all <- c(params$rankers, params$subsets)
  folds <- createFolds(target, k = params$kOut, returnTrain = TRUE)
  
  for (name in intersect(names(all), c("lasso", "EN", "LR"))){
    
    trControl <- caret::trainControl(method = "cv", number = params$kOut, 
                                     classProbs = T, index = folds,
                                     selectionFunction = "oneSE")
    model <- caret::train(x = X, y = target, family = "binomial",
                          method = all[[name]], trControl = trControl, 
                          tuneGrid = params$tuneGrids[[name]], 
                          metric = "Accuracy")
    if (name == "LR") {
      ranks[,name] <- coef(summary(model))[-1,"Pr(>|z|)"]
      } else substs[[name]] <- 
      names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    acc[[name]] <- model$results[rownames(model$bestTune), "Accuracy"]
  }
  
  trControl <- caret::trainControl(method = "cv", number = params$kInn,
                                   selectionFunction = "oneSE")
  
  for (name in intersect(names(all), c("rfRFE", "rbfSvmRFE"))){
    
    ctrlRFE <- caret::rfeControl(functions = unlist(ifelse(grepl("^rf", name),
                                                    list(rfFuncs), 
                                                    list(caretFuncs))), 
                                 method = "cv", number = params$kOut, 
                                 index = folds)
    resRFE <- caret::rfe(x = X, y = target,
                         sizes = params$sizes,
                         rfeControl = ctrlRFE,
                         trControl = trControl,
                         tuneGrid = params$tuneGrids[[all[[name]]]],
                         method = all[[name]])
    substs[[name]] <- predictors(resRFE)
    acc[[name]] <- resRFE$results[resRFE$results$Variables==resRFE$bestSubset,"Accuracy"]
  }
  
  for (name in intersect(names(all), (c("rfSA", "linSvmSA", "rbfSvmSA")))){
    
    ctrlSA <- caret::safsControl(functions = unlist(ifelse(grepl("^rf", name),
                                                           list(rfSA), 
                                                           list(caretSA))),
                                 method = "cv", number = params$kOut, 
                                 index = folds, improve = 10)
    resSA <- caret::safs(x = X, y = target,
                         iters = 10,
                         safsControl = ctrlSA,
                         method = all[[name]],
                         tuneGrid = NULL, # too long otherwise
                         trControl = trControl,
                         differences = F)
    substs[[name]] <- resSA$optVariables
    acc[[name]] <- resSA$averages[resSA$optIter, "Accuracy"]
  } 
  
  filt_eval <- FSinR::filterEvaluator('giniIndex')
  hyb_search <- FSinR::whaleOptimization(10, 50)
  res <- hyb_search(cbind(X, target), "target", featureSetEval = filt_eval)
  substs[["hyb"]] <- dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]

  ids <- rep(NA, params$N)
  for (x in 1:params$kOut) ids[setdiff(1:params$N, folds[[x]])] <- x
  res <- sparseSVM::cv.sparseSVM(as.matrix(X), target, alpha = 1, gamma = 0.1, nlambda=100,
                                 lambda.min = 0.01, screen = "ASR", max.iter = 1000, eps = 1e-5,
                                 fold.id = ids, nfolds = params$kOut)
  lambda_res <- lambdaSE(res)
  substs[["l1SVM"]] <- 
    names(res$fit$weights[-1,lambda_res[[1]]])[res$fit$weights[-1,lambda_res[[1]]] != 0]
  acc[["l1SVM"]] <- 1 - lambda_res[[2]]
  
  ranks[,"Boruta"] <- (Boruta::Boruta(target~X))$finalDecision
  
  # check sample size parameter
  ranks[,"ReliefF"] <- FSelectorRcpp::relief(x=X, y=target, sampleSize = 10)[,2]
  
  substs <- c(substs, setNames(lapply(params$fcbcThres, 
                                      function(x) Biocomb::select.fast.filter(
                                        cbind(X, target),
                                        disc.method = "MDL",
                                        threshold = x,
                                        attrs.nominal = which(grepl("(.*Bin)|(.*Ord)", 
                                                                    colnames(X))))),
                                paste0("FCBC_", params$fcbcThres)))
  
  saveRDS(list(ranks, substs, acc), paste0("mech_", params$mechanism,
                                           "_pm_", params$pm,
                                           "_corrPred_", params$corrPred,
                                           "_pr_", params$pr,
                                           "_rep_", rpt))
}

###--------------------------------------------------------------------------###
              
simRep <- functio(fixed, varied, seed, rp){
  
  for (i in 1:nrow(varied)){
    conds <- c(fixed, varied[i,])
    # if change in X / Y simulate else keep the ones from before
    X <- simX(params)
    Y <- simY(params, X[,grepl("^rel", colnames(X))])
    R <- factor(simR(params, X, Y), labels = c("c", "m"))
    res <- fitAlgo(params, X, R)
    saveRDS(res, paste0(""))
  }
}

###--------------------------------------------------------------------------###
