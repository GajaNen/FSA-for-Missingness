###--------------------------------------------------------------------------###

fitAlgo <- function(func, name, params, folds, X, Y){
  
  regs <- name %in% c("lasso", "EN", "LR")
  kTr <- regs * params$kOut + (!regs) * params$kInn
  nTR <- 1 + (name == 'EN')
  if (kTr == params$kOut) foldsTr <- folds else foldsTr <- NULL
  tuneGrid <- params$tuneGrids[[name]]
  trControl <- caret::trainControl(method = "cv", number = kTr, 
                                   classProbs = T, index = foldsTr,
                                   selectionFunction = "oneSE")
  if (regs){
    
    model <- caret::train(x = X, y = Y, family = "binomial",
                          method = func, trControl = trControl, 
                          tuneGrid = tuneGrid, metric = "Accuracy")
    if (name == "LR") {
      nam <- names(which((coef(summary(model))[-1,"Pr(>|z|)"]) < params$pGLM))
    } else nam <- names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    
    indices <- if (length(name)) extrReg(nam) else NULL
    # mean cv pf value
    acc <- model$results[rownames(model$bestTune), "Accuracy"]
    
  } else if (name %in% c("RF-RFE", "SVM-RFE")) {
    
    # RFE
    ctrlRFE <- caret::rfeControl(functions = params[[paste0(name, "rfe")]], 
                                 method = "cv", number = params$kOut, 
                                 index = folds)
    resRFE <- caret::rfe(x = Xrf, y = Y,
                         sizes = params$sizes,
                         rfeControl = ctrlRFE,
                         trControl = trControl,
                         tuneGrid = tuneGrid,
                         method = func)
    preds <- predictors(resRFE)
    acc <- resRFE$results[resRFE$results$Variables==resRFE$bestSubset,"Accuracy"]
  } else if (name %in% c("RF-SA", "SVM-SA")){
    
    ctrlSA <- caret::safsControl(functions = params[[paste0(name, "sa")]], 
                                 method = "cv",number = params$kOut, 
                                 index = folds, improve = 50)
    resSA <- caret::safs(x = X, y = Y,
                         iters = 100,
                         safsControl = ctrlSA,
                         method = func,
                         tuneGrid = tuneGrid,
                         trControl = trControl)
    preds <- resSA$optVariables
    acc <- resSA$averages[resSA$optIter, "Accuracy"]
  } else if (name == "ACO"){
    
    filt_eval <- FSinR::filterEvaluator('ReliefFeatureSetMeasure')
    aco_search <- FSinR::antColony()
    res <- aco_search(cbind(X, Y), "Y", featureSetEval = filt_eval)
    preds <- paste0("X",
                       extrReg(dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]),
                       "D")
  } else if (name == "SVM-L1"){
    
    ids <- rep(NA, 1:params$N)
    for (x in 1:params$kOut) ids[setdiff(1:params$N, folds[[x]])] <- x
    res <- sparseSVM::cv.sparseSVM(as.matrix(XsvmDumm), Y, alpha = 1, gamma = 0.1, nlambda=100,
                                  lambda.min = 0.01, screen = "ASR", max.iter = 1000, eps = 1e-5,
                                  fold.id = ids, nfolds = params$kOut)
    lambda_res <- lambdaSE(res)
    predsSSVM <- extrReg(
      names(res$fit$weights[-1,lambda_res[[1]]])[res$fit$weights[-1,lambda_res[[1]]] != 0])
    preds <- if (length(predsSSVM)) paste0("X", predsSSVM, "D") else NULL
    acc <- 1 - lambda_res[[2]]
  } else (name == "Boruta"){
    
    if (name == "RF") {
      # add accuracy of mod + Boruta 
      mod <- Boruta::Boruta(Y~X)
      predB <- names(mod$finalDecision[mod$finalDecision == "Confirmed"])
      if (length(predB)) preds <- paste0("X", extrReg(predB), "D") else NULL
    }
    
  } else if (name == "ReliefF"){
    
    ranks <- FSelectorRcpp::relief(x=X, y=Y, sampleSize = 50) # works with factors
    # consider only those which are not smaller than mean - 2*var
    preds <- 
      ranks[ranks$importance > (mean(ranks$importance) - 2*var(ranks$importance)),
            "attributes"]
    # or take the largest gap as cut-off
    ord.ranks <- ranks[order(ranks$importance, decreasing = T),]
    # min & negative diffs -- in case of 0 diff, this would be selected
    preds <- ord.ranks[0:which.max(abs(diff(ord.ranks$importance))), "attributes"]
  } 
  
  return(list(acc=acc, preds=preds, R2r=R2r))
  
}

###--------------------------------------------------------------------------###
