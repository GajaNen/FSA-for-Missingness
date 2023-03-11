###--------------------------------------------------------------------------###

fitAlgo <- function(func, name, params, folds, X, Y){
  
  #R2r <- NA
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
      # nam: coefficients of the final optimal model/HP fitted to the whole train set
      #R2r <- performance::r2_mckelvey(model$finalModel)
      nam <- names(which((coef(summary(model))[-1,"Pr(>|z|)"]) < params$pGLM))
      
    } else nam <- names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    
    indices <- extrReg(nam)
    # mean cv pf value
    acc <- model$results[rownames(model$bestTune), "Accuracy"]
    
  } else if (name %in% c("RF", "SVM")) {
    
    # change that to predefined structure
    preds <- rep(list(), 3)
    acc <- list()
    
    # RFE
    # to modify how rfe is fitted to the data modify selectVar function
    # also e.g. to remove the outer loop or not, to make it faster
    # modify other funcs
    # instead of guessing how is done, decide how you want to make it done
    # and modify the existing code
    ctrlRFE <- caret::rfeControl(functions = params[[paste0(name, "rfe")]], 
                                 method = "cv", number = params$kOut, 
                                 index = folds)
    resRFE <- caret::rfe(x = Xrf, y = Y,
                         sizes = params$sizes,
                         rfeControl = ctrlRFE,
                         trControl = trControl,
                         tuneGrid = tuneGrid,
                         method = func)
    preds$rfe <- predictors(resRFE)
    acc$rfe <- resRFE$results[resRFE$results$Variables==resRFE$bestSubset,"Accuracy"]
    
    # SA
    ctrlSA <- caret::safsControl(functions = params[[paste0(name, "sa")]], 
                                 method = "cv",number = params$kOut, 
                                 index = folds, improve = 50)
    resSA <- caret::safs(x = X, y = Y,
                         iters = 100,
                         safsControl = ctrlSA,
                         method = func,
                         tuneGrid = tuneGrid,
                         trControl = trControl)
    preds$sa <- resSA$optVariables
    acc$sa <- resSA$averages[resSA$optIter, "Accuracy"]
    # ACO
    filt_eval <- FSinR::filterEvaluator('ReliefFeatureSetMeasure')
    aco_search <- FSinR::antColony()
    res <- aco_search(cbind(Xrf, Y), "Y", featureSetEval = filt_eval)
    preds$ACO<- paste0("X",
                       extrReg(dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]),
                       "D")
    
    # Boruta for RF
    if (name == "RF") {
      # for "embedded method", e.g. 
      mod <- ranger::ranger(Y~., data=Xrf, importance = "impurity")
      # add accuracy of mod + Boruta 
      mod <- Boruta::Boruta(Y~X)
      predB <- names(mod$finalDecision[mod$finalDecision == "Confirmed"])
      if (length(predB)) preds$Boru <- paste0("X", predB, "D") else NULL
      # to get accuracy, a model with the variables selected by B
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
