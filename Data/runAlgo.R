###--------------------------------------------------------------------------###

fitAlgo <- function(func, name, params, folds, X, Y){
  
  #R2r <- NA
  regs <- name %in% c("lasso", "EN", "LR")
  kTr <- regs * params$kOut + (!regs) * params$kInn
  nTR <- 1 + (name == 'EN')
  if (kTr == params$kOut) foldsTr <- folds else foldsTr <- NULL
  tuneGrid <- params$tuneGrids[[name]]
  trControl <- caret::trainControl(method = "cv", number = kTr, 
                                   classProbs = T, index = foldsTr)
  if (regs){
    
    model <- caret::train(x = X, y = Y, family = "binomial",
                          method = func, trControl = trControl, 
                          tuneGrid = tuneGrid)
    if (name == "LR") {
      # nam: coefficients of the final optimal model/HP fitted to the whole train set
      #R2r <- performance::r2_mckelvey(model$finalModel)
      nam <- names(which((coef(summary(model))[-1,"Pr(>|z|)"]) < params$pGLM))
      
    } else nam <- names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    
    indices <- extrReg(nam)
    # mean cv pf value
    acc <- model$results[rownames(model$bestTune), "Accuracy"]
    
  } else {
    
    # change that to predefined structure
    preds <- rep(list(), 3)
    acc <- list()
    
    # RFE
    ctrlRFE <- caret::rfeControl(functions = params[[paste0(name, "rfe")]], 
                                 method = "cv", number = params$kOut, 
                                 index = folds)
    resRFE <- caret::rfe(x = X, y = Y,
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
    
    # Boruta for RF
    if (name == "RF") {
      mod <- Boruta(Y~X)
      predB <- names(mod$finalDecision[mod$finalDecision == "Confirmed"])
      if (length(predB)) preds$Boru <- paste0("X", predB, "D") else NULL
      # mean cross-validated pf
    }
    
  }
  
  return(list(acc=acc, preds=preds, R2r=R2r))
  
}

###--------------------------------------------------------------------------###
