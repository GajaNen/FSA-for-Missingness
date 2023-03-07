###--------------------------------------------------------------------------###

fitAlgo <- function(func, name, params, folds, X, Y){
  
  R2r <- NA
  regs <- name %in% c("lasso", "EN", "LR")
  kTr <- regs * params$kOut + (!regs) * params$kInn
  nTR <- 1 + (name == 'EN')
  if (kTr == params$kOut) foldsTr <- folds else foldsTr <- NULL
  tuneGrid <- params$tuneGrids[[name]]
  trControl <- caret::trainControl(method = "cv", number = kTr, 
                                   classProbs = T, index = foldsTr)
  if (regs){
    
    model <- caret::train(x=X, y=Y, family = "binomial",
                          method = func, trControl = trControl, 
                          tuneGrid = tuneGrid, tuneLength = 40)
    if (name == "LR") {
      # nam: coefficients of the final optimal model/HP fitted to the whole train set
      R2r <- performance::r2_mckelvey(model$finalModel)
      nam <- names(which((coef(summary(model))[-1,"Pr(>|z|)"]) < params$pGLM))
      
    } else nam <- names(which(coef(model$finalModel, model$bestTune$lambda)[-1,] != 0))
    
    indices <- extrReg(nam)
    # mean cv pf value
    acc <- model$results[rownames(model$bestTune), "Accuracy"]
    
  } else {
    
    ctrl <- caret::rfeControl(functions = params[[name]], method = "cv",
                              number = params$kOut, index = folds)
    res <- caret::rfe(x=X, y=Y,
                      sizes = params$sizes,
                      rfeControl = ctrl,
                      trControl = trControl,
                      tuneGrid=tuneGrid,
                      method = func)
    indices <- extrReg(predictors(res))
    # mean cross-validated pf
    acc <- res$results[res$results$Variables==res$bestSubset,"Accuracy"]
    
  }
  
  preds <- if (length(indices)) paste0("X", indices, "D") else NULL
  return(list(acc=acc, preds=preds, R2r=R2r))
  
}

###--------------------------------------------------------------------------###
