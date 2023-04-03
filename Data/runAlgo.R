###--------------------------------------------------------------------------###

# this is an ugly and long function
# but making a smaller one and than using a form of looping, e.g. (x)apply, for
# each algorithm I want to apply is the same, if not worse
# open to other suggestion though 

fitAlgo <- function(params, X, target){
  
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
  
  for (name in intersect(names(all), c("rfRFE", "linSvmRFE", "rbfSvmRFE"))){
    
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
    # these ranks are only when variables were selected
    # I want to obtain averaged ranking over 5 folds obtained before fitting 
    # each subset size
    #rk <- aggregate(resRFE$variables[,"Overall"], list(resRFE$variables$var), mean)
    #ranks[,name] <- rk[as.integer(sapply(colnames(X),function(x) which(x == rk[,1]))),"x"]
    acc[[name]] <- resRFE$results[resRFE$results$Variables==resRFE$bestSubset,"Accuracy"]
  }
  
  for (name in intersect(names(all), (c("rfSA", "rbfSvmRFE")))){
    
    ctrlSA <- caret::safsControl(functions = unlist(ifelse(grepl("^rf", name),
                                                           list(rfSA), 
                                                           list(caretSA))),
                                 method = "cv",number = params$kOut, 
                                 index = folds)
    resSA <- caret::safs(x = X, y = target,
                         iters = 50,
                         safsControl = ctrlSA,
                         method = all[[name]],
                         tuneGrid = NULL,
                         trControl = trControl)
    substs[[name]] <- resSA$optVariables
    acc[[name]] <- resSA$averages[resSA$optIter, "Accuracy"]
  } 
    
  #filt_eval <- FSinR::filterEvaluator('ReliefFeatureSetMeasure')
  #aco_search <- FSinR::antColony()
  #res <- aco_search(cbind(X, target), "target", featureSetEval = filt_eval)
  #preds <- dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]

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
  #ranks_stir <- stir::stir(X) maybe if not too slow
  
  # check how to add nominal: e.g. if you change name, ^.bin, ^.ord
  # and save for each threshold a subset 
  # FCBC_thres1
  # and so add the number of subsets for (n-1)+n_thresholds_fcbc
  lapply(params$fcbcThres, 
         function(x) Biocomb::select.fast.filter(cbind(X, target),
                                                      disc.method = "MDL",
                                                      threshold = x,
                                                      attrs.nominal = c(9:24, 57:120)))
  return(list(ranks, substs))
  
}

###--------------------------------------------------------------------------###

simCond(params, seed, rpt){
  
  # add magic with seeds
  X <- simX(params)
  Y <- simY(params, X[,grepl("^rel", colnames(X))])
  R <- factor(simR(params, X, Y), labels = c("c", "m"))
  res <- fitAlgo(params, X, R)
  saveRDS(res, paste0(""))
  
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

# make a loop over conditions (or rather rows of variedParams)
# and take care of the RNG -- simCond not even needed if I use foreach
# and then for over varied and within this

###--------------------------------------------------------------------------###