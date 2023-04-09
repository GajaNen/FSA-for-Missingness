###--------------------------------------------------------------------------###

fitAlgo <- function(params, dat, rpt, tp){
  
  ranks <- data.table::data.table(array(
    data = NA, dim = c(params$Ntotal, length(params$rankers))))
  data.table::setnames(ranks, names(params$rankers))
  substs <- vector(mode = "list", length = length(params$subsets))
  names(substs) <- names(params$subsets)
  acc <- list()
  all <- c(params$rankers, params$subsets)
  folds <- caret::createFolds(dat[,target], k = params$kOut, returnTrain = TRUE)
  
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
  
  trControl <- caret::trainControl(method = "cv", number = params$kInn,
                                   selectionFunction = "oneSE")
  
  for (name in intersect(names(all), c("rfRFE", "linSvmRFE", "rbfSvmRFE"))){
    
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
  
  for (name in intersect(names(all), (c("rfSA", "linSvmSA", "rbfSvmSA")))){
    
    ctrlSA <- caret::safsControl(functions = unlist(ifelse(grepl("^rf", name),
                                                           list(rfSA), 
                                                           list(caretSA))),
                                 method = "cv", number = params$kOut, 
                                 index = folds, improve = 10)
    resSA <- caret::safs(x=dat[,.SD,.SDcols = !"target"], y=dat[,target],
                         iters = 10,
                         safsControl = ctrlSA,
                         method = all[[name]],
                         tuneGrid = NULL, # too long otherwise
                         trControl = trControl,
                         differences = F)
    resSA <- dat[]
    substs[[name]] <- resSA$optVariables
    acc[[name]] <- resSA$averages[resSA$optIter, "Accuracy"]
  } 
  
  filt_eval <- FSinR::filterEvaluator('giniIndex')
  hyb_search <- FSinR::whaleOptimization(10, 50)
  res <- hyb_search(dat, "target", featureSetEval = filt_eval)
  substs[["hyb"]] <- dimnames(res$bestFeatures)[[2]][res$bestFeatures==1]

  ids <- rep(NA, params$N)
  for (x in 1:params$kOut) ids[setdiff(1:params$N, folds[[x]])] <- x
  res <- sparseSVM::cv.sparseSVM(X=as.matrix(dat[,.SD,.SDcols = !"target"]), 
                                 y=dat[,target], alpha = 1, gamma = 0.1, nlambda=100,
                                 lambda.min = 0.01, screen = "ASR", max.iter = 1000, eps = 1e-5,
                                 fold.id = ids, nfolds = params$kOut)
  lambda_res <- lambdaSE(res)
  substs[["l1SVM"]] <- 
    names(res$fit$weights[-1,lambda_res$bestIndex])[res$fit$weights[-1,lambda_res$bestIndex] != 0]
  acc[["l1SVM"]] <- 1 - lambda_res$bestME
  
  ranks[, Boruta := (Boruta::Boruta(target~., data = dat))$finalDecision]
  ranks[, ReliefF := FSelectorRcpp::relief(target~., data=dat, 
                                           sampleSize = 10)[,2]]
  substs <- c(substs, setNames(lapply(params$fcbcThres, 
                                      function(x) Biocomb::select.fast.filter(
                                        dat,
                                        disc.method = "MDL",
                                        threshold = x,
                                        attrs.nominal = grep(".*Bin", 
                                                             names(dat)))),
                                paste0("FCBC_", params$fcbcThres)))
  
  saveRDS(list(ranks, substs, acc), paste0("mech_", params$mechanism,
                                           "_pm_", params$pm,
                                           "_corrPred_", params$corrPred,
                                           "_pr_", params$pr,
                                           "_rep_", rpt,
                                           "try2"))
}

###--------------------------------------------------------------------------###

simRep <- function(fixed, varied, rpt, ncond=4){
  
  prev <- rep("none", ncond)
  for (i in 1:24){
    conds <- c(fixed, varied[i,])
    changes <- varied[i, 1:ncond] != prev 
    names(changes) <- colnames(varied)[1:ncond]
    if (changes["pr"] || changes["corrPred"]) dat <- simDat(conds)
    mssng <- simR(conds, dat)
    dat[, target := factor(mssng$R, labels = c("c", "m"))]
    fitAlgo(params = conds, 
            dat = dat[,.SD,.SDcols = !c("Y")], 
            tp = mssng$preds,
            rpt = rpt)
    message(paste0("Condition ", i, "out of ", nrow(varied), 
                   "in repetition ", rpt, "finished!"))
  }
}

###--------------------------------------------------------------------------###

library(doParallel)
library(doRNG)
ncores <- detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

options(error = dump.frames)
a <- Sys.time()
set.seed(1813544)

# make message printed
x2 <- foreach(nmc=1:10, .packages=c("stats",
                                    "mvnfast",
                                    "data.table",
                                    "caret",
                                    "ranger",
                                    "kernlab",
                                    "Boruta",
                                    "FSelectorRcpp",
                                    "FSinR",
                                    "sparseSVM",
                                    "Biocomb")) %dorng% {
                                      simRep(fixedParams, variedParams,nmc)
                                      }
b <- Sys.time()
print(b-a)

stopCluster(cl)
registerDoSEQ()

# results are fully reproducible and seeds are the same

lapply(1:10, function(x) attr(x2, 'rng')[[n]] == attr(x, 'rng')[[n]])

readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1")[[1]] ==
  readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1try2")[[1]]

lapply(1:5, function(n) 
readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1")[[3]][[n]] ==
    readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1try2")[[3]][[n]])

lapply(1:9, function(n) 
  readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1")[[2]][[i]] ==
    readRDS("mech_mcar_pm_0.1_corrPred_0_pr_0.05_rep_1try2")[[2]][[i]])
