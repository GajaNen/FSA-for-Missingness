# each condition: a row of variedParams + all fixedParams

# fixed parameters
# FOR TESTING PURPOSES ONLY
# @dir: char, output directory for saving results
# @N: num, number of samples
# @kInn: int, number of inner level folds for algorithms which use two levels of resampling scheme
# @kOut: int. number of folds for outer (or only) level of resampling
# @sizes: vector of ints, sizes considered in RFE
# @tuneGrids: a list of named vectors, tuning grids for all functions from caret package
# @deg: int, degree of b-splines
# @fcbcThres: vector or list of nums, thresholds for FCBC (fast correlation based filter)
# @rankers: list of names:function (RFE, SA, glmnet, i.e. all from caret) or names of ranker FSA (others)
# @subsest: list of name:NULL for methods which output subsets
fixedParams <- list(dir=file.path("Simulation 2", "Results"),
                    N=299, 
                    kInn=3, 
                    kOut=5, 
                    sizes = c(0:6),
                    tuneGrids = list(lasso=expand.grid(alpha=1,
                                                       lambda=seq(0.1,4,0.1)),
                                     EN=expand.grid(alpha=seq(0.1, 0.99, 0.1),
                                                    lambda=seq(0.1,4,0.5)),
                                     ranger=expand.grid(min.node.size=8:10,
                                                        mtry=1:3,
                                                        splitrule="gini",
                                                        importance="impurity"),
                                     svmLinear=expand.grid(C=seq(1,10,1)),
                                     svmRadial=expand.grid(C=seq(1,10,1),
                                                           sigma=seq(0.1,1,0.2))),
                    Ntotal = 120,
                    map.funcs = list("normParam"=qnorm, "gamParam"=qgamma,"binParam"=qbinom,
                                     "logParam"=qlogis, "paramY"=qbeta),
                    deg = 3,
                    fcbcThres = list(0, 0.2, 0.3),
                    rankers = list("LR"="glm",  
                                   "ReliefF"=NULL,
                                   "Boruta"=NULL,
                                   "RF"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet", 
                                   "rbfSvmSA"="svmRadial",
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL,
                                   "hyb"=NULL)
)


# varied factors

variedParams <- data.table::setDT(expand.grid(
  list(mechanism = c("mcar", "marr", "marc", "mart", "mnar"),
       pm = c(0.1, 0.5),
       R2r = c(0.3, 0.6)))
)


variedParams <- variedParams[!(mechanism=="mcar" & R2r == 0.6)]
variedParams <- variedParams[mechanism=="mcar", R2r := NA]

# coefs for splines
variedParams[, theta := list(c(0.3, 0.5, 0.6, 0.1, 0.2))]
