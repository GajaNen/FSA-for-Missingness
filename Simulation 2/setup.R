# glossary:
# FSA: feature selection algorithm
# R: non-response indicator (target of FSAs)
# Y: incomplete variable
# R2r: explained variance in the latent variable for R
# PM: percentage of missing values (0.1, 0.5)
# mechanism: missingness mechanism (MCAR, MAR, MNAR)

# dependencies: data.table, sfsmisc

###--------------------------------------------------------------------------###

# each condition: a row of variedParams + all fixedParams

# fixed parameters (don't change over conditions)
# @dir: char, output directory for saving results
# @N: num, number of samples to be taken with replacement from the heart failure data
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
                    seed=1813544,
                    streams=500,
                    kInn=3, 
                    kOut=5, 
                    sizes = c(0, 1, 2, 5),
                    tuneGrids = list(lasso=expand.grid(alpha=1,
                                                       lambda=lseq(0.004,4,40)),
                                     EN=expand.grid(alpha=lseq(0.3, 0.9, 2),
                                                    lambda=lseq(0.004,4,30)),
                                     ranger=expand.grid(min.node.size=c(1, 3),
                                                        mtry=c(floor(sqrt(10)), 4),
                                                        splitrule="gini",
                                                        importance="impurity"),
                                     svmRadial=expand.grid(C=lseq(2**(-5),2**(15),10), # was not tuned due to computational limitations
                                                           sigma=lseq(2**(-5),2**(15),5))), 
                    deg = 3,
                    Ntotal = 10,
                    fcbcThres = list(0, 0.1, 0.3, 0.5),
                    rankers = list("LR"="glm",  
                                   "ReliefF"=NULL,
                                   "Boruta"=NULL,
                                   "RF"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet", 
                                   "rbfSvmSA"="svmRadial", # was not tuned due to computational limitations
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL,
                                   "hyb"=NULL)
)

###--------------------------------------------------------------------------###

# varied factors (change over conditions)

variedParams <- data.table::setDT(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       R2r = c(0.3, 0.6)))
)

###--------------------------------------------------------------------------###

# don't vary R2r for MCAR
variedParams <- variedParams[!(mechanism=="mcar" & R2r == 0.6)]
variedParams <- variedParams[mechanism=="mcar", R2r := NA]

###--------------------------------------------------------------------------###
