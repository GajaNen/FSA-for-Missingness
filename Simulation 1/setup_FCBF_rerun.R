# glossary:
# FSA: feature selection algorithm
# Xrel: relevant variables (i.e. true predictors)
# Xirrel: irrelevant variables (i.e. non-predictors)
# R: non-response indicator (target of FSAs)
# Y: incomplete variable
# PR: percentage of relevant variables (0.05, 0.2)
# PC: predictor correlation (Yes, No) - the corresponding variable is corrPred
# PM: percentage of missing values (0.1, 0.5)
# mechanism: missingness mechanism (MCAR, MAR, MNAR)


# dependencies: clusterGeneration, data.table, sfsmisc

###--------------------------------------------------------------------------###

# each condition: a row of variedParams + all fixedParams

# fixed parameters (don't change over conditions)
# @dir: char, output directory for saving results
# @N: num, number of samples
# @addY: logical (T or F, not NULL!) simulate Y with a copula or not (specifically with Xrel or not)
# @kInn: int, number of inner level folds for algorithms which use two levels of resampling scheme
# @kOut: int. number of folds for outer (or only) level of resampling
# @sizes: vector of ints, sizes considered in RFE
# @tuneGrids: a list of named vectors, tuning grids for all functions from caret package
# @Ntotal: int, N of all variables, must be a multiple of 3
# @map.funcs: named list, mapping from the names of parameters to quantile function to be applied for Gaussian copula
# @deg: int, degree of b-splines to be applied to continuous relevant variables in generating R with a GAM
# @fcbcThres: vector or list of nums, thresholds for FCBF (fast correlation based filter)
# @rankers: list of names:function (RFE, SA, glmnet, i.e. all from caret) or names of ranker FSA (others)
# @subsets: list of name:NULL for methods which output subsets
fixedParams <- list(simN=1,
                    dir=file.path(getwd(), "Simulation 1", "rerunFCBF_sim1"),
                    seed=1813544,
                    streams=500,
                    N=1000, 
                    addY=T,
                    kInn=3, 
                    kOut=5, 
                    sizes = c(0, 8, 24, 60),
                    tuneGrids = list(lasso=expand.grid(alpha=1,
                                                       lambda=lseq(0.004,4,40)),
                                     EN=expand.grid(alpha=lseq(0.3, 0.9, 2),
                                                    lambda=lseq(0.004,4,30)),
                                     ranger=expand.grid(min.node.size=c(1, 10),
                                                        mtry=c(floor(sqrt(120)), 24),
                                                        splitrule="gini",
                                                        importance="impurity")),
                    Ntotal = 120,
                    map.funcs = list("normParam"=qnorm, "gamParam"=qgamma,"binParam"=qbinom,
                                     "logParam"=qlogis, "paramY"=qbeta),
                    deg = 3,
                    fcbcThres = list(0, 0.1, 0.3, 0.5),
                    rankers = list("LR"="glm",  
                                   "ReliefF"=NULL,
                                   "Boruta"=NULL,
                                   "RF"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet", 
                                   "rbfSvmSA"="svmRadial",
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL)
)

###--------------------------------------------------------------------------###

# population correlation matrices for Gaussian copula 

n <- fixedParams$addY==T
Nt <- fixedParams$Ntotal

# correlation matrix for PR=0.05 of Xrel and Y when PC=Yes
cm_lp <- rcorrmatrix(Nt*0.05 + n, 1)

# correlation matrix for PR=0.2 of Xrel and Y when PC=Yes
cm_hp <- rcorrmatrix(Nt*0.2 + n, 1)

# correlation matrices (high and low PR) when PC=No

R2y <- 0.2

## PR = 0.05
lp_diag <- diag(1,Nt*0.05) # correlations between Xrel
# arbitrarily generate pearson's rho correlations Xrel~Y 
# such that explained variance in Y is 0.2
# just to obtain some non-zero Kendall's correlations and a valid
# correlation matrix
cors <- runif(Nt*0.05, 0.4, 0.6)
a <- sqrt(R2y/(t(cors) %*% solve(lp_diag) %*% cors))[1] 
cors_lp_diag <- a*cors

## PR = 0.2
hp_diag <- diag(1,Nt*0.2) # correlations beteween Xrel
# arbitrarily generate pearson's rho correlations Xrel~Y 
# such that explained variance in Y is 0.2
# just to obtain some non-zero Kendall's correlations and a valid
# correlation matrix
cors <- runif(Nt*0.2, 0.4, 0.6)
a <- sqrt(R2y/(t(cors) %*% solve(hp_diag) %*% cors))[1] 
cors_hp_diag <- a*cors

cm_lp_diag <- cbind(rbind(lp_diag, cors_lp_diag), c(cors_lp_diag,1))
cm_hp_diag <- cbind(rbind(hp_diag, cors_hp_diag),c(cors_hp_diag,1))

###--------------------------------------------------------------------------###


# varied factors (change over conditions)
variedParams <- data.table::setDT(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       corrPred = c(0, 1),
       pr = c(0.05, 0.2)))
)

# don't vary PR when mechanism=MCAR
variedParams <- variedParams[!(mechanism=="mcar" & pr == 0.2)]

###--------------------------------------------------------------------------###

# transformations for ordinal Xrel for generating R under MAR and MNAR with a GAM

HPtrans <- c(sin, exp, exp, sin, log, tan, exp, cospi) # PR=0.2
LPtrans <- c(sin, exp) # PR=0.05

variedParams[, trans := fifelse(pr==0.2, list(HPtrans), list(LPtrans))]

###--------------------------------------------------------------------------###

# fix pseudo R^2 in missingness indicators
variedParams[mechanism!="mcar", R2r := 0.6]

###--------------------------------------------------------------------------###

# marginals for Xrel and Xirrel
## always 1/3 Continuous, Binary, Ordinal in this order 
## (but within each type, distributions can vary)

logsp <- list(c(location = 0, scale = 1)) # marginal for ordinal variables (always the same)

## marginals for Xrel

paramy <- list(c(shape1=2,shape2=5)) # marginal for Y

### PR = 0.2
norms_hp_r <- c(rep(list(c(mean = 0, sd = 1)), 2), rep(list(c(mean = 2, sd = 5)), 2))
gams_hp_r <- c(rep(list(c(shape = 0.2, rate = 1)), 2), rep(list(c(shape = 2, rate=3)), 2))
binsp_hp_r <- c(rep(list(c(size=1, prob=0.3)), 2), rep(list(c(size=1, prob=0.7)), 2),
                rep(list(c(size=1, prob=0.5)), 2), rep(list(c(size=1, prob=0.1)), 2))

### PR = 0.05
norms_lp_r  <- list(c(mean = 0, sd = 1))
gams_lp_r <- list(c(shape = 0.2, rate = 1))
binsp_lp_r <- c(rep(list(c(size=1, prob=0.3)), 1), rep(list(c(size=1, prob=0.5)), 1))


variedParams[, paramsRel := 
               lapply(pr, function(x){
                 if(x == 0.2){
                   c(norms_hp_r, gams_hp_r, binsp_hp_r, rep(logsp, (Nt)*x/3),paramy)
                 } else {
                   c(norms_lp_r, gams_lp_r, binsp_lp_r, rep(logsp, (Nt)*x/3),paramy)
                 }
               })
]


### population cumulative probabilities of each category for ordinal Xrel
variedParams[, popProbsRel := 
               lapply(pr, function(x)
                 c(rep(list(c(0.2, 0.1, 0.3, 0.2, 0.2)), Nt*x/6), 
                   rep(list(c(0.4, 0.1, 0.1, 0.2, 0.2)), Nt*x/6))
               )
]

## marginals for Xirrel
Ni_lp <- Nt-Nt*0.05
Ni_hp <- Nt-Nt*0.2

### PR = 0.2
norms_hp_i <- c(rep(list(c(mean = 0, sd = 1)), Ni_hp/12), rep(list(c(mean = 2, sd = 5)), Ni_hp/12))
gams_hp_i <- c(rep(list(c(shape = 0.2, rate = 1)), Ni_hp/12), rep(list(c(shape = 2, rate=3)), Ni_hp/12))
binsp_hp_i <- c(rep(list(c(size=1, prob=0.3)), Ni_hp/12), rep(list(c(size=1, prob=0.7)), Ni_hp/12),
                rep(list(c(size=1, prob=0.5)), Ni_hp/12), rep(list(c(size=1, prob=0.1)), Ni_hp/12))

### PR = 0.05
norms_lp_i <- rep(list(c(mean = 0, sd = 1)), Ni_lp/6)
gams_lp_i <- rep(list(c(shape = 0.2, rate = 1)), Ni_lp/6)
binsp_lp_i <- c(rep(list(c(size=1, prob=0.3)), Ni_lp/6), 
                rep(list(c(size=1, prob=0.5)), Ni_lp/6))


variedParams[, paramsIrrel := 
               lapply(pr, function(x){
                 if (x == 0.05){
                   c(norms_lp_i, gams_lp_i, binsp_lp_i, rep(logsp, (Nt)*(1-x)/3))
                 } else {
                   c(norms_hp_i, gams_hp_i, binsp_hp_i, rep(logsp, (Nt)*(1-x)/3))
                 }
               })]

variedParams[, popProbsIrrel := 
               lapply(pr, function(x)
                 c(rep(list(c(0.2, 0.1, 0.3, 0.2, 0.2)), Nt*(1-x)/6), 
                   rep(list(c(0.4, 0.1, 0.1, 0.2, 0.2)), Nt*(1-x)/6))
               )
]

## set repetitions of each type of marginal distribution (e.g. normal, gamma, etc.)
variedParams[, repsRel := 
               lapply(pr, function(x) c(Nt*x/6, Nt*x/6, Nt*x/3, Nt*x/3, 1))]
variedParams[, repsIrrel := 
               lapply(pr, function(x) c(Nt*(1-x)/6, Nt*(1-x)/6, Nt*(1-x)/3, Nt*(1-x)/3))]


## set names of desired distributions to be used for mapping to appropriate quantile function
variedParams[, distsRel := 
               lapply(repsRel,
                      function(x) unlist(mapply(rep, names(fixedParams$map.funcs), x,
                                                USE.NAMES = F)))]

variedParams[, distsIrrel := 
               lapply(repsIrrel,
                      function(x) unlist(mapply(rep, names(fixedParams$map.funcs)[seq_along(x)], x,
                                                USE.NAMES = F)))]

###--------------------------------------------------------------------------###

# define correlation matrix of Xrel and Y, which differ across levels of PR

## for PC=Yes
variedParams[pr==0.05 & corrPred, corMatRel := list(cm_lp)]
variedParams[pr!=0.05 & corrPred, corMatRel := list(cm_hp)]

## for PC=No
variedParams[pr==0.05 & corrPred==0, corMatRel := list(cm_lp_diag)]
variedParams[pr!=0.05 & corrPred==0, corMatRel := list(cm_hp_diag)]

# generate diagonal matrices for Xirrel
variedParams[, corMatIrrel := lapply(pr, function(x) diag(1, Nt*(1-x)))]

###--------------------------------------------------------------------------###

# coefficient for basis function in splines (generating R under MAR and MNAR)
variedParams[, theta := list(runif(fixedParams$deg + 2, 0.1, 0.9))]
