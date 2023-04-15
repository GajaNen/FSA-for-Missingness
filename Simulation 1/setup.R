
# each condition: a row of variedParams + all fixedParams

# fixed parameters
# FOR TESTING PURPOSES ONLY
# @dir: char, output directory for saving results
# @N: num, number of samples
# @addY: logical (T or F, not NULL!) simulate Y with a copula or not (specifically with X relevant)
# @kInn: int, number of inner level folds for algorithms which use two levels of resampling scheme
# @kOut: int. number of folds for outer (or only) level of resampling
# @sizes: vector of ints, sizes considered in RFE
# @tuneGrids: a list of named vectors, tuning grids for all functions from caret package
# @Ntotal: int, N of all variables, must be a multiple of 3
# @map.funcs: names list, mapping from name of parameters to quantile function to be applied for Gaussian copula
# @deg: int, degree of b-splines
# @fcbcThres: vector or list of nums, thresholds for FCBC (fast correlation based filter)
# @rankers: list of names:function (RFE, SA, glmnet, i.e. all from caret) or names of ranker FSA (others)
# @subsest: list of name:NULL for methods which output subsets
fixedParams <- list(dir=file.path("Simulation 1", "Results"),
                    seed=1813544,
                    streams=500,
                    N=1000, 
                    addY=T,
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

# always generate X relevant and Y together, irrelevant separately
# TWO LEVELS: COMPLETELY INDEPEDENT RELEVANT VARIABLES WITH SOME DEPENDENCY WITH Y
# AND SOME DEPENDECY BETWEEN RELEVANT VARIABLES AND WITH Y
# IRRELEVANT VARIABLES ARE ALWAYS INDEPENDENT FROM EACH OTHER AND X RELEVANT AND Y

# THE SAID DEPENDENCY IS NON-PARAMETRIC (E.G. KENDALL'S CORRELATION)
# WHICH FOR GAUSSIAN COPULA HAS A BIJECTIVE TRANSFORMATION TO PEARSON'S RHO
# this associated rho matrix then IS
# USED IN GENERATING NORMALS FOR COPULA --- KENDALL'S CORRELATION PRESERVED 
# WITH TRANSFORMATIONS TO OTHER DISTRIBUTIONS

# IN RELATION TO THE picture in the thesis draft: X->Y, Y->Xrel & no Z &
# additional set of independent (between each other & with Xrel and Y) irrelevant variables

# population correlation matrices for Gaussian copula 

# THIS IS A PROVISIONAL WAY TO GENERATE COR MATRIX JUST FOR TESTING PURPOSES
# in the actual population matrix used for the simulation, there will be dependency
# Between Xrel and Xrel-Y or only Xrel-Y


#cor mat low percentage relevant variables and Y when corrPred==1
n <- fixedParams$addY==T
cm_lp <- cov2cor(genPositiveDefMat(dim=fixedParams$Ntotal*0.05+n, covMethod = "eigen")$Sigma)
#cor mat high percentage relevant variables (Xrel and Y) when corrPred==1
cm_hp <- cov2cor(genPositiveDefMat(dim=fixedParams$Ntotal*0.2+n, covMethod = "eigen")$Sigma)


Nt <- fixedParams$Ntotal

# THIS IS A PROVISIONAL WAY TO GENERATE CORS XREL-Y JUST
# TO MAKE SURE THERE'S SOME DEPENDENCY
R2y <- 0.15

# matrices (high and low percentage of rel vars) for corrPred==0

lp_diag <- diag(1,Nt*0.05)
cors <- runif(Nt*0.05, 0.4, 0.6)
a <- sqrt(R2y/(t(cors) %*% solve(lp_diag) %*% cors))[1] 
cors_lp_diag <- a*cors
# cors_lp_diag %*% solve(lp_diag) %*% cors_lp_diag

hp_diag <- diag(1,Nt*0.2)
cors <- runif(Nt*0.2, 0.4, 0.6)
a <- sqrt(R2y/(t(cors) %*% solve(hp_diag) %*% cors))[1] 
cors_hp_diag <- a*cors
#cors_hp_diag %*% solve(hp_diag) %*% cors_hp_diag

cm_lp_diag <- cbind(rbind(lp_diag, cors_lp_diag), c(cors_lp_diag,1))
cm_hp_diag <- cbind(rbind(hp_diag, cors_hp_diag),c(cors_hp_diag,1))


# varied factors
# PM: MISSINGNESS PERCENTAGE, PR: PERCENTAGE OF RELEVANT VARS
# CORRPRED: WHETHER RELEVANT VARIABLES HAVE SOME DEPENDENCY OR NOT
variedParams <- data.table::setDT(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       corrPred = c(0, 1),
       pr = c(0.05, 0.2)))
)


# ALL THESE PARAMETERS ARE JUST FOR TESTING PURPOSES
# PR WON'T BE VARIED IN CASES WHEN MECHANISM IS MCAR TO REDUCE RUN TIME
# SINCE INDICATOR IS GENERATED RANDOMLY ANYWAY
# another way to reduce number of conditions is to use only Y as a predictor in mnar
# this way at least PR doesn't have to be varied (and maybe even corrPred), which means
# at least 4 conditions less 

variedParams <- variedParams[!(mechanism=="mcar" & pr == 0.2)]

# TRANSFORMATIONS FOR DISCRETE VARS HIGH AND LOW PR VERSIONS
HPtrans <- c(sin, exp, cos, sin, sin, sin, sin, sin)
LPtrans <- c(sin, exp)

variedParams[, trans := fifelse(pr==0.2, list(HPtrans), list(LPtrans))]
# fix pseudo R^2 in missingness indicators (maybe per mechanism)
variedParams[mechanism!="mcar", R2r := 0.6]

# always 1/3 Cont, Bin, Ord in this order (but within each type, distributions can var)
# set some test parameters for normal, gamma, binomial, logistic distributions (relevant)
# and population baseline probs for ordinal vars + Y
Nt <- fixedParams$Ntotal
norms <- list(c(mean = 0,sd = 1))
gams <- list(c(shape = 0.2, rate = 1))
binsp <- list(c(size=1, prob=0.3))
logsp <- list(c(location = 0, scale = 1))
paramy <- list(c(shape1=2,shape2=5))
pops <- c(0.2, 0.1, 0.3, 0.2, 0.2)

variedParams[, paramsRel := lapply(pr, function(n) c(rep(norms, (Nt)*n/6),
                                                     rep(gams, (Nt)*n/6),
                                                     rep(binsp, (Nt)*n/3),
                                                     rep(logsp, (Nt)*n/3),
                                                     paramy))]

variedParams[, popProbsRel := 
               lapply(pr, function(n) lapply(1:((Nt)*n/3),
                                             function(x) pops))]

# set some params for normal, gamma, binomial distributions (irrelevant)
variedParams[, paramsIrrel := lapply(pr, function(n) c(rep(norms, (Nt) * (1 - n) / 6),
                                                       rep(gams, (Nt) * (1 - n) / 6),
                                                       rep(binsp, (Nt) * (1 - n) / 3),
                                                       rep(logsp, (Nt) * (1 - n) / 3)))]

variedParams[, popProbsIrrel := 
               lapply(pr,function(n) lapply(1:((Nt) * (1 - n) / 3),
                                            function(x) pops))]

# set repetitions of each type of var
variedParams[, repsRel := lapply(pr, function(x) c(Nt*x/6, Nt*x/6, Nt*x/3, Nt*x/3, 1))]
variedParams[, repsIrrel := lapply(pr, function(x) c(Nt*(1-x)/6, Nt*(1-x)/6, Nt*(1-x)/3, Nt*(1-x)/3))]
# set distributions to be used for mapping to appropriate quantile function
variedParams[, distsRel := lapply(repsRel, 
                                  function(x) unlist(mapply(rep, names(fixedParams$map.funcs), x,
                                                            USE.NAMES = F)))]

variedParams[, distsIrrel := lapply(repsIrrel, 
                                    function(x) unlist(mapply(rep, names(fixedParams$map.funcs)[seq_along(x)], x,
                                                              USE.NAMES = F)))]

# define correlation matrix which differs between percentage of relevant variables
# for the correlated predictors case and uncorrelated case

variedParams[pr==0.05 & corrPred, corMatRel := list(cm_lp)]
variedParams[pr!=0.05 & corrPred, corMatRel := list(cm_hp)]

variedParams[pr==0.05 & corrPred==0, corMatRel := list(cm_lp_diag)]
variedParams[pr!=0.05 & corrPred==0, corMatRel := list(cm_hp_diag)]

# for irrel variables, generate diagonal matrices
variedParams[, corMatIrrel := lapply(pr, function(x) diag(1, Nt*(1-x)))]


# coefs for splines
variedParams[, theta := list(c(0.3, 0.5, 0.6, 0.1, 0.2))]
