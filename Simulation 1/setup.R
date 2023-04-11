# fixed parameters
# FOR TESTING PURPOSES ONLY
# @dir: output directory
# @N: samples
# @addY: simulate Y with a copula or not (specifically with X relevant)
# @kInn: number of folds for algorithms which use two levels of resampling scheme
# @kOut: number of folds for outer (or only) level of resampling
# @sizes: sizes considered in RFE
# @tuneGrids: tuning grids for all functions from caret package
# @Ntotal: N of all variables
# @map.funcs: mapping from name of parameters to function to be applied
# @R2y: used to generate some depedency between Xrel and Y (as a population matrix)
# @deg: degree of b-splines
# @fcbcThres: thresholds for FCBC (fast correlation based filter)
# @rankers: list of names:function (RFE, SA, glmnet, i.e. all from caret) or names of ranker FSA (others)
# @subsest: same but for subsets
fixedParams <- list(dir="Results",
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
                    R2y = 0.15,
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


# THIS IS A PROVISIONAL WAY TO GENERATE COR MATRIX JUST FOR TESTING PURPOSES
# always generate X relevant and Y together, irrelevant separately
# TWO LEVELS: COMPLETELY INDEPEDENT RELEVANT VARIABLES WITH SOME DEPENDENCY WITH Y
# AND SOME DEPENDECY BETWEEN RELEVANT VARIABLES AND WITH Y
# THE SAID DEPENDENCY IS NON-PARAMETRIC (E.G. KENDALL'S CORRELATION)
# WHICH FOR GAUSSIAN COPULA HAS A BIJECTIVE TRANSFORMATION TO PEARSON'S RHO
# this associated rho matrix then 
# USED IN GENERATING NORMALS FOR COPULA --- KENDALL'S CORRELATION PRESERVED 
# WITH TRANSFORMATIONS TO OTHER DISTRIBUTIONS
# IRRELEVANT VARIABLES ARE ALWAYS COMPLETELY INDEPEDENT
# so basically it's the picture in the thesis draft, with X->Y, Y->Xrel & no Z &
# additional set of indepedent (between & with Xrel and Y) irrelevant variables

# population correlation matrices for Gaussian copula 

#cor mat low percentage relevant variables (Xrel and Y)
cm_lp <- genPositiveDefMat(dim=fixedParams$Ntotal*0.05, covMethod = "eigen")$Sigma
cm_lp <- cov2cor(cm_lp)
cors <- runif(fixedParams$Ntotal*0.05, 0.4, 0.6)
# some depedency with Y
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                             solve(cm_lp) %*% 
                             cors))[1] # theoretical bounds for cor?
cors_lp <- a*cors
cm_lp <- cbind(rbind(cm_lp, cors_lp), c(cors_lp,1))
#cors_lp %*% solve(cm_lp) %*% cors_lp

#cor mat high percentage relevant variables (Xrel and Y)
cm_hp <- genPositiveDefMat(dim=fixedParams$Ntotal*0.2, covMethod = "eigen")$Sigma
cm_hp <- cov2cor(cm_hp)
cors <- runif(fixedParams$Ntotal*0.2, 0.4, 0.6)
#some depedency between Xrel and Y
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                                 solve(cm_hp) %*% 
                                 cors))[1] # theoretical bounds for cor?
cors_hp <- a*cors
cm_hp <- cbind(rbind(cm_hp, cors_hp), c(cors_hp,1))
#cors_hp %*% solve(cm_hp) %*% cors_hp

Nt <- fixedParams$Ntotal

# matrices (high and low percentage of rel vars) for corrPred==0
# with the same dependency with Y (linear, which is not preserved with a copula, but just
# to create some form of depenndce structure)

lp_diag <- diag(1,Nt*0.05)
cors <- runif(fixedParams$Ntotal*0.05, 0.4, 0.6)
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                             solve(lp_diag) %*% 
                             cors))[1] # theoretical bounds for cor?
cors_lp_diag <- a*cors

hp_diag <- diag(1,Nt*0.2)

cors <- runif(fixedParams$Ntotal*0.2, 0.4, 0.6)
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                             solve(hp_diag) %*% 
                             cors))[1] # theoretical bounds for cor?
cors_hp_diag <- a*cors

cm_lp_diag <- cbind(rbind(lp_diag, cors_lp_diag), c(cors_lp_diag,1))
cm_hp_diag <- cbind(rbind(hp_diag, cors_hp_diag),c(cors_hp_diag,1))


# varied factors
# PM: MISSINGNESS PERCENTAGE, PR: PERCENTAGE OF RELEVANT VARS
# CORRPRED: WHETHER RELEVANT VARIABLES HAVE SOME DEPENDENCY OR NOT
variedParams <- data.table::as.data.table(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       corrPred = c(0, 1),
       pr = c(0.05, 0.2)))
)


# ALL THESE PARAMETERS ARE JUST FOR TESTING PURPOSES
# PR WON'T BE VARIED IN CASES WHEN MECHANISM IS MCAR TO REDUCE RUN TIME
# SINCE INDICATOR IS GENERATED RANDOMLY ANYWAY

variedParams <- variedParams[!(mechanism=="mcar" & pr == 0.2)]

# TRANSFORMATIONS FOR DISCRETE VARS HIGH AND LOW PR VERSIONS
HPtrans <- c(sin, exp, cos, sin, sin, sin, sin, sin,
             sin, exp, exp, exp, sin, exp, exp, exp)
LPtrans <- c(sin, exp, cos, sin)

variedParams[, trans := fifelse(pr==0.2, list(HPtrans), list(LPtrans))]
# fix pseudo R^2 in missingness indicators (maybe per mechanism)
variedParams[mechanism!="mcar", R2r := 0.6]

# set some test parameters for normal, gamma, binomial, logistic distributions (relevant)
# and population baseline probs for ordinal vars
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
