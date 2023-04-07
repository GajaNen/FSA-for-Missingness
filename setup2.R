options(scipen=999)
lapply(list("Data/simData.R", "Data/simMiss.R", "helpers.R", "depend.R"), source)
source("Data/modRfFuncs.R")

# fixed parameters
fixedParams <- list(N=1000, 
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
                    paramY = c(shape1=2,shape2=5),
                    margY = "beta",
                    R2y = 0.15,
                    deg = 3,
                    fcbcThres = list(0.2, 0.3, 0.4),
                    rankers = list("LR"="glm",  
                                   "ReliefF"=NULL,
                                   "Boruta"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet",
                                   "rbfSvmRFE"="svmRadial", 
                                   "rbfSvmSA"="svmRadial",
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL,
                                   "hyb"=NULL)
)


# population correlation matrices for Gaussian copula
library(clusterGeneration)
cm_lp <- genPositiveDefMat(dim=fixedParams$Ntotal*0.05, covMethod = "eigen")$Sigma
cm_lp<-cov2cor(cm_lp)
cors <- runif(fixedParams$Ntotal*0.05, 0.4, 0.6)
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                             solve(cm_lp) %*% 
                             cors))[1] # theoretical bounds for cor?
cors_lp <- a*cors
cm_lp <- cbind(rbind(cm_lp, cors_lp), c(cors_lp,1))
#cors_lp %*% solve(cm_lp) %*% cors_lp
cm_hp <- genPositiveDefMat(dim=fixedParams$Ntotal*0.2, covMethod = "eigen")$Sigma
cm_hp <- cov2cor(cm_hp)
cors <- runif(fixedParams$Ntotal*0.2, 0.4, 0.6)
a <- sqrt(fixedParams$R2y/(t(cors) %*% 
                                 solve(cm_hp) %*% 
                                 cors))[1] # theoretical bounds for cor?
cors_hp <- a*cors
cm_hp <- cbind(rbind(cm_hp, cors_hp), c(cors_hp,1))
#cors_hp %*% solve(cm_hp) %*% cors_hp

Nt <- fixedParams$Ntotal

cm_lp_diag <- unlist(sapply(1:(Nt*0.05), function(x) c(rep(0, Nt*0.05-x), cors_lp[x])))
cm_hp_diag <- unlist(sapply(1:(Nt*0.2), function(x) c(rep(0, Nt*0.2-x), cors_hp[x])))


# varied factors
variedParams <- data.table::as.data.table(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       corrPred = c(0, 1),
       pr = c(0.05, 0.2)))
)


# don't vary pr for mcar remove with datatable magic
#variedParams[(mechanism=="mcar" & pr == 0.2) := NULL,]

HPtrans <- c(sin, exp, cos, sin, sin, sin, sin, sin,
             sin, exp, exp, exp, sin, exp, exp, exp)
LPtrans <- c(sin, exp, cos, sin)

variedParams[, trans := fifelse(pr==0.2, list(HPtrans), list(LPtrans))]
# fix pseudo R^2 in missingness indicators per mechanism
variedParams[mechanism!="mcar", R2r := 0.6]
# fix R^2 in incomplete variable (maybe per mechanism)
#variedParams[, R2y := 0.3]

# set some test parameters for normal, gamma, binomial distributions (relevant)
Nt <- fixedParams$Ntotal
norms <- list(c(mean = 0,sd = 1))
variedParams[, normParamRel := 
               lapply(pr,function(n) rep(norms, (Nt)*n/6))]

gams <- list(c(shape = 0.2, rate = 1))
variedParams[, gamParamRel := 
               lapply(pr,function(n) rep(gams, (Nt)*n/6))]

binsp <- list(c(size=1, prob=0.3))
variedParams[, binParamRel := 
               lapply(pr,function(n) rep(binsp, (Nt)*n/3))]

logsp <- list(c(location = 0, scale = 1))
variedParams[, logParamRel := 
               lapply(pr,function(n) rep(logsp, (Nt)*n/3))]

pops <- c(0.2, 0.1, 0.3, 0.2, 0.2)
variedParams[, popProbsRel := 
               lapply(pr,function(n) matrix(rep(pops, (Nt)*n/3),
                                            nrow=(Nt)*n/3,
                                            byrow=T))]

# set some toy parameter for normal, gamma, binomial distributions (irrelevant)
variedParams[, normParamIrrel := 
               lapply(pr,function(n) rep(norms, (Nt) * (1 - n) / 6))]

variedParams[, gamParamIrrel := 
               lapply(pr,function(n) rep(gams, (Nt) * (1 - n) / 6))]

variedParams[, binParamIrrel := 
               lapply(pr,function(n) rep(binsp, (Nt) * (1 - n) / 3))]

variedParams[, logParamIrrel := 
               lapply(pr,function(n) rep(logsp, (Nt) * (1 - n) / 3))]


variedParams[, popProbsIrrel := 
               lapply(pr,function(n) matrix(rep(pops, (Nt) * (1 - n) / 3),
                                            nrow=(Nt) * (1 - n) / 3,
                                            byrow=T))]

# define correlation matrix which differs between percentage of relevant variables
# for the correlated predictors case and uncorrelated case

variedParams[pr==0.05 & corrPred, corMatRel := list(cm_lp[lower.tri(cm_lp)])]
variedParams[pr!=0.05 & corrPred, corMatRel := list(cm_hp[lower.tri(cm_hp)])]

variedParams[pr==0.05 & corrPred==0, corMatRel := list(cm_lp_diag)]
variedParams[pr!=0.05 & corrPred==0, corMatRel := list(cm_hp_diag)]

variedParams[, corMatIrrel := lapply(pr, function(x) rep(0, 
                                                         (Nt*(1-x))*(Nt*(1-x)-1)/2))]

variedParams[, theta := list(c(0.3, 0.5, 0.6, 0.1, 0.2))]

params <- c(fixedParams, variedParams[1,])


a <- Sys.time()
# check if it works
for (i in 1:24){
  parms <- c(fixedParams, variedParams[i,])
  dat <- simDat(parms)
}
b <- Sys.time()
print(b-a)
