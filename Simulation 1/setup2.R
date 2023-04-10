options(scipen=999)
lapply(list("Data/simData.R", "Data/simMiss.R", "helpers.R", "depend.R"), source)
source("Data/modRfFuncs.R")

# fixed parameters
fixedParams <- list(N=1000, 
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
                    fcbcThres = list(0.2, 0.3, 0.4),
                    rankers = list("LR"="glm",  
                                   "ReliefF"=NULL,
                                   "Boruta"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet", 
                                   "rbfSvmSA"="svmRadial",
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL,
                                   "hyb"=NULL)
)


# population correlation matrices for Gaussian copula
library(clusterGeneration)
cm_lp <- genPositiveDefMat(dim=fixedParams$Ntotal*0.05, covMethod = "eigen")$Sigma
cm_lp <- cov2cor(cm_lp)
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

#cm_lp_diag <- unlist(sapply(1:(Nt*0.05), function(x) c(rep(0, Nt*0.05-x), cors_lp[x])))
#cm_hp_diag <- unlist(sapply(1:(Nt*0.2), function(x) c(rep(0, Nt*0.2-x), cors_hp[x])))

cm_lp_diag <- cbind(rbind(diag(1,Nt*0.05), cors_lp), c(cors_lp,1))
cm_hp_diag <- cbind(rbind(diag(1,Nt*0.2), cors_hp),c(cors_hp,1))


# varied factors
variedParams <- data.table::as.data.table(expand.grid(
  list(mechanism = c("mcar", "mar", "mnar"),
       pm = c(0.1, 0.5),
       corrPred = c(0, 1),
       pr = c(0.05, 0.2)))
)


# don't vary pr for mcar remove with datatable magic
variedParams <- variedParams[!(mechanism=="mcar" & pr == 0.2)]

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

# set some toy parameter for normal, gamma, binomial distributions (irrelevant)
variedParams[, paramsIrrel := lapply(pr, function(n) c(rep(norms, (Nt) * (1 - n) / 6),
                                                      rep(gams, (Nt) * (1 - n) / 6),
                                                      rep(binsp, (Nt) * (1 - n) / 3),
                                                      rep(logsp, (Nt) * (1 - n) / 3)))]

variedParams[, popProbsIrrel := 
               lapply(pr,function(n) lapply(1:((Nt) * (1 - n) / 3),
                                            function(x) pops))]


variedParams[, repsRel := lapply(pr, function(x) c(Nt*x/6, Nt*x/6, Nt*x/3, Nt*x/3, 1))]
variedParams[, repsIrrel := lapply(pr, function(x) c(Nt*(1-x)/6, Nt*(1-x)/6, Nt*(1-x)/3, Nt*(1-x)/3))]

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

variedParams[, corMatIrrel := lapply(pr, function(x) diag(1, Nt*(1-x)))]

variedParams[, theta := list(c(0.3, 0.5, 0.6, 0.1, 0.2))]

#params <- c(fixedParams, variedParams[1,])


a <- Sys.time()
#prevdat <- simDat(c(fixedParams, variedParams[i,]))
# check if it works
for (i in 1:nrow(variedParams)){
  parms <- c(fixedParams, variedParams[i,])
  dat <- simDat(parms)
  #print(identical(prevdat, newdat))
  #prevdat <- newdat
  res <- simR(parms, dat)
  dat[, target := factor(res$R, labels = c("c", "m"))]
  #d <- dat
  #dat <- dat[,.SD,.SDcols = !c("Y")]
  #X <- dat[,.SD,.SDcols = !c("Y")]
  #dat[, target := factor(res$R, labels = c("c", "m"))]
}
b <- Sys.time()
print(b-a)


# profiling
Rprof(tmp <- tempfile())
#simDat(parms)
params <- parms
replicate(100, simCorMix(params = params, Nvar = 24, prfx = "Rel"))
#fitAlgo(params, dat[,!"Y"], 1, c("a","b"))
Rprof()
summaryRprof(tmp)
