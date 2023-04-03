options(scipen=999)
lapply(list("Data/simData.R", "Data/simMiss.R", "depend.R", "helper.R"), source)
source("Data/modRfFuncs.R")

# fixed parameters
fixedParams <- list(N=1000, 
                    trans=c(sin, exp, cos, identity),
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
                    deg = 4,
                    fcbcThres = list(0.2, 0.3, 0.4),
                    rankers = list("LR"="glm", 
                                   "linSvmRFE"="svmLinear", 
                                   "rbfSvmRFE"="svmRadial", 
                                   "rfRFE"="ranger", 
                                   "ReliefF"=NULL, 
                                   "FCBF"=NULL, 
                                   "Boruta"=NULL),
                    subsets = list("lasso"="glmnet",
                                   "EN"="glmnet", 
                                   "linSvmRFE"="svmLinear", 
                                   "linSvmSA"="svmLinear",
                                   "rbfSvmRFE"="svmRadial", 
                                   "rbfSvmSA"="svmRadial", 
                                   "rfSA"="ranger", 
                                   "rfRFE"="ranger", 
                                   "l1SVM"=NULL)
                    )


# population correlation matrices for Gaussian copula
cm_lp <- rlkjcorr(1, fixedParams$Ntotal*0.05, 0.0001)
diag(cm_lp) <- 1
cm_hp <- rlkjcorr(1, fixedParams$Ntotal*0.2, 0.0001)
diag(cm_hp) <- 1
#plot((cm_hp[lower.tri(cm_hp)][order(cm_hp[lower.tri(cm_hp)])]))

# varied factors
variedParams <- expand.grid(list(mechanism = c("mcar", "mar", "mnar"),
                                 pm = c(0.1, 0.5),
                                 corrPred = c(0, 1),
                                 pr = c(0.05, 0.2)))

variedParams <- as.data.table(variedParams)

# don't vary pr for mcar remove with datatable magic
variedParams

# fix pseudo R^2 in missingness indicators per mechanism
variedParams$R2r <-  ifelse(variedParams$mechanism == "mcar", 0.1,
                            ifelse(variedParams$mechanism == "mar", 
                                   0.6, 0.6))
# fix R^2 in incomplete variable per mechanism
variedParams$R2y <- ifelse(variedParams$mechanism == "mcar", 0.1,
                           ifelse(variedParams$mechanism == "mar", 
                                  0.3, 0.3))

# set some toy parameter for normal, gamma, binomial distributions (relevant)

variedParams$normParamRel <- lapply((fixedParams$Ntotal * variedParams$pr) / 6,
                                    function(x) lapply(1:x, 
                                                       function(y) c(mean = 0, 
                                                                     sd = 1)))

variedParams$gamParamRel <- lapply((fixedParams$Ntotal * variedParams$pr) / 6,
                                   function(x) lapply(1:x, 
                                                      function(y) c(shape = 0.5, 
                                                                    rate = 1)))

variedParams$binParamRel <- lapply((fixedParams$Ntotal * variedParams$pr) / 3,
                                    function(x) lapply(1:x, 
                                                       function(y) c(size = 1, 
                                                                     prob = 0.3)))

variedParams$popProbsRel <- lapply((fixedParams$Ntotal * variedParams$pr) / 3,
                                    function(x) matrix(rep(c(0.2, 0.1, 0.3, 0.2, 0.2), x),
                                                       nrow = x,
                                                       byrow = T)
                                   )

# set some toy parameter for normal, gamma, binomial distributions (irrelevant)

variedParams$normParamIrrel <- lapply(fixedParams$Ntotal * (1 - variedParams$pr) / 6,
                                      function(x) lapply(1:x, 
                                                         function(y) c(mean = 0, 
                                                                       sd = 1)))

variedParams$gamParamIrrel <- lapply(fixedParams$Ntotal * (1 - variedParams$pr) / 6,
                                     function(x) lapply(1:x, 
                                                        function(y) c(shape = 0.5, 
                                                                      rate = 1)))

variedParams$binParamIrrel <- lapply(fixedParams$Ntotal * (1 - variedParams$pr) / 3,
                                     function(x) lapply(1:x, 
                                                        function(y) c(size = 1, 
                                                                      prob = 0.3)))

variedParams$popProbsIrrel <- lapply(fixedParams$Ntotal * (1 - variedParams$pr) / 3,
                                     function(x) matrix(rep(c(0.2, 0.1, 0.3, 0.2, 0.2), x),
                                                        nrow = x,
                                                        byrow = T)
                                     )
# define correlation matrix which differs between percentage of relevant variables
# for the correlated predictors case

variedParams$cormat <-  ifelse(variedParams$pr == 0.05 & variedParams$corrPred, 
                               lapply(1:24, 
                                      function(x) cm_lp[lower.tri(cm_lp)]), 
                               ifelse(variedParams$corrPred, 
                                      lapply(1:24, 
                                             function(z) cm_hp[lower.tri(cm_hp)]), 
                                      rep(NA, 24)))

variedParams$marginals <- lapply(fixedParams$Ntotal * variedParams$pr,
                                 function(x) c(rep("norm", (x / 6)), 
                                               rep("gamma", (x / 6)),
                                               rep("binom", (x / 3)), 
                                               rep("logis", (x / 3)))
                                 )

parms <- c(fixedParams, variedParams[24,])
X <- simX(parms)
Y <- simY(parms, X[,grepl("^rel", colnames(X))])
summary(lm(Y~., data=X[,grepl("^rel", colnames(X))]))
R <- as.factor(rbinom(1000, 1, 0.1))
target <- R
target <- factor(target, labels = c("c", "m"))

params <- parms

# smth wrong with tuning parameters of linear svm - works with no tuning grid
# check for rbf
# the same was last time, the problem could be C==0 -- no support vectors found
# especially for mcar 

# SA is super slow --- a way to parallelize it
# or remove the ones which are less important (e.g. en, statistical relief)

# bayesian methods and other suggestion
# rfe-sve and rfe-rf are very popular and so nice to keep it
# least angle regression and others not really

#for SA got this and there's no results for sa rbf svm
# Variable differences could not be computed: 
# Not enough results to compute differences

###############################################
#################### check ####################
###############################################

#all <- c(fixedParams, variedParams)

difs <- replicate(24, matrix(NA, nrow=24, ncol=24))
meansBin1 <- matrix(NA, nrow=24, ncol=8)
meansBin2 <- matrix(NA, nrow=24, ncol=38)
ord1 <- replicate(24, matrix(NA, nrow=8, ncol=5))
ord2 <- replicate(24, matrix(NA, nrow=38, ncol=5))

for (i in 1:24){
  
  parms <- c(fixedParams, variedParams[i,])
  X <- simX(parms)
  Y <- simY(parms, X[,grepl("^rel", colnames(X))])
  if (parms$corrPred){
    m <- matrix(NA,nrow=parms$Ntotal*parms$pr,ncol=parms$Ntotal*parms$pr)
    m[lower.tri(m)] <- unlist(parms$cormat)
    diag(m) <- 1
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    if (parms$pr == 0.05){
      difs[1:(parms$Ntotal*parms$pr),1:(parms$Ntotal*parms$pr),i] <- 
        2/pi*asin(cm_lp) - cor(X[,1:(parms$Ntotal*parms$pr)], method = "kendall")
    }
    else {
      difs[1:(parms$Ntotal*parms$pr),1:(parms$Ntotal*parms$pr),i] <- 
        2/pi*asin(cm_hp) - cor(X[,1:(parms$Ntotal*parms$pr)], method = "kendall")
      }
  }
  
  Nrel <- parms$Ntotal * parms$pr
  Nrest <- parms$Ntotal - Nrel
  # check percentages of binary variables
  relBin <- X[,(Nrel / 3 + 1):(Nrel / 3 *2)]
  meansBin1[i,1:ncol(relBin)] <- (0.3 - apply(relBin, 2, mean))
  irrelBin <- X[,(Nrel+(Nrest / 3 + 1)):(Nrel+(Nrest / 3 * 2))]
  meansBin2[i,1:ncol(irrelBin)] <- (0.3 - apply(irrelBin, 2, mean))
  # check percentages of categorical variables 
  relOrd <- X[,(2*Nrel / 3 + 1):Nrel]
  ord1[1:ncol(relOrd),,i] <- 
    parms$popProbsRel[[1]] - t(apply(relOrd, 2, function(x) table(x)/length(x)))
  irrelOrd <- X[,(Nrel+(2*Nrest / 3 + 1)):(Nrel+Nrest)]
  ord2[1:ncol(irrelOrd),,i] <- 
    parms$popProbsIrrel[[1]] - t(apply(irrelOrd, 2, function(x) table(x)/length(x)))
}

# also check each aspect but with replication and average to see how it behaves in the limit
# if it approach pop values

difs[1:6,1:6,which(variedParams$corrPred == 1 & variedParams$pr == 0.05)]
difs[1:24,1:24,which(variedParams$corrPred == 1 & variedParams$pr == 0.2)]
meansBin1
meansBin2
ord1 # 8th repetition, 2nd variable is suspicious? two times two biases opposite
# same 14th repetition, 2nd variable
ord2

#length(params$cormat) == (24*23)/2
