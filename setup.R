options(scipen=999)
lapply(list("Data/simData.R", "Data/simMiss.R", "helpers.R"), source)
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
                    deg = 4,
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
cm_lp <- rethinking::rlkjcorr(1, fixedParams$Ntotal*0.05, 0.0001)
diag(cm_lp) <- 1
cm_hp <- rethinking::rlkjcorr(1, fixedParams$Ntotal*0.2, 0.0001)
diag(cm_hp) <- 1
#plot((cm_hp[lower.tri(cm_hp)][order(cm_hp[lower.tri(cm_hp)])]))

# varied factors
variedParams <- expand.grid(list(mechanism = c("mcar", "mar", "mnar"),
                                 pm = c(0.1, 0.5),
                                 corrPred = c(0, 1),
                                 pr = c(0.05, 0.2)))

variedParams <- data.table::as.data.table(variedParams)

# don't vary pr for mcar remove with datatable magic
#variedParams[(mechanism=="mcar" & pr == 0.2) := NULL,]


variedParams$trans <- ifelse(variedParams$pr==0.2,
                             lapply(1:24, function(x) 
                               c(sin, exp, cos, sin, sin, sin, sin, sin,
                                            sin, exp, exp, exp, sin, exp, exp, exp)),
                             lapply(1:24, function(x) c(sin, exp, cos, sin)))

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
nms <- names(X)[grepl("^Rel", names(X))]
Y <- simY(parms, X[, ..nms])
#Y <- simY(parms, X[,grepl("^rel", colnames(X))])
summary(lm(V1~., data=cbind(X[, ..nms],Y)))
R <- simR(parms, X, Y)
target <- R$R
mean(target) 
parms$pm
target <- factor(target, labels = c("c", "m"))
require(performance)
rels <- copy(X[, ..nms])
rels[,c(9:24) := mapply({function(f, x) f(x)}, parms$trans, .SD, SIMPLIFY = FALSE),
     .SDcols = names(X)[grepl("(^RelBin)|(^RelOrd)", names(X))]]
rels[,c(1:8) := lapply(.SD, simSpline, deg = parms$deg),
     .SDcols = names(X)[grepl("^RelCont", names(X))]]
m <- glm(target~., data = as.data.frame(cbind(rels,Y,target)), family = binomial("probit"))
r2_mckelvey(m)
parms$R2r

## check McKelvey's R^2 ##

meansDT <- data.table::data.table(rep(NA, 1000),rep(NA, 1000),rep(NA, 1000),rep(NA, 1000))
r2rsDT <- data.table::data.table(rep(NA, 1000),rep(NA, 1000),rep(NA, 1000),rep(NA, 1000))

meansDT2 <- data.table::data.table(rep(NA, 200),rep(NA, 200),rep(NA, 200),rep(NA, 200))
r2rsDT2 <- data.table::data.table(rep(NA, 200),rep(NA, 200),rep(NA, 200),rep(NA, 200))
# rows 

i1 <- c(fixedParams, variedParams[pr==0.05 & mechanism=="mar" & pm==0.1 & corrPred==1,])
i2 <-  c(fixedParams, variedParams[pr==0.05 & mechanism=="mar" & pm==0.5 & corrPred==1,])
i3 <-  c(fixedParams, variedParams[pr==0.05 & mechanism=="mnar" & pm==0.1 & corrPred==1,])
i4 <-  c(fixedParams, variedParams[pr==0.05 & mechanism=="mnar" & pm==0.5 & corrPred==1,])

countr <- 1

test.conds <- expand.grid(list(pm=c(0.1, 0.5),
                               mechanism=c("mar", "mnar")))
for (i in 1:4){
  mech <- test.conds[i,"mechanism"]
  perc <- test.conds[i,"pm"]
  parms <- c(fixedParams, variedParams[pr==0.05 & mechanism==(mech) 
                                       & pm==(perc) & corrPred==1,])
  
  means <- rep(NA, 1000)
  r2r <- rep(NA, 1000)
  for (j in 1:1000){
    X <- simX(parms)
    nms <- names(X)[grepl("^Rel", names(X))]
    Y <- simY(parms, X[, ..nms])
    R <- simR(parms, X, Y)
    target <- R$R
    means[j] <- mean(target) 
    target <- factor(target, labels = c("c", "m"))
    #require(performance)
    is.mnar <- params$mechanism == "mnar"
    Xt <- cbind(data.table::copy(X), Y[,..is.mnar]) 
    contNms <- c(names(Xt)[grep("^RelCont", names(Xt))], names(Y)[is.mnar])
    Xt[, (contNms) := lapply(.SD, simSpline, deg = params$deg, coefSpl = params$theta),
           .SDcols = contNms]
    catNms <- names(Xt)[grep("(^RelBin)|(^RelOrd)", names(Xt))]
    Xt[, (catNms) := mapply({function(f, x) f(x)}, unlist(params$trans), .SD, SIMPLIFY = FALSE), 
           .SDcols = catNms]
    m <- glm(target~., data = cbind(as.data.frame(Xt),target), 
             family = binomial("probit"))
    r2r[j] <- performance::r2_mckelvey(m)
    #parms$R2r
    
  }
  
  meansDT[, (paste0("V", countr)) := means]
  r2rsDT[, (paste0("V", countr)) := r2r]
  countr <- countr + 1
}

countr <- 1
for (i in 1:4){
  mech <- test.conds[i,"mechanism"]
  perc <- test.conds[i,"pm"]
  parms <- c(fixedParams, variedParams[pr==0.05 & mechanism==(mech) 
                                       & pm==(perc) & corrPred==1,])
  
  means <- rep(NA, 200)
  r2r <- rep(NA, 200)
  for (j in 1:200){
    X <- simX(parms)
    nms <- names(X)[grepl("^Rel", names(X))]
    Y <- simY(parms, X[, ..nms])
    R <- simR(parms, X, Y)
    target <- R$R
    means[j] <- mean(target) 
    target <- factor(target, labels = c("c", "m"))
    #require(performance)
    is.mnar <- params$mechanism == "mnar"
    Xt <- cbind(data.table::copy(X), Y[,..is.mnar]) 
    contNms <- c(names(Xt)[grep("^RelCont", names(Xt))], names(Y)[is.mnar])
    Xt[, (contNms) := lapply(.SD, simSpline, deg = params$deg, coefSpl = params$theta),
       .SDcols = contNms]
    catNms <- names(Xt)[grep("(^RelBin)|(^RelOrd)", names(Xt))]
    Xt[, (catNms) := mapply({function(f, x) f(x)}, unlist(params$trans), .SD, SIMPLIFY = FALSE), 
       .SDcols = catNms]
    m <- glm(target~., data = cbind(as.data.frame(Xt),target), 
             family = binomial("probit"))
    r2r[j] <- performance::r2_mckelvey(m)
    #parms$R2r
    
  }
  
  meansDT2[, (paste0("V", countr)) := means]
  r2rsDT2[, (paste0("V", countr)) := r2r]
  countr <- countr + 1
}


saveRDS(meansDT[, lapply(.SD, mean)], "mean pm over 1000 repetitions")
saveRDS(r2rsDT[,lapply(.SD, mean)], "mean R2r over 1000 repetitions")

saveRDS(meansDT[, lapply(.SD, mean)], "mean pm over 1000 repetitions scaled preds")
saveRDS(r2rsDT[,lapply(.SD, mean)], "mean R2r over 1000 repetitions scaled preds")

hist(means)
plot(density(means))
mean(means)
sd(means)
hist(r2rs)
plot(density(r2rs))
mean(r2rs)

params <- parms


#logisZ <- relt[, (Nrel - Nrel / 3 + 1): Nrel, drop = F]
#N <- params$N
#nvar <- (Nrel / 3)


a <- Sys.time()
# check if it works
for (i in 1:24){
  params <- c(fixedParams, variedParams[i,])
  X <- simX(params)
  nms <- names(X)[grepl("^Rel", names(X))]
  Y <- simY(params, X[,..nms])
  R <- simR(params, X, Y)
}
b <- Sys.time()
print(b-a)

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
  nms <- names(X)[grepl("^Rel", names(X))]
  Y <- simY(parms, X[, .SD, ,.SDcols = nms])
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
  

}

# also check each aspect but with replication and average to see how it behaves in the limit
# if it approach pop values

difs[1:6,1:6,which(variedParams$corrPred == 1 & variedParams$pr == 0.05)]
difs[1:24,1:24,which(variedParams$corrPred == 1 & variedParams$pr == 0.2)]



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
meansBin1
meansBin2
ord1 # 8th repetition, 2nd variable is suspicious? two times two biases opposite
# same 14th repetition, 2nd variable
ord2

#length(params$cormat) == (24*23)/2
