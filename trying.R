meansDT <- data.table::data.table(rep(NA, 1000),rep(NA, 1000),rep(NA, 1000),rep(NA, 1000))
r2rsDT <- data.table::data.table(rep(NA, 1000),rep(NA, 1000),rep(NA, 1000),rep(NA, 1000))
r2ymnar <- data.table::data.table(rep(NA, 1000),rep(NA, 1000))

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
  r2mnar <- rep(NA, 1000)
  for (j in 1:1000){
    X <- simX(parms)
    nms <- names(X)[grepl("^Rel", names(X))]
    Y <- simY(parms, X[,..nms])
    R <- simR(parms, X, Y)
    target <- R$R
    means[j] <- mean(target) 
    target <- factor(target, labels = c("c", "m"))
    #require(performance)
    is.mnar <- parms$mechanism == "mnar"
    Xt <- cbind(data.table::copy(X), data.table::copy(Y[,..is.mnar]))
    contNms <- c(names(Xt)[grep("^RelCont", names(Xt))], names(Y)[is.mnar])
    Xt[, (contNms) := lapply(.SD, simSpline, deg = parms$deg, coefSpl = unlist(parms$theta)),
       .SDcols = contNms]
    catNms <- names(Xt)[grep("(^RelBin)|(^RelOrd)", names(Xt))]
    Xt[, (catNms) := mapply({function(f, x) f(x)}, unlist(parms$trans), .SD, SIMPLIFY = FALSE), 
       .SDcols = catNms]
    allnms <- c(contNms, catNms)
    m <- glm(target~., data = cbind(as.data.frame(Xt[,..allnms]),target), 
             family = binomial("probit"))
    r2r[j] <- performance::r2_mckelvey(m)
    if (mech=="mnar"){
      m2 <- glm(target~Y, data = cbind(as.data.frame(Xt[,Y]),target), 
               family = binomial("probit"))
      r2mnar[j] <- performance::r2_mckelvey(m2)
    }
    #parms$R2r
    
  }
  
  meansDT[, (paste0("V", countr)) := means]
  r2rsDT[, (paste0("V", countr)) := r2r]
  if (mech == "mnar") r2rsDT[, (paste0("V", countr)) := r2mnar]
  countr <- countr + 1
}


meansDT[, lapply(.SD, mean)]
r2rsDT[,lapply(.SD, mean)]



means <- rep(NA, 1000)
r2rs <- rep(NA, 1000)

for (i in 1:1000){
  X <- simX(parms)
  nms <- names(X)[grepl("^Rel", names(X))]
  Y <- simY(parms, X[, ..nms])
  R <- simR(parms, X, Y)
  target <- R$R
  means[i] <- mean(target) 
  target <- factor(target, labels = c("c", "m"))
  #require(performance)
  rels <- copy(X[, ..nms])
  catnms <- names(X)[grepl("(^RelBin)|(^RelOrd)", names(X))]
  rels[, catnms := mapply({function(f, x) f(x)}, parms$trans, .SD, SIMPLIFY = FALSE),
       .SDcols = catnms]
  rels[,c(1:8) := lapply(.SD, simSpline, deg = parms$deg),
       .SDcols = names(X)[grepl("^RelCont", names(X))]]
  m <- glm(target~., data = as.data.frame(cbind(rels,Y,target)), family = binomial("probit"))
  r2rs[i] <- r2_mckelvey(m)
  #parms$R2r
  
}
