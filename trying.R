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
  r2mnar <- matrix(NA, nrow=1000,ncol=2)
  for (j in 1:1000){
    d <- simDat(parms)
    R <- simR(parms, d)
    target <- R$R
    prds <- R$preds
    means[j] <- mean(target) 
    target <- factor(target, labels = c("c", "m"))
    #require(performance)
    Xt <- data.table::copy(d[,..prds])
    contNms <- prds[grep("(^RelCont)|(Y)", prds)]
    Xt[, (contNms) := lapply(.SD, simSpline, deg = parms$deg, coefSpl = unlist(parms$theta)),
       .SDcols = contNms]
    catNms <- prds[grep("(^RelBin)|(^RelOrd)", prds)]
    Xt[, (catNms) := mapply({function(f, x) f(x)}, unlist(parms$trans), .SD, SIMPLIFY = FALSE), 
       .SDcols = catNms]
    m <- glm(target~., data = cbind(Xt,target), 
             family = binomial("probit"))
    r2r[j] <- performance::r2_mckelvey(m)
    if (mech=="mnar"){
      m2 <- glm(target~Y, data = cbind(Xt[,"Y"],target), 
               family = binomial("probit"))
      r2mnar[j,1] <- performance::r2_mckelvey(m2)
      m3 <- glm(target~., data = cbind(Xt[,.SD,.SDcols=!c("Y")],target), 
                family = binomial("probit"))
      r2mnar[j,2] <- performance::r2_mckelvey(m3)
    }
    #parms$R2r
    
  }
  
  meansDT[, (paste0("V", countr)) := means]
  r2rsDT[, (paste0("V", countr)) := r2r]
  if (mech == "mnar") r2rsDT[, (paste0("V", 1:2)) := r2mnar]
  countr <- countr + 1
}


meansDT[, lapply(.SD, mean)]
r2rsDT[,lapply(.SD, mean)]
r2rsDT[,lapply(.SD, mean)]

means <- rep(NA, 1000)
r2rs <- rep(NA, 1000)

for (i in 1:1000){
  d <- simDat(parms)
  R <- simR(parms, d)
  target <- R$R
  means[i] <- mean(target) 
  prds <- R$preds
  target <- factor(target, labels = c("c", "m"))
  #require(performance)
  rels <- data.table::copy(d[,..prds])
  catNms <- prds[grep("(^RelBin)|(^RelOrd)", prds)]
  contNms <- prds[grep("(^RelCont)|(Y)", prds)]
  rels[, (catNms) := mapply({function(f, x) f(x)}, unlist(parms$trans), .SD, SIMPLIFY = FALSE),
       .SDcols = catNms]
  rels[, (contNms) := lapply(.SD, simSpline, deg = parms$deg, coefSpl = unlist(parms$theta)),
       .SDcols = contNms]
  m <- glm(target~., data = cbind(rels,target), family = binomial("probit"))
  r2rs[i] <- r2_mckelvey(m)
  #parms$R2r
  
}
