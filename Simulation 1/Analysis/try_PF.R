
conds <- c(fixedParams, variedParams[1,])
d <- simDat(conds)
nam <- copy(names(d[,!"Y"]))
x <- getPF(c(fixedParams, variedParams[1,]), 500,nam)
conds <- c(fixedParams, variedParams[20,])
d <- simDat(conds)
nam <- copy(names(d[,!"Y"]))
x2 <- getPF(c(fixedParams, variedParams[20,]), 500,nam)
conds <- c(fixedParams, variedParams[17,])
d <- simDat(conds)
nam <- copy(names(d[,!"Y"]))
x3 <- getPF(c(fixedParams, variedParams[17,]), 500,nam)

acc_obj <- lapply(list(x, x2, x3), function(z) list(z[["accMean"]], z[["accSD"]]))


# do this first for accuracies
accsM <- setDT(lapply(1:6, function(x) numeric(20)))
accsSD <- setDT(lapply(1:6, function(x) numeric(20)))
pfM <- setDT(lapply(1:8, function(x) numeric(20)))
pfSD <- setDT(lapply(1:8, function(x) numeric(20)))
setnames(accsM, c("l1SVM ", "rfRFE", "LR", "lasso" , "EN", "RF"))
setnames(accsSD, c("l1SVM ", "rfRFE", "LR", "lasso" , "EN", "RF"))
setnames(pfM,  c("lasso", "EN" , "rfRFE", "l1SVM", "FCBC_0",
                 "FCBC_0.1", "FCBC_0.3", "FCBC_0.5"))
setnames(pfSD, c("lasso", "EN" , "rfRFE", "l1SVM", "FCBC_0",
                 "FCBC_0.1", "FCBC_0.3", "FCBC_0.5"))
for (i in 1:nrow(variedParams)){
  cds <- c(fixedParams, variedParams[i,])
  d <- simDat(cds)
  nam <- names(d[,!"Y"])
  rs <- getSubs(cds,500,nam)
  accsM[i, names(accsM) := rs[["accMean"]]]
  accsSD[i, names(accsSD) := rs[["accSD"]]]
  pfM[i, names(pfM) := rs[["PFM"]]]
}

x <- cbind(accsM, conds)

objAcc <- melt(cbind(accsM, variedParams[,1:4]), measure.vars = c("l1SVM ", "rfRFE", "LR", "lasso" , "EN", "RF"),
     variable.name = "algorithm", value.name = "performance")

objSD <- melt(accsSD, measure.vars = c("l1SVM ", "rfRFE", "LR", "lasso" , "EN", "RF"),
              variable.name = "algorithm", value.name = "performance")

obj <- melt(cbind(pfM[,c("lasso", "EN" , "rfRFE", "l1SVM")], variedParams[,1:4]), 
            measure.vars = c("lasso", "EN" , "rfRFE", "l1SVM"),
            variable.name = "algorithm", value.name = "performance")

obj[, cond := ifelse(pm==0.1&corrPred==0,
                     1, ifelse(pm==0.5&corrPred==1,
                               2, ifelse(pm==0.1&corrPred==0, 3,
                                         ifelse(pm==0.05&corrPred==1, 4, NA))))]

#obj[, c("pm","pr") :=]
# accsM[, algo := list(names(accsM))]
# accsSD[, algo := list(names(accsSD))]
obj <-  accsM
#obj$algo <- rownames(obj)

conds <- variedParams[,1:4]
d <- simDat(conds)
nam <- copy(names(d[,!"Y"]))
nrpt <- 500
nms <- nam
