d <- simDat(c(fixedParams, variedParams[1,]))
nam <- copy(names(d))

x <- getPF(c(fixedParams, variedParams[1,]), 500,nam)
x2 <- getPF(c(fixedParams, variedParams[20,]), 500,nam)
x3 <- getPF(c(fixedParams, variedParams[17,]), 500,nam)

