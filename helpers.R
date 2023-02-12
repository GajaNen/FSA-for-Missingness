creatDumm <- function(dat) data.frame(apply(dat, 2,function(x) {model.matrix(~factor(x))[,-1]}))


extrReg <- function(rstr, symb="\\d+") as.numeric(unique(unlist(regmatches(rstr, regexpr(symb, rstr))))) 
