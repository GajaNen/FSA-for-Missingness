###--------------------------------------------------------------------------###

# generate ordinal variable from quantiles of logistic distribution given
# baseline (population) probabilities of each category

genOrd <- function(logisZ, lP, N){
  
  cprop <- cumsum(lP)
  if (max(cprop) != 1){
    stop("Population probabilities of each ordinal variable must sum to 1.")
  }
  quant <- stats::qlogis(cprop)
  matlp <- matrix(rep(quant, N),
                  ncol = length(quant),
                  byrow = TRUE)
  grp <- (logisZ > cbind(-Inf, matlp))
  cats <- apply(grp, 1, sum)
  return(cats)
}

###--------------------------------------------------------------------------###

# generate (ir)relevant variables of mixed types and a given correlation structure
# with a Gaussian copula (only non-parametric correlation structure retained)

# this function can have a side effect of modifying the original DT if this
# DT is an input together with the names of relevant variables
# otherwise a new DT is constructed within function

simCorMix <- function(params, Nvar, prfx, addY=NULL, dts=NULL, nms=NULL){
  
  namY <- c("Y")[addY]
  S <- params[[paste0("corMat",prfx)]][[1]]
  if (!is.null(addY)) Nsim <- Nvar + 1 else Nsim <- Nvar
  if (is.null(dts))  {
    dts <- data.table::setDT(lapply(1:Nsim, rep, NA, 1000))
    nms <- names(dts)
  }
  dts[, (nms) := 
        data.table::as.data.table(pnorm(mvnfast::rmvn(1000, rep(0,Nsim), S)))]
  reps <- c(Nvar/6, Nvar/6, Nvar/3, Nvar/3, 1*addY)
  dists <- unlist(mapply(rep, names(params$map.funcs)[seq_along(reps)], reps),
                  use.names = FALSE)
  dts[, (nms) := mapply(function(x,y,z) params$map.funcs[[y]](x, z[1],z[2]), 
                            .SD, dists, c(params[[paste0("params",prfx)]][[1]],
                                          list(params$paramY)[addY]), 
                            SIMPLIFY = FALSE),
      .SDcols = nms]
  cm.reps <- cumsum(reps)
  ords <- (cm.reps[3]+1):cm.reps[4]
  probs <- lapply(1:(Nvar/3), function(x) params[[paste0("popProbs", prfx)]][[1]][x,])
  dts[, (ords) := mapply(genOrd, .SD, probs, MoreArgs = list(N=params$N), SIMPLIFY = F), 
       .SDcols = ords]
  return(dts)
}

###--------------------------------------------------------------------------###

# simulate correlated or independent relevant and irrelevant predictors

simDat <- function(params){
  
  Nrel <- params$pr * params$Ntotal
  Nirrl <- params$Ntotal - Nrel
  if ((Nrel %% 6) || (Nirrl %% 6)) stop("Number of the relevant and
                                        irrelevant features each must
                                        be a multiple of 6. Adjust Ntotal
                                        or pr.")
  out <- data.table::setDT(lapply(1:(params$Ntotal+1),
                                  function (x) rep(NA, 1000)))
  data.table::setnames(out, c(paste0("Rel", rep(c("Cont", "Bin", "Ord"),each=Nrel/3),1:(Nrel/3)),
                              "Y", # this is ugly but works
                              paste0("Irrel", rep(c("Cont", "Bin", "Ord"),each=Nirrl/3),1:(Nirrl/3))))
  nmsRel <- names(out)[grep("(^Rel)|Y", names(out))]
  simCorMix(params = params, Nvar = Nrel, prfx = "Rel", addY=TRUE, dts=out, nms=nmsRel)
  nmsIrrl <- names(out)[grep("^Irrel", names(out))]
  simCorMix(params = params, Nvar = Nirrl, prfx = "Irrel", dts = out, nms=nmsIrrl)
  return(out)

  
  
}

###--------------------------------------------------------------------------###

