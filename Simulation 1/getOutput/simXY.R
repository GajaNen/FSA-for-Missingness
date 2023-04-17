###--------------------------------------------------------------------------###

# generate ordinal variable from quantiles of logistic distribution given
# baseline (population) probabilities of each category

#@logisZ a vector of log. quantiles
#@lP: population probabiltiies
#@N: number of samples
#returns a numeric vector of categories

simOrd <- function(logisZ, lP, N){
  
  cprop <- cumsum(lP)
  if (max(cprop) != 1){
    stop("Population probabilities of each ordinal variable must sum to 1.")
  }
  quant <- stats::qlogis(cprop)
  matlp <- matrix(rep(quant, N),ncol = length(quant),byrow = TRUE)
  grp <- (logisZ > cbind(-Inf, matlp))
  cats <- apply(grp, 1, sum)
  return(cats)
}

###--------------------------------------------------------------------------###

# generate (ir)relevant variables (and Y if addY)
# of mixed types with a given correlation structure
# and marginal distributions with a Gaussian copula 
# (only non-parametric correlation structure retained)

# this function can have a side effect of modifying the original DT if this
# DT is an input together with the names of relevant variables
# otherwise a new DT is constructed within function

#@prfx: used to access parameters (Rel or Irrel)
#@Nsim: number of vars to be generated (at least one of Nsim or nms must be specified)
#@dts: DT to be modified by reference, if NULL a new DT of params$N * Nsim created
#@nms: colnames of variables from DT to be considered, determines Nsim if Nsim missing
# if existing DT is passed, they either refer to existing columns or new cols are 
# created, if not they are the automatically created names (V,1:Nsim)

#return DT

simCorMix <- function(params, prfx, Nsim=NULL, dts=NULL, nms=NULL){
  
  #useful checks when specifying population parameters (when OK no need to check it in each rep)
  # if (!identical(dim(params[[paste0("corMat",prfx)]][[1]])[1],length(params[[paste0("params",prfx)]][[1]]),
  #     length(params[[paste0("dists",prfx)]][[1]]))) {
  #   stop(pasteo("Cormat, MargParam and DistNames don't match for ", prfx, " vars"))
  # }
  if (is.null(Nsim)&&is.null(nms)) {
    stop("Provide number or names of vars to be generated!")
  } else if (is.null(Nsim)) Nsim <- length(nms)
  # if (length(params[[paste0("popProbs", prfx)]][[1]]) != (Nsim %/% 3)){
  #   stop(paste0("number of probability vectors for odinal must be a third of ",prfx,"evant vars."))
  # }
  if (is.null(dts)) { # if no DT provided, create a new one
    dts <- data.table::setDT(lapply(1:Nsim, rep, NA, params$N))
    nms <- names(dts)
  } #else if (is.null(nms)) stop("Provide names for variables to be modified in dts.")
  
  dts[, (nms) := # get uniforms from MV standard normals
        data.table::as.data.table(stats::pnorm(
          mvnfast::rmvn(n=1000,mu=rep(0,Nsim),
                        sigma=params[[paste0("corMat",prfx)]][[1]])))
  ]
  dts[, (nms) := # apply appropriate quantile funcs to the uniforms to get desired distr
        mapply(function(x,y,z) params$map.funcs[[y]](x, z[1],z[2]), 
               .SD,  params[[paste0("dists",prfx)]][[1]], 
               params[[paste0("params",prfx)]][[1]],
               SIMPLIFY = FALSE),
      .SDcols = nms]
  # locate logistic variables (may not be named with Ord, but must be last 1/3 of vars)
  n.pt <- Nsim %/% 3 # number of X vars per type 
  ords <- nms[(n.pt*2+1):(n.pt*3)]
  dts[, (ords) := mapply(simOrd, .SD, params[[paste0("popProbs", prfx)]][[1]], 
                         MoreArgs = list(N=params$N), SIMPLIFY = F), 
      .SDcols = ords] # generate ordinal variables
  return(dts)
}

###--------------------------------------------------------------------------###

# simulate correlated or independent relevant and irrelevant predictors (and Y if addY)
# return DT of all variables (Xrel, Xirrel & Y if addY)

simDat <- function(params){
  
  Nrel <- params$pr * params$Ntotal
  Nirrl <- params$Ntotal - Nrel 
  if ((Nrel %% 3) || (Nirrl %% 3)) {
    stop("Adjust Ntotal or pr so that Nrel and Nirrl are a multiple of 3.")
  } # always 1/3 Cont, Bin, Ord each!
  out <- data.table::setDT(lapply(1:(params$Ntotal+params$addY),function(x) numeric(params$N)))
  nms <-  mapply(function(x,y) paste0(x, rep(c("Cont", "Bin", "Ord"),each=y/3),1:(y/3)),
                 c("Rel", "Irrel"), c(Nrel, Nirrl)) # name rel & irrel vars
  data.table::setnames(out, c(nms$Rel, "Y"[params$addY], nms$Irrel))
  simCorMix(params = params, prfx = "Rel", dts=out, nms=c(nms$Rel,"Y"[params$addY])) # sim rel vars + y if addY
  simCorMix(params = params, prfx = "Irrel", dts = out, nms=nms$Irrel) # sim irrel vars
  return(out)
}

###--------------------------------------------------------------------------###