# options(scipen=999)
# rm(list=ls(all=T))
# gc()
if (!grepl("(.*Thesis-main$)|(.*Thesis$)", getwd())){
  stop("change WD to Thesis dir!")
}

set.seed(1813544)

lapply(c(list.files(file.path("Simulation 1", "getOutput"),pattern=".*R$",full.names = T),
         file.path("Simulation 1", "setup.R")), 
       source)


ncores <- detectCores()
#Sys.setenv("OMP_THREAD_LIMIT"=ncores-1)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)
#registerDoRNG(1813544)

#options(error = dump.frames)

x <- foreach(nMc=1:10, .packages=c("mvnfast",
                                   "data.table",
                                   "splines",
                                   "caret",
                                   "ranger",
                                   "Boruta",
                                   "FSelector",
                                   "FSinR",
                                   "sparseSVM",
                                   "FCBF",
                                   "rlecuyer")) %dopar% {
                                     simRep(fixedParams, variedParams[20,], nMc)
                                   }

#res1 <- readRDS("Simulation 1/Results/mech_mnar_pm_0.5_corrPred_1_pr_0.2_rep_2.RDS")


stopCluster(cl)
registerDoSEQ()
#closeAllConnections()