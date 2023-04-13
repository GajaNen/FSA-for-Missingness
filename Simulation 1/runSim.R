#options(scipen=999)
#rm(list=ls(all=T))
#gc()
if (!grepl("(.*Thesis-main$)|(.*Thesis$)", getwd())){
  stop("change WD to Thesis dir!")
}
lapply(c(list.files(file.path("Simulation 1", "getOutput"),pattern=".*R$",full.names = T),
         file.path("Simulation 1", "setup.R")), 
       source)


ncores <- detectCores()
#Sys.setenv("OMP_THREAD_LIMIT"=ncores-1)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)
#registerDoRNG(1813544)

#options(error = dump.frames)
set.seed(1813544)

x <- foreach(nmc=1:10, .packages=c("mvnfast",
                                   "data.table",
                                   "splines",
                                   "caret",
                                   "ranger",
                                   "Boruta",
                                   "FSelector",
                                   "FSinR",
                                   "sparseSVM",
                                   "Biocomb")) %dorng% {
                                     simRep(fixedParams, variedParams,nmc)
                                   }

stopCluster(cl)
registerDoSEQ()
#closeAllConnections()