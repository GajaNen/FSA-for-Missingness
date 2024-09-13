# options(scipen=999)
rm(list=ls(all=T))
gc()
if (!grepl("(.*Thesis-main$)|(.*Thesis$)", getwd())){
  stop("change WD to Thesis dir!")
}

set.seed(1813544)

lapply(c(file.path("Simulation 2", "Preprocessing", "preprocess.R"),
         file.path("Simulation 2", "Scripts", "depend.R"),
         list.files(file.path("Simulation 2", "Scripts", "getOutput"),pattern=".*R$",full.names = T),
         file.path("Simulation 2", "Scripts", "setup.R")), 
       source)

# if the directory for results doesn't exist, create one
if (!dir.exists(fixedParams$dir)) dir.create(fixedParams$dir)

ncores <- detectCores()
Sys.setenv("OMP_THREAD_LIMIT"=1)
cl <- makeCluster(ncores-2)
clusterEvalQ(cl = cl, {
  library(data.table) # make sure we load the package on the cluster
  setDTthreads(1)
})
registerDoParallel(cl)

a <- Sys.time()

# modify nMC for a smaller number of repetitions, e.g. nMc=1:10
x <- foreach(nMc=1:500, .packages=c("data.table",
                                    "caret",
                                    "ranger",
                                    "Boruta",
                                    "FSelector",
                                    "FSinR",
                                    "sparseSVM",
                                    "FCBF",
                                    "rlecuyer")) %dopar% {
                                      simRep(dt = datPop, 
                                             fixed = fixedParams, 
                                             varied = variedParams, 
                                             rpt = nMc)
                                    }

b <- Sys.time()

print(b-a)

stopImplicitCluster()
stopCluster(cl)
registerDoSEQ()
#closeAllConnections()
