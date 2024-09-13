# options(scipen=999)
rm(list=ls(all=T))
gc()
if (!grepl("(.*Thesis-main$)", getwd())){
  stop("change WD to Thesis-main dir!")
}

set.seed(1813544)

lapply(c(file.path("Simulation 1", "Scripts", "depend.R"),
         file.path("Simulation 1", "Scripts", "getOutput", "helpers.R"),
         file.path("Simulation 1", "Scripts", "getOutput", "simMiss.R"),
         file.path("Simulation 1", "Scripts", "getOutput", "simXY.R"),
         file.path("Simulation 1", "Scripts", "getOutput", "runAlgo_FCBF_rerun.R"),
         file.path("Simulation 1", "Scripts", "setup.R")), 
       source)

fixedParams$dir <- file.path("Simulation 1", "Results_rerun")
fixedParams$prevdir <- file.path("Simulation 1", "Results")

# create the directory for the results if it doesn't exist yet
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

# modify nMc to run a subset of repetitions, e.g. nMc=1:10
x <- foreach(nMc=1:500, .packages=c("mvnfast",
                                    "data.table",
                                    "splines",
                                    "FCBF",
                                    "rlecuyer")) %dopar% {
                                      simRep(fixedParams, variedParams, nMc)
                                    }

b <- Sys.time()

print(b-a)

stopImplicitCluster()
stopCluster(cl)
registerDoSEQ()
#closeAllConnections()
