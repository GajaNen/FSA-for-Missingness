datPop <- read.csv(file.path("Simulation 2", "heart_failure.csv"))
data.table::setDT(datPop)
datPop[, time := NULL]
data.table::setnames(datPop, c("RelCont1", "RelBin1", "RelCont2", "RelBin2",
                              "RelCont3", "RelBin3", "RelCont4", "RelCont5", 
                              "RelCont6", "RelBin4", "RelBin5", "Y"))
