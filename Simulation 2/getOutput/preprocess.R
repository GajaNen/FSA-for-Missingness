datPop <- data.table::fread(file=file.path("Simulation 2", "heart_failure.csv"))
datPop[, c("time", "DEATH_EVENT") := NULL]
data.table::setnames(datPop, c("IrrelCont1", "IrrelBin1", "IrrelCont2", "IrrelBin2",
                              "IrrelCont3", "IrrelBin3", "IrrelCont4", "IrrelCont5", 
                              "IrrelCont6", "RelBin1", "Y"))
length(datPop)
