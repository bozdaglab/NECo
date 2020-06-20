rm(list=ls())

source("GecoFunctions.R")

t0 <<- Sys.time() 

nbhd <- generateNbhd("input_files", 40)

trimTopNbhd(nbhd, "nbhd_ic")
cat(format(Sys.time()-t0), " to finish all!\n")
