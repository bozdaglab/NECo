rm(list=ls())

source("GecoFunctions.R")

t0 <<- Sys.time() 
inputFileName <- "input_files_ic"
nbhd <- generateNbhd(inputFileName, 40)

trimTopNbhd(nbhd, "nbhd_ic")
cat(format(Sys.time()-t0), " to finish all!\n")
