source("GecoFunctions.R")

# "input_files.txt" and number of cores for parallel processing
# generate walk matrix and runs RWR to generate neighborhood with RWR scores
nbhd <- generateNbhd("input_files", 40)

# Generates different top N neighborhoods using nbhd data frame, 2nd parameter is the prefix for output file name
# for different top N parameters provide following parameters
# top_Ggs=c(150, 250), top_Gps=c(150, 250), top_Pgs=c(150, 250), top_Pps=c(150, 250)
trimTopNbhd(nbhd, "nbhd")
