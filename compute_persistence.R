#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)
library(stringi)

source("TDAsignal.R")

home_path = args[1]
dimofstudy = 1

simulations_with_curr_params = list.files(path=home_path, pattern = "*.csv$")
#print(length(simulations_with_curr_params))

for (sim in simulations_with_curr_params) {
  #print(paste(home_path, sim, sep=""))
  simulated_data = read.csv(paste(home_path, sim, sep=""))
  
  # Since rows represent patients we have to transpose since TDA signal was built
  # for columns representing patients
  TDAsignaloutput <- apply(t(simulated_data), 2, TDAsignal)

   diagrams = vector(mode = "list", 2 * length(TDAsignaloutput))
        
        for (i in 1:(length(diagrams)/2)) {
          currTDASignalOutput = TDAsignaloutput[[i]][[3]]
          currTDASignalOutput = currTDASignalOutput[currTDASignalOutput[, 1] == dimofstudy,]
          diagrams[[2*i]] = currTDASignalOutput
          diagrams[[2*i - 1]] = c(0, 0, 0)
        }
  
  all_diags = do.call(rbind, diagrams)
  persistence_path = paste(home_path, "Persistence/Dim1/", "persistence", sim, sep = "")
  write.csv(all_diags, persistence_path)
  
}
