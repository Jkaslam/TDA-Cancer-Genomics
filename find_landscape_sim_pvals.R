#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)
library(stringi)

source("generate_simulation_data.R")
source("TDAsignal.R")
source("persistentcurvemaker.R")
source("statisticsofcurves.R")
source("persistencelandscapefromoutputTDA.R")

# The parameters for the simulations we want to run
# How many times we want to run the simulation
simulation_num = 50
# List of standard deviations we want to run simulations for
std_devs = seq(from = 0.2, to = 0.5, by = .02)
# List of means for aberrations we want to run simulations for
means = c(-1, .6, 1)
# List of lengths of aberrations we want to run simulations for
lengths = c(1, 2, 3, 5, 10, 15)
# The location we want the aberrations in the test set to start at
# assumes aber_start + any length in lengths is not > num_probes
aber_start = 1
# Num of probes in a profile
num_probes = 20
# Num of patients in a simulation (assumes number is even)
num_patients = 120
#Always have half test and half control
num_test_patients = num_patients/2
#Dimension of study
dimofstudy = 0
#Number of persistence landscape functions to test
max_k = 4
#
mix = as.numeric(args[3])
print(mix)
# Generates the simulation data, comment out if you already have it
#generate_sim_data(simulation_num, std_devs, means, lengths, aber_start, num_probes, num_patients, mix)

maxfiltrationvalue=2
filtrationvector = seq(from=0.01,to=maxfiltrationvalue,by=0.01)

sim_path = args[1]
# Loop through all combinations of standard deviations, means and lengths of
# the aberration 
for (sd in std_devs) {
  for (mean in means) {
    for (length in lengths) {
      # Make a p-value table for each persistent landscape function
      pvaluetables = vector(mode = "list", max_k)
      for (i in 1:max_k) {
        pvaluetable <- data.frame(matrix(ncol = 11, nrow = 0))
        colnames(pvaluetable) <- c("Simulation", "Mean", "Standard Deviation", "Length", "L1 P-Val", "L2 P-Val", "FDR L1", "FDR L2")        
        pvaluetables[[i]] = pvaluetable
      }
      
      # Finds the files with the data for the current simulations
      curr_pattern = paste("sim", "\\d+", "stdd", toString(sd), "mean", toString(mean), "len", toString(length), ".csv", sep = "")
      simulations_with_curr_params = list.files(path= sim_path, pattern = curr_pattern)
      
      for (sim in simulations_with_curr_params) {
	#print(sim)
        curr_sim_num = stri_extract_first_regex(sim, "\\d+")
        simulated_data = read.csv(paste(sim_path, sim, sep=""))
        test_cases = 1:num_test_patients
        control_cases = (num_test_patients + 1):nrow(simulated_data)
        
        # Since rows represent patients we have to transpose since TDA signal was built
        # for columns representing patients
        TDAsignaloutput <- apply(t(simulated_data), 2, TDAsignal)
        diagrams = vector(mode = "list", length(TDAsignaloutput))
        landscape_matrices = vector(mode = "list", length(TDAsignaloutput))
        
        # For each persistent diagram, find up to max_k persistent landscape functions
        # from that persistent diagram. 
        for (i in 1:length(diagrams)) {
          diagrams[[i]] = TDAsignaloutput[[i]][[3]]
          landscape_matrices[[i]] = persistencelandscapefromoutputTDA(diagrams[[i]], curve_dimension = dimofstudy, maxscale = 2, par_discretization = filtrationvector, Kmax = max_k)
          #print(landscape_matrices[[i]])
        }
        
        for (i in 1:max_k) {
          landscapes = matrix(0,nrow=length(diagrams),ncol=length(filtrationvector))
          for (k in 1:length(landscape_matrices)) {
            landscapes[k, ] = landscape_matrices[[k]][i,]
          }
          #print(landscapes)
          statofcurvesoutput = statisticsofcurves(landscapes, "betti", test_cases, control_cases)
          l1andl2pvals = c(statofcurvesoutput[[1]], statofcurvesoutput[[2]])
          num_rows = nrow(pvaluetables[[i]])
          pvaluetables[[i]][num_rows + 1, ] = list(curr_sim_num, mean, sd, length, l1andl2pvals[1], l1andl2pvals[2], 0, 0)
        }
        
      }
      
      for (curr_landscape in 1:max_k)
      {
        #Apply FDR to the p-values
        pvaluetables[[curr_landscape]][, 7] = p.adjust(pvaluetables[[curr_landscape]][, 5], method="fdr")
        pvaluetables[[curr_landscape]][, 8] = p.adjust(pvaluetables[[curr_landscape]][, 6], method="fdr")
        
        path = paste(args[2], "mix", mix, "landscapek", curr_landscape, "mean", mean, "stdd", sd, "len", length, "pvals", ".csv", sep="")
        write.csv(pvaluetables[[curr_landscape]], path)
      }
      
    }
  }
}
