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
# What percentage of patients in the test set have an aberration
mix = as.numeric(args[3])
print(mix)
# Generates the simulation data
sim_path = args[1]
#generate_sim_data(simulation_num, std_devs, means, lengths, aber_start, num_probes, num_patients, mix, sim_path)
print("Finished generating simulated data")

maxfiltrationvalue= 2
filtrationvector = seq(from=0.01,to=maxfiltrationvalue,by=0.01)

for (sd in std_devs) {
  for (mean in means) {
    for (length in lengths) {
      # Initialize the matrix containing the p-values for the simulations with current parameters
      pvaluetable <- data.frame(matrix(ncol = 8, nrow = 0))
      colnames(pvaluetable) <- c("Simulation", "Mean", "Standard Deviation", "Length", "L1 P-Val", "L2 P-Val", "FDR L1", "FDR L2")
      
      # Finds the files with the data for the current simulations
      curr_pattern = paste("sim", "\\d+", "stdd", toString(sd), "mean", toString(mean), "len", toString(length), ".csv", sep = "")
      simulations_with_curr_params = list.files(path=sim_path, pattern = curr_pattern)
      
      for (sim in simulations_with_curr_params) {
        curr_sim_num = stri_extract_first_regex(sim, "\\d+")
        simulated_data = read.csv(paste(sim_path, sim, sep=""))
        test_cases = 1:num_test_patients
        control_cases = (num_test_patients + 1):nrow(simulated_data)

        # Since rows represent patients we have to transpose since TDA signal was built
        # for columns representing patients
        TDAsignaloutput <- apply(t(simulated_data), 2, TDAsignal)
        list_profiles<-lapply(TDAsignaloutput, persistentcurvemaker, par_discretization=filtrationvector, dimstudy = dimofstudy)

        # initialize persistence curve matrix
        bettiprofilematrix <- matrix(data=0, nrow=num_patients, ncol=length(filtrationvector))
        lifespanmatrix <- bettiprofilematrix

        # fill the betti and lifespan matrices
        for (i in 1:length(list_profiles)) {
          #print(i)
          bettiprofilematrix[i,] = list_profiles[[i]]$betti1
          lifespanmatrix[i,] = list_profiles[[i]]$lifespan
        }

        #curvetype <- "betti0"
        #curvematrix<-bettiprofilematrix
        curvetype <- "lifespan"
        curvematrix<-lifespanmatrix
        statofcurvesoutput = statisticsofcurves(curvematrix, curvetype, test_cases, control_cases)
        l1andl2pvals = c(statofcurvesoutput[[1]], statofcurvesoutput[[2]])
        
        # COMMENT THE SUCCEEDING BLOCK OF CODE OUT IF YOU DON'T WANT TO PLOT BOTH THE TEST AND CONTROL CURVES
        # Plot the average test and control curves to determine whether the test curve is on top or not
        # allowing us to eliminate regions detected as significant where the control set caused the
        # detection of significance.
        #plotpath = paste("/Users/jkaslam/Desktop/Simulations/betti0-plots/", substr(sim, 1, nchar(sim) - 4), "plot", ".jpeg", sep="")
        # plotpath = paste("/Users/jkaslam/Desktop/Simulations/lifespan0-plots/", substr(sim, 1, nchar(sim) - 4), "plot", ".jpeg", sep="")
        # 
        # #print(plotpath)
        # jpeg(file = plotpath)
        # #Test in red
        # plot(filtrationvector, statofcurvesoutput[[3]], type="l", col="red")
        # # Control in blue
        # lines(filtrationvector, statofcurvesoutput[[4]], type="l",col="blue")
        # legend(1.5, .1, legend=c("Test", "Control"), col=c("red", "blue"), lty = 1:1)
        # title(paste("Lifespan Curve", "stdd", toString(sd), "mean", toString(mean), "len", toString(length), sep = " "))
        # #title(paste("Betti Curve", "stdd", toString(sd), "mean", toString(mean), "len", toString(length), sep = " "))
        # dev.off()
        
        # Save the p-values
        pvaluetable[nrow(pvaluetable) + 1, ] = list(curr_sim_num, mean, sd, length, l1andl2pvals[1], l1andl2pvals[2], 0, 0)
        
        # Apply FDR to the p-values
        pvaluetable[, 7] = p.adjust(pvaluetable[, 5], method="fdr")
        pvaluetable[, 8] = p.adjust(pvaluetable[, 6], method="fdr")
      
      }
      
      # Write p-values to a file
      pval_path = paste(args[2], "mix", mix, "mean", mean, "stdd", sd, "len", length, "pvalslifespan", ".csv", sep="")
      write.csv(pvaluetable, pval_path)
      #print("wrote a csv file")
    }
  }
}
