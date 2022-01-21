#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Script to read in the p-values from TAaCGH simulations with different curve types
# and find the sensitivity for each of the methods, outputting this in a summary csv file.

# Where to read the simulation p-value files from
input_path = args[1]
print(input_path)
#input_path = "/share/tdacancer/jkaslam/Output/Simulation_Results/Lifespan_Curves/mix0.8/"

# When to determine if a p-value is significant (usually .05)
significance_cutoff = .05

output = data.frame(matrix(ncol = 4, nrow = 0))
colnames(output) = c("Mean", "Stdd", "Length", "Sensitivity")

simulation_file_names = list.files(path=input_path)
print(length(simulation_file_names))

for (sim_file_name in simulation_file_names) {
  sim_result = read.csv(paste(input_path, sim_file_name, sep=""))
  # Finds all simulation values that are significant and divides by the total number of simulations
  sensitivity = length(which(sim_result[,9] < significance_cutoff))/nrow(sim_result)
  output[nrow(output)+1, ] = append(sim_result[1, 3:5], sensitivity)
  output = output[order(output[, 1], output[, 2], output[, 3], decreasing = FALSE),]
}

write.csv(output, paste(input_path, "summary.csv", sep=""))
