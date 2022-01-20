source("simulate_profiles.R")

generate_sim_data <- function(simulation_num, std_devs, means, lengths, aber_start, num_probes, num_patients, mix) {
  for (sim in 1:simulation_num) {
    for (sd in std_devs) {
      for (mean in means) {
        for (length in lengths) {
          simulation_data = generate_data_set(mean, sd, aber_start, length, num_probes, num_patients, mix)
          file_path = paste("/Users/jkaslam/Desktop/Simulations/", "sim", sim, "stdd", sd, "mean", mean, "len", length, ".csv", sep="")
          write.csv(simulation_data, file_path, row.names = FALSE)
        }
      }
    }
  }
}