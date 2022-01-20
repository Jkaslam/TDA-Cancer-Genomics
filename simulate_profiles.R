# Generates a profile with an aberration with given mean and standard deviation. 
# Aberration starts at aber_start and is as long as aber_length. 
generate_test_profile <- function(mean, standard_deviation, num_probes, aber_start, aber_length) {
  normal_section_1 = rnorm(aber_start - 1, mean = 0, sd = standard_deviation)
  aber_section = rnorm(aber_length, mean = mean, sd = standard_deviation)
  normal_section_2 = rnorm(num_probes - (aber_start + aber_length - 1), mean = 0, sd = standard_deviation)
  return(c(normal_section_1, aber_section, normal_section_2))
}

# Generates a control profile with default mean of 0 and given standard deviation.
generate_control_profile <- function(mean = 0, standard_deviation, num_probes) {
  return(rnorm(num_probes, mean = mean, sd = standard_deviation))
}

# Generates a data set with given mean for aberrations. Half of dataset are from the test set
# and the other half are from the control set. Only mix % of profiles in the test set are aberrant profiles. 
generate_data_set <- function(mean, standard_deviation, aber_start, aber_length, num_probes, num_patients, mix) {
  patients <- matrix(data=0, nrow = num_patients, ncol = num_probes)
  aberrant_profs_in_test = as.integer(num_patients/2 * mix)
  
  for (i in 1:num_patients){
    if (i <= aberrant_profs_in_test) {
      patients[i,] = generate_test_profile(mean, standard_deviation, num_probes, aber_start, aber_length)
    }
    else {
      patients[i,] = generate_control_profile(standard_deviation = standard_deviation, num_probes = num_probes)
    }
  }
  
  return(patients)
}
  
