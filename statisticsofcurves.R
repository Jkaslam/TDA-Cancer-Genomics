# Takes in a matrix of persistence curves and indices of rows that are the 
# case and control persistence curves. Computes average persistence curves for
# cases and controls, then computes the p-value for the L1 and L2 norm
# of the mean and control curves. 
statisticsofcurves <- function(curvematrix, curvetype, cases, controls) {
  numpatients = length(cases) + length(controls)
  
  # Separates the persistence curves into cases and controls
  casesmatrix <- curvematrix[cases,]
  controlsmatrix <- curvematrix[controls,]
  
  # Calculates the average persistence curve for cases and control
  meancases=colMeans(casesmatrix)
  meancontrols=colMeans(controlsmatrix)
  
  ymax=max(max(meancases),max(meancontrols))
  
  # Computes the original statistics (L1 and L2 norms) values for comparison 
  #to compute p-value. 
  original_statistic=sum(abs(meancases-meancontrols))
  original_statistic2=sum((meancases-meancontrols)^2)
  
  # Calculates the p-values of the the L1 and L2 norm statistics. 
  permutation_statistic=rep(0,10000)
  permutation_statistic2=rep(0,10000)
  
  for (i in 1:10000){
    new_cases=sample(numpatients, length(cases))
    new_controls=which(!(1:numpatients %in% new_cases))
    
    mean_new_cases <- colMeans(curvematrix[new_cases,])
    mean_new_controls <- colMeans(curvematrix[new_controls,])
    
    permutation_statistic[i] <- sum(abs(mean_new_cases-mean_new_controls))
    permutation_statistic2[i] <- sum((mean_new_cases-mean_new_controls)^2)
  }
  
  pvalue <- mean(permutation_statistic > original_statistic)
  pvalue2 <- mean(permutation_statistic2 > original_statistic2)

  return(list(pvalue, pvalue2, meancases, meancontrols))
}

