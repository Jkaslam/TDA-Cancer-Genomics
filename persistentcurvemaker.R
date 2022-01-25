# Computes the Betti and lifespan curves from the given interval matrix. 
persistentcurvemaker <- function(oneTDAsignaloutput, dimstudy, par_discretization=seq(from=0,to=2,by=0.01)) {
  IntervalMatrix <-  oneTDAsignaloutput$intervals
  bettiprofile <- rep(0, length(par_discretization))
  lifespanprofile <- rep(0,length(par_discretization))
  lifespan_vector <- IntervalMatrix[,"Death"]-IntervalMatrix[,"Birth"]
  total_lifespan <- sum(lifespan_vector)
  if (nrow(IntervalMatrix>0)){
    # We will do calculations only between the minimum value of the left hand side of the
    # intervals and the maximum value of the right hand side.
    if (length(IntervalMatrix[,1]) == 0)
    {
      minimum_value <- 1
    }
    else{
      minimum_value <- max(1, which(par_discretization < min(IntervalMatrix[,1])))
      #adding 1 in the max here is necessary for intervals starting at 0 like in betti0, like length of (par_dicretization).
    }
    maximum_value <- min(length(par_discretization), min(which(par_discretization>max(IntervalMatrix[,2]))))
    for (i in minimum_value:maximum_value){
      intervals_containing <- which((!(IntervalMatrix[,1]>par_discretization[i])) & par_discretization[i]<IntervalMatrix[,2])
      bettiprofile[i] <- length(intervals_containing)
      if (dimstudy == 0) {
        if (sum(lifespan_vector[intervals_containing]) - max(par_discretization) < 0) {
          print(i)
          lifespanprofile[i] = 0
        } else {
          lifespanprofile[i] <- sum(lifespan_vector[intervals_containing]) - max(par_discretization)
        }
      }
      else {
        lifespanprofile[i] <- sum(lifespan_vector[intervals_containing])
      }
    }
  }
  else 
  { 
    total_lifespan <- 0 
  }

  profile_list <- list(betti1=bettiprofile,lifespan=lifespanprofile,totallifespan=total_lifespan)
  return(profile_list)
}