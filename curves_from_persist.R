betti_from_persist <- function(interval, par_discretization=seq(from=0,to=2,by=0.01), dimstudy) {
  interval = interval[, 2:3, drop = FALSE]
  
  bettiprofile <- rep(0, length(par_discretization))
  if (length(interval) == 0)
  {
    minimum_value <- 1
  }
  else{
    minimum_value <- max(1, which(par_discretization < min(interval[,1])))
    #adding 1 in the max here is necessary for intervals starting at 0 like in betti0, like length of (par_dicretization).
  }
  maximum_value <- min(length(par_discretization), min(which(par_discretization>max(interval[,2]))))
  for (i in minimum_value:maximum_value){
    intervals_containing <- which((!(interval[,1]>par_discretization[i])) & par_discretization[i]<interval[,2])
    bettiprofile[i] <- length(intervals_containing)
    # if (dimstudy == 0 & bettiprofile[i] == 1) {
    #   break
    # }
    
  }
  
  return(bettiprofile)
}

lifespan_from_persist <- function(interval, par_discretization=seq(from=0,to=2,by=0.01), dimstudy) {
  interval = interval[, 2:3, drop = FALSE]
  #print(interval)
  lifespan_vector <- interval[,2]-interval[,1]
  
  lifespanprofile <- rep(0, length(par_discretization))
  if (length(interval) == 0)
  {
    minimum_value <- 1
  }
  else{
    minimum_value <- max(1, which(par_discretization < min(interval[,1])))
    #adding 1 in the max here is necessary for intervals starting at 0 like in betti0, like length of (par_dicretization).
  }
  maximum_value <- min(length(par_discretization), min(which(par_discretization>max(interval[,2]))))
  for (i in minimum_value:maximum_value){
    intervals_containing <- which((!(interval[,1]>par_discretization[i])) & par_discretization[i]<interval[,2])
    if (dimstudy == 0) {
      #print(paste("par", par_discretization[length(par_discretization)]))
      lifespanprofile[i] <- sum(lifespan_vector[intervals_containing]) - par_discretization[length(par_discretization)]
      lifespanprofile[length(lifespanprofile)] = 0
    }
    else {
      lifespanprofile[i] <- sum(lifespan_vector[intervals_containing])
    }
    
  }
  return(lifespanprofile)
}

landscape_from_persist <- function(interval, par_discretization=seq(from=0, to=2, by=0.01), dimstudy, Kmax=4) {
  if (ncol(interval) == 3) {
    IntervalMatrix <- interval[interval[,1]==dimstudy, 2:3]
  } else{
    IntervalMatrix <- interval
  }
  
  persistencelandscapematrix <- matrix(0,nrow=Kmax,ncol=length(par_discretization))
  
  if (!(length(IntervalMatrix) == 0) && length(dim(IntervalMatrix)) > 1) {
    minimum_value <- max(1, which(par_discretization<min(IntervalMatrix[,1])))
    maximum_value <- min(length(par_discretization), min(which(par_discretization > max(IntervalMatrix[,2]))))
    #print(minimum_value)
    #print(maximum_value)
    for (i in minimum_value:maximum_value){
      intervals_containing <- which((!(IntervalMatrix[,1]>par_discretization[i])) & par_discretization[i]<IntervalMatrix[,2])
      num_intervalos <- length(intervals_containing)
      #print(num_intervalos)
      if(num_intervalos>0){
        #print("hi")
        orderedptsinlandscapes <- sort(
          pmin(
            par_discretization[i]-IntervalMatrix[intervals_containing,1], # pmin vectorial min(t-b,d-t)
            IntervalMatrix[intervals_containing,2]-par_discretization[i]),
          decreasing=TRUE)
        persistencelandscapematrix[1:min(num_intervalos,Kmax),i] <- orderedptsinlandscapes[1:min(num_intervalos,Kmax)]
      }
    }  
  }
  
  if (length(dim(IntervalMatrix)) == 1) {
    minimum_value <- max(1,which(par_discretization<min(IntervalMatrix[1])))
    maximum_value <- min(length(par_discretization),min(which(par_discretization>max(IntervalMatrix[2]))))
    for (i in minimum_value:maximum_value){
      intervals_containing <- which((!(IntervalMatrix[1]>par_discretization[i])) & par_discretization[i]<IntervalMatrix[2])
      num_intervalos <- length(intervals_containing)
      if(num_intervalos>0){
        orderedptsinlandscapes <- sort(
          pmin(
            par_discretization[i]-IntervalMatrix[1], # pmin vectorial min(t-b,d-t)
            IntervalMatrix[2]-par_discretization[i]),
          decreasing=TRUE)
        persistencelandscapematrix[1:min(num_intervalos,Kmax),i] <- orderedptsinlandscapes[1:min(num_intervalos,Kmax)]
      }
    }  
  }
  
  return(persistencelandscapematrix)
}

rank_from_persist <- function(interval, par_discretization=seq(from=0,to=2,by=0.01), dimstudy) {
  interval = interval[interval[,1]==dimstudy, 2:3, drop = FALSE]
  #print(dimstudy)
  if (dimstudy == 0) {
    #print("hey")
    rank_func_profile = rep(0, length(par_discretization))
    #print(rank_func_profile)
    for (i in 1:length(par_discretization)) {
      intervals_containing <- which(interval[,2] >= par_discretization[i])
      rank_func_profile [i] = length(intervals_containing)
    }
    return(rank_func_profile)
  } else  {
    #print("hi")
    rankfuncmatrix = matrix(data=0, 
                            nrow=length(par_discretization),
                            ncol=length(par_discretization))
    if (!is.null(interval) & nrow(interval) > 0) {
      for (i in 1:length(par_discretization)) {
        for (j in 1:length(par_discretization)) {
          intervals_containing <- which((interval[,1] <= par_discretization[i]) & interval[,2] >= par_discretization[j])
          rankfuncmatrix[i,j] = length(intervals_containing)
        }
      }
    }
    return(rankfuncmatrix)
  }
 

  
}



