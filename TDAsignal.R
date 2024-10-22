# Takes in time series data and computes the sliding window point cloud from
# this data and returns the persistent homology of this cloud as a list of
# birth-death intervals with corresponding dimensions. 

#This auxiliary function takes as input a single matrix (cycle location input) 
#and returns a vector with the non repeated corresponding set of values 
# (x_i, x_{i+1}) that appears in a signal
vector_generator <- function(x, windowsize, indicesextended) {
  aux <- unique(as.vector(x))
  for (i in 1:(windowsize-1)) {aux <- c(aux, aux+i)}
  locations <-  indicesextended[sort(unique(aux))]
  return(locations)
}

######## The main function TDAsignal, return a list with a possibly empty 
# intervals as first component and a possibly empty list of signal location 
# values for each interval.
TDAsignal <- function(signal, windowsize=2, maxfiltrationvalue=.8, dimstudy=0) {

    numpoints <- length(signal)
    indicesaux <- c((1:numpoints),(1:(windowsize-1)))
    
    # CREATE CLOUD
    # Wraps around the extra values 
    signalaux <- c(signal, signal[1:(windowsize - 1)])
    
    cloud <- matrix(0, nrow = numpoints, ncol = windowsize)
    for (i in 1:numpoints) {
      cloud[i,] <- signalaux[i:(i + windowsize - 1)]
    }
    
  # Distance matrix.
  distancia <- as.matrix(dist(cloud,upper = TRUE))
  
  # CALCULATION OF PERSISTENCE DIAGRAMS + GENERATORS
  outputTDA <- ripsDiag(
    distancia,
    maxdimension = 1,
    maxscale = maxfiltrationvalue,
    dist = "arbitrary",
    library = "Dionysus",
    location = TRUE
  )

  indicesgenerators <- which(outputTDA$diagram[,"dimension"]==dimstudy)
  numgenerators <- length(indicesgenerators)
  intervals <- matrix(outputTDA$diagram[indicesgenerators,2:3], ncol=2)
  dimnames(intervals) <- list(NULL,c("Birth","Death"))
  if (numgenerators>0) {
     spangenerators <- intervals[,"Death"]-intervals[,"Birth"]
     totalspan <- sum(spangenerators)
   }
  else {
     spangenerators <- 0
     totalspan <- 0
  }
  
  generators = vector(mode = "list", 1)
  outputTDAsignal <- list(intervals=intervals, generators = generators, diagram = outputTDA$diagram)
  
  return(outputTDAsignal)  
}
