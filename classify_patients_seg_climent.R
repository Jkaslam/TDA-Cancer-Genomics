### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)

source("TDAsignal.R")
source("curves_from_persist.R")
source("statisticsofcurves.R")

maxfiltrationvalue=2
par_discretization = seq(from=0.01,to=maxfiltrationvalue,by=0.01)
dimofstudy = 0
max_k = 4

### READING THE DATA

climentdata <- read.delim("./DATA/climent_segments_data_full.txt")
climent_phen <- read.delim("./DATA/climent_segments_phen.txt")


phenotype = "ER"
curvetypes = c("betti", "lifespan", "landscape2", "landscape3", "landscape4")

for (curve in curvetypes) {
  #significant_segments <- read.csv(paste("./DATA/HER2 ", curve, " Non-Int Segments.csv"))
  significant_segments <- read.csv("./DATA/ER landscape4 Non-Int Segments.csv")
  num_segments = nrow(significant_segments)
  if (phenotype == "HER2")  {
    phen_patients = which(climent_phen$HER2==1)
  }
  else {
    if (phenotype == "ER") {
      phen_patients = which(climent_phen$ER == 1)
    }
    else {
      phen_patients = which(climent_phen$Molecular_Subtypes == phenotype)
    }
  }
  
  num_phen_patients = length(phen_patients)
  numpatients = ncol(climentdata)-6
  patientlist = 1:numpatients
  controls = patientlist[!(patientlist%in%phen_patients)]
  
  segment_classifications = data.frame(matrix(
    nrow=num_segments,
    ncol=numpatients + 1))
  # Loops through significant segments for the phenotype
  for (seg in 1:num_segments) {
    chromosomeofstudy = significant_segments[seg, 1]
    armofstudy = significant_segments[seg, 2]
    segmentofstudy = significant_segments[seg, 3]
    
    persistence_path = "/Users/jkaslam/Desktop/ClimentPersistence/Dim0/"
    persistence = read.csv(paste(persistence_path, chromosomeofstudy, armofstudy, segmentofstudy,".csv",sep=""))
    persistence = persistence[,2:4]
    delimiting_row_indices = which(rowSums(persistence) == 0.0)
    
    intervals = vector(mode = "list", nrow(persistence))
    int_row_to_fill = 1
    for (i in 1:(length(delimiting_row_indices) - 1)) {
      delim1 = delimiting_row_indices[[i]]
      delim2 = delimiting_row_indices[[i+1]]
      
      if (delim2 - delim1 > 1) {
        intervals[[int_row_to_fill]] = persistence[(delim1 + 1):(delim2 - 1),]
        
      }
      else {
        intervals[[int_row_to_fill]] = matrix(data=0, nrow = 1, ncol = 3)
      }
      int_row_to_fill = int_row_to_fill + 1
    }
    if (delimiting_row_indices[length(delimiting_row_indices)] != nrow(persistence)) {
      intervals[[int_row_to_fill]] = persistence[(delimiting_row_indices[length(delimiting_row_indices)]+1):nrow(persistence),]
    }
    else {
      intervals[[int_row_to_fill]] = intervals[[int_row_to_fill]] = matrix(data=0, nrow = 1, ncol = 3)
    }
    
    intervals = intervals[!sapply(intervals, is.null)]
    
    if (substr(curve, 1, 4) == "land") {
      landscape_profiles <- lapply(intervals, landscape_from_persist, par_discretization = par_discretization, dimstudy = dimofstudy, Kmax = max_k)
      landscapes = matrix(0,nrow=numpatients,ncol=length(par_discretization))
      landtype = strtoi(substr(curve, nchar(curve), nchar(curve)))
      for (k in 1:length(landscape_profiles)) {
        landscapes[k, ] = landscape_profiles[[k]][landtype,]
      }
      curvematrix = landscapes
    }
    if (curve == "betti") {
      betti_profiles <-lapply(intervals, betti_from_persist,
                              par_discretization=par_discretization, dimstudy = dimofstudy)
      bettiprofilematrix <- matrix(data=0, 
                                   nrow=numpatients,
                                   ncol=length(par_discretization))
      for (j in 1:numpatients){
        bettiprofilematrix[j,] = betti_profiles[[j]]
      }
      curvematrix = bettiprofilematrix
    }
    if (curve == "lifespan") {
      lifespan_profiles <- lapply(intervals, lifespan_from_persist,
                                  par_discretization=par_discretization, dimstudy = dimofstudy)
      lifespanprofilematrix <- matrix(data=0, 
                                      nrow=numpatients,
                                      ncol=length(par_discretization))
      for (j in 1:numpatients) {
        lifespanprofilematrix[j, ] = lifespan_profiles[[j]]
      }
      curvematrix = lifespanprofilematrix
    }
    curr_classifications = rep(0, numpatients)
    for (patient in 1:numpatients) {
      pat_curve = curvematrix[patient,]
      test_curve_no_pat = colMeans(curvematrix[phen_patients[!phen_patients%in%c(patient)],])
      ctrl_curve_no_pat = colMeans(curvematrix[controls[!controls%in%c(patient)],])
      test_stat = sum((pat_curve - test_curve_no_pat)^2)
      ctrl_stat = sum((pat_curve - ctrl_curve_no_pat)^2)
      
      if (test_stat < ctrl_stat) {
        curr_classifications[patient] = 1
      }
    }
    segment_classifications[seg, ] = c(paste("seg", chromosomeofstudy, armofstudy, segmentofstudy, sep=""), curr_classifications)
  }
  colnames(segment_classifications) = c("Segments", colnames(climentdata[,7:ncol(climentdata)]))
  rownames(segment_classifications) = segment_classifications[,1]
  segment_classifications = segment_classifications[,2:ncol(segment_classifications)]
  normal_like = c("X187", "X219", "X256", "X370")
  segment_classifications = segment_classifications[,!(colnames(segment_classifications)%in%normal_like)]
  # Change the path based on where you want the file to end up
  path = paste("./DATA/LogisticData/", phenotype, curve, "climentsegclass.csv", sep="")
  write.csv(segment_classifications, path)
  
}