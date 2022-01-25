# Plots the persistence curves (Betti, lifespan, landscape) of all patients
# from the Horlings dataset to specified directories.

### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)
library(tidyverse)
library(gplots)

source("TDAsignal.R")
source("curves_from_persist.R")
source("statisticsofcurves.R")

### READING THE DATA
maxfiltrationvalue=2
horlingsdata <- read.delim("./DATA/horlings_update_data_full.txt")
horlings_dict <- read.delim("./DATA/horlings_sect_dict_cyto.txt")

# Paths of where to put plots (user changes)
life_path = ""
betti_path = ""
landscape_path = ""

numpatients = ncol(horlingsdata) - 5

# Choose chromosome of study including arm and segment
chromosomeofstudy = 11
armofstudy = "q"
segmentofstudy = 6

# The dimension of homology to consider
dimofstudy = 0
# The number of persistence landscapes to use
maxk = 4

conditionChromosome <- horlings_dict$Chrom==chromosomeofstudy & horlings_dict$Arm==armofstudy
conditionSegment <- horlings_dict$Segment %in% segmentofstudy

horlings_dictbis <- horlings_dict[conditionChromosome & conditionSegment,]

# Positions are 0-based and R vectors are 1-based.
seqstart <- min(horlings_dictbis$Beg)+1 
seqend <- max(horlings_dictbis$End)+1

#EXTRACT RELEVANT DATA. PATIENTS START AT COLUMN 6.
patients_data <- horlingsdata[seqstart:seqend,6:ncol(horlingsdata)]

# All the chromosome info for this segment.
probes_data <- horlingsdata[seqstart:seqend,1:4] 
numpatients <- ncol(patients_data)
numprobes <- nrow(patients_data)
colnames(patients_data) <- paste0("patient_",1:numpatients)

cases= 1:numpatients

# Path of the precomputed persistence for the dataset. 
persistence_path = ""
persistence = read.csv(paste(persistence_path, chromosomeofstudy, armofstudy, segmentofstudy,".csv",sep=""))
persistence = persistence[,2:4]
delimiting_row_indices = which(rowSums(persistence) == 0.0)

# Turn the persistence information into a list of birth-death intervals
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
} else {
  intervals[[int_row_to_fill]] = intervals[[int_row_to_fill]] = matrix(data=0, nrow = 1, ncol = 3)
}

intervals = intervals[!sapply(intervals, is.null)]
par_discretization=seq(from=0, to=maxfiltrationvalue, by=0.01)

TDAsignaloutput <- apply(patients_data, 2, TDAsignal, maxfiltrationvalue = maxfiltrationvalue, dimstudy = dimofstudy)     

# Compute Betti, lifespan and persistence landscape curves from intervals
betti_profiles<-lapply(intervals, betti_from_persist,
                       par_discretization=par_discretization, dimstudy = dimofstudy)
lifespan_profiles <-lapply(intervals, lifespan_from_persist,
                           par_discretization=par_discretization, dimstudy = dimofstudy)
landscape_profiles <- lapply(intervals, landscape_from_persist, par_discretization = par_discretization, dimstudy = dimofstudy, Kmax = maxk)


bettiprofilematrix <- matrix(data=0, 
                             nrow=numpatients,
                             ncol=length(par_discretization))
lifespanprofilematrix <- bettiprofilematrix



for (i in 1:numpatients){
  bettiprofilematrix[i,] = betti_profiles[[i]]
  lifespanprofilematrix[i,] =lifespan_profiles[[i]]
}

# Plots all persistence curves and saves them to give directories
par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, mar = c(5.1, 4.4, 4.1, 2.1))
for (patient in cases) {
  curr_patient = colnames(horlingsdata)[patient+5]
  for (i in 1:maxk) {
    landscapes = matrix(0,nrow=numpatients,ncol=length(par_discretization))
    for (k in 1:length(landscape_profiles)) {
      landscapes[k, ] = landscape_profiles[[k]][i,]
    }
    
    plotpath = paste(landscape_path, "Landscape", i, "Patient", curr_patient, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
    jpeg(file = plotpath)
    plot(par_discretization, landscapes[patient, ], type="l", col="red", lwd = 3, xlab = "Filtration Parameter", ylab = "Landscape")
    
    title(paste("Patient", curr_patient, paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep=""), "Landscape", paste("K =", i, sep=""), sep=" "))
    dev.off()
  }
  
  if (dimofstudy == 0) {
    plotpath = paste(betti_path, "Betti Patient", curr_patient, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
    jpeg(file = plotpath)
    plot(par_discretization, bettiprofilematrix[patient, ], type="l", col="red", lwd = 3, xlab = "Filtration Parameter", ylab = "Connected Components")
    title(paste("Patient", curr_patient, paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep=""), "Betti-0 Curve", sep=" "))
    dev.off()
  }
  else {
    plotpath = paste(betti_path, "Betti Patient", curr_patient, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
    jpeg(file = plotpath)
    plot(par_discretization, bettiprofilematrix[patient, ], type="l", col="red", lwd = 3, xlab = "Filtration Parameter")
    title(paste("Patient", curr_patient, paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep=""), "Betti-1 Curve",  sep=" "))
    dev.off()
  }
  plotpath = paste(life_path, "Life Patient", curr_patient, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
  jpeg(file = plotpath)
  plot(par_discretization, lifespanprofilematrix[patient, ], type="l", col="red", lwd = 3, xlab = "Filtration Parameter", ylab = "Lifespan")
  title(paste("Patient", curr_patient, paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep=""), paste("Lifespan-", dimofstudy, sep=""), "Curve", sep=" "))
  dev.off()
  
}

