### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)
library(tidyverse)

source("TDAsignal.R")
source("persistentcurvemaker.R")
source("generator_list_maker.R")
source("statisticsofcurves.R")

# Computes the persistence of the Horlings dataset for each patient. This data
# is saved in CSV files labeled by chromosome arms. 
dimofstudy = 0

### READING THE DATA
horlingsdata <- read.delim("./DATA/horlings_update_data_full.txt")
horlings_phen <- read.delim("./DATA/horlings_segments_phen.txt")
horlings_dict <- read.delim("./DATA/horlings_sect_dict_cyto.txt")

numpatients = ncol(horlingsdata)-5

# The segment lengths of the arms of each chromosome in order
segmentlengths <- c(11, 6, 7, 12, 8, 10, 3, 12, 3, 11, 4, 9, 3, 9, 3, 8, 4, 8, 3, 8, 4, 8, 2, 9, 8, 8, 6, 2, 4, 2, 6, 4, 2, 2, 3, 1, 3, 3, 7)
nopchromosomoes <- c(13, 14, 15, 18, 20, 21, 22)

  segmentindex = 0;
  # Loops through the chromosomes
  for (chromosomeofstudy in 1:23) {
    # Determines the current arm
    j = 1
    while (j <= 2) {
      segmentindex = segmentindex + 1
      if (j == 1 && !(chromosomeofstudy %in% nopchromosomoes)) {
        armofstudy = "p"
      }
      else {
        j = j + 1
        armofstudy = "q"
      }
      
      # Loops through the segments in the current chromosome arm
      for (segmentofstudy in 1:segmentlengths[segmentindex]) {
        conditionChromosome <- horlings_dict$Chrom==chromosomeofstudy & horlings_dict$Arm==armofstudy
        conditionSegment <- horlings_dict$Segment %in% segmentofstudy
        
        horlings_dictbis <- horlings_dict[conditionChromosome & conditionSegment,]
        
        seqstart <- min(horlings_dictbis$Beg)+1 #positions are 0-based and R vectors are 1-based.
        seqend <- max(horlings_dictbis$End)+1
        
        #EXTRACT RELEVANT DATA. PATIENTS START AT COLUMN 6.
        patients_data <- horlingsdata[seqstart:seqend,6:ncol(horlingsdata)]
        colnames(patients_data) <- paste0("patient_",1:numpatients)
        
        
        TDAsignaloutput <- apply(patients_data, 2, TDAsignal)
        diagrams = vector(mode = "list", 2 * length(TDAsignaloutput))
        
        for (i in 1:(length(diagrams)/2)) {
          currTDASignalOutput = TDAsignaloutput[[i]][[3]]
          currTDASignalOutput = currTDASignalOutput[currTDASignalOutput[, 1] == dimofstudy,]
          diagrams[[2*i]] = currTDASignalOutput
          diagrams[[2*i - 1]] = c(0, 0, 0)
        }
        
        all_diags = do.call(rbind, diagrams)
        
        # Set this string to be equal to the path where you want the persistence
        # CSV files. 
        output_path = ""
        persistence_path = paste(output_path, chromosomeofstudy, armofstudy, segmentofstudy,".csv", sep="")
        write.csv(all_diags, persistence_path)
      }
      
      j = j + 1
      
    }
    
  }
  