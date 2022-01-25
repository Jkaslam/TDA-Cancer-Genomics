# Computes the p-values for each breast cancer subtype and each segment of
# the chromosome using landscape curves. The output is saved in a CSV file. 
# The plots of the average landscape curves for the test and control of each 
# segment are plotted as well. These plots are saved in JPEG files.

### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)

source("TDAsignal.R")
source("curves_from_persist.R")
source("statisticsofcurves.R")

maxfiltrationvalue=.5
dimofstudy = 0
# Maximum number of persistence landscapes to compute
max_k = 4
# The path where the user wants p-value CSV files to go.
pval_path = ""
# The path where the user wishes landscape plots to go.
landscape_plot_path = ""

### READING THE DATA
horlingsdata <- read.delim("./DATA/horlings_update_data_full.txt")
horlings_phen <- read.delim("./DATA/horlings_segments_phen.txt")
horlings_dict <- read.delim("./DATA/horlings_sect_dict_cyto.txt")

numpatients = ncol(horlingsdata)-5

# The segment lengths of the arms of each chromosome in order
segmentlengths <- c(11, 6, 7, 12, 8, 10, 3, 12, 3, 11, 4, 9, 3, 9, 3, 8, 4, 8, 3, 8, 4, 8, 2, 9, 8, 8, 6, 2, 4, 2, 6, 4, 2, 2, 3, 1, 3, 3, 7)
nopchromosomoes <- c(13, 14, 15, 18, 20, 21, 22)

# Loops through the subtypes we care about
typeNames = c("Luminal A", "Luminal B", "Basal", "HER2+", "ER")
for (ind in 5) {
  pvaluetables = vector(mode = "list", max_k)
  for (i in 1:max_k) {
    pvaluetable <- data.frame(matrix(ncol = 11, nrow = 0))
    colnames(pvaluetable) <- c("Test", "Control", "Chrom", "Arm", "Segment", "L1 P-Val", "L2 P-Val", "bp Range", "Cytoband Range", "FDR L1", "FDR L2")
    pvaluetables[[i]] = pvaluetable
  }
  
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
        
        # Positions are 0-based and R vectors are 1-based.
        seqstart <- min(horlings_dictbis$Beg)+1 
        seqend <- max(horlings_dictbis$End)+1
        
        # EXTRACT RELEVANT DATA. PATIENTS START AT COLUMN 6.
        patients_data <- horlingsdata[seqstart:seqend,6:ncol(horlingsdata)]
        
        probes_data <- horlingsdata[seqstart:seqend,1:4] #all the chromosome info for this segment.
        numpatients <- ncol(patients_data)
        numprobes <- nrow(patients_data)
        colnames(patients_data) <- paste0("patient_",1:numpatients)
        
        ####################################
        ## SELECT CASES AND CONTROLS
        
        typeName= typeNames[ind]
        if (typeName == "HER2+") {
          cases=which(horlings_phen$HER2==1)
        }
        else {
          cases = which(horlings_phen$Molecular_Subtypes == typeName)
        }
        
        patientlist = 1:numpatients
        cases = which(horlings_phen$ER == 1)
        controls=patientlist[!(patientlist%in%cases)]
        if (typeName == "Luminal B") {
          basal = which(horlings_phen$Molecular_Subtypes == "Basal")
          both = c(cases, basal)
          controls=patientlist[!(patientlist%in%both)]
        }
        #################################
        
        # The path where precomputed Horlings persistence is located
        persistence_path = ""
        persistence = read.csv(paste(persistence_path, chromosomeofstudy, armofstudy, segmentofstudy,".csv",sep=""))
        persistence = persistence[,2:4]
        delimiting_row_indices = which(rowSums(persistence) == 0.0)
        
        # Compute the list of interval matrices from the Horlings persistence
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
        par_discretization=seq(from=0, to=maxfiltrationvalue, by=0.01)
        landscape_profiles <- lapply(intervals, landscape_from_persist, par_discretization = par_discretization, dimstudy = dimofstudy, Kmax = max_k)
        
        # Retrieves the base pair range for the current chromosome arm segment
        currchromarmsegcond <- horlings_dict$Chrom==chromosomeofstudy & horlings_dict$Arm==armofstudy & horlings_dict$Segment==segmentofstudy
        currbpstart = horlings_dict[currchromarmsegcond, 7]
        currbpend = horlings_dict[currchromarmsegcond, 8]
        currbprange = paste(currbpstart, currbpend,  sep ="-")
        
        # Retrieves the cytoband range for the current chromosome arm segment
        currcytostart = horlings_dict[currchromarmsegcond, 9]
        currcytoend = horlings_dict[currchromarmsegcond, 10]
        currcytorange = paste(currcytostart, currcytoend,  sep ="-")
        
        for (i in 1:max_k) {
          landscapes = matrix(0,nrow=numpatients,ncol=length(par_discretization))
          for (k in 1:length(landscape_profiles)) {
            landscapes[k, ] = landscape_profiles[[k]][i,]
          }
          
          statofcurvesoutput = statisticsofcurves(landscapes, "land", cases, controls)
          l1andl2pvals = c(statofcurvesoutput[[1]], statofcurvesoutput[[2]])
          
          par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, mar = c(5.1, 4.4, 4.1, 2.1))
          plotpath = paste(landscape_plot_path, typeName, chromosomeofstudy, armofstudy, segmentofstudy, "k", i, ".jpeg", sep="")
          jpeg(file = plotpath)
          #Test in red
          plot(par_discretization, statofcurvesoutput[[3]], type="l", col="red", lwd = 3, xlab = "Filtration Parameter", ylab = "", ylim = c(0, .05))
          
          # Control in blue
          lines(par_discretization, statofcurvesoutput[[4]], type="l",col="blue", lwd = 3)
          legend("topright", legend=c("Test", "Control"), col=c("red", "blue"), lty = 1:1, lwd = 3)
          title(paste("Landscape K =", i, paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep=""), segmentofstudy, sep = " "))
      
          dev.off()
          
          num_rows = nrow(pvaluetables[[i]])
          pvaluetables[[i]][num_rows + 1, ] = list(typeName, "Rest", chromosomeofstudy, armofstudy, segmentofstudy, l1andl2pvals[1], l1andl2pvals[2], currbprange, currcytorange, 0, 0)
        }
        
      }
      
      j = j + 1
      
    }
    
  }
}
for (curr_landscape in 1:max_k)
{
  #Apply FDR to the p-values
  pvaluetables[[curr_landscape]][, 10] = p.adjust(pvaluetables[[curr_landscape]][, 6], method="fdr")
  pvaluetables[[curr_landscape]][, 11] = p.adjust(pvaluetables[[curr_landscape]][, 7], method="fdr")
  
  path = paste(pval_path, typeName, "landscapek", curr_landscape, "dim0.csv", sep="")
  write.csv(pvaluetables[[curr_landscape]], path)
}