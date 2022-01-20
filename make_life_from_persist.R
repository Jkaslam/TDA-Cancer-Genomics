### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)

source("TDAsignal.R")
source("curves_from_persist.R")
source("statisticsofcurves.R")

maxfiltrationvalue=.5
filtrationvector = seq(from=0.01,to=maxfiltrationvalue,by=0.01)
dimofstudy = 0

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
  pvaluetable <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(pvaluetable) <- c("Test", "Control", "Chrom", "Arm", "Segment", "L1 P-Val", "L2 P-Val", "bp Range", "Cytoband Range", "FDR L1", "FDR L2")
  
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
        # print(segmentlengths[segmentindex])
        conditionChromosome <- horlings_dict$Chrom==chromosomeofstudy & horlings_dict$Arm==armofstudy
        conditionSegment <- horlings_dict$Segment %in% segmentofstudy
        
        horlings_dictbis <- horlings_dict[conditionChromosome & conditionSegment,]
        #print(horlings_dictbis)
        
        seqstart <- min(horlings_dictbis$Beg)+1 #positions are 0-based and R vectors are 1-based.
        seqend <- max(horlings_dictbis$End)+1
        
        #EXTRACT RELEVANT DATA. PATIENTS START AT COLUMN 6.
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
        
        persistence_path = "/Users/jkaslam/Desktop/Horlings_Persistence/Dim0/"
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
        par_discretization=seq(from=0, to=maxfiltrationvalue, by=0.01)
        
        life_profiles <-lapply(intervals, lifespan_from_persist,
                                par_discretization=par_discretization, dimstudy = dimofstudy)
        
        lifeprofilematrix <- matrix(data=0, 
                                     nrow=numpatients,
                                     ncol=length(par_discretization))
        
        for (i in 1:numpatients){
          lifeprofilematrix[i,] = life_profiles[[i]]
        }
        
        curvetype <- "life"
        curvematrix<-lifeprofilematrix
        statofcurvesoutput = statisticsofcurves(curvematrix, curvetype, cases, controls)
        l1andl2pvals = c(statofcurvesoutput[[1]], statofcurvesoutput[[2]])
        
        # COMMENT THE SUCCEEDING BLOCK OF CODE OUT IF YOU DON'T WANT TO PLOT BOTH THE TEST AND CONTROL CURVES
        # Plot the average test and control curves to determine whether the test curve is on top or not
        # allowing us to eliminate regions detected as significant where the control set caused the
        # detection of significance.
        #HpcStor/home/jkaslam
        plotpath = paste("/Users/jkaslam/Desktop/P-Vals-L2/lifespan0-plots/", typeName, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
        #plotpath = paste("/HpcStor/home/jkaslam/TDA/betti0-plots/", typeName, chromosomeofstudy, armofstudy, segmentofstudy, ".jpeg", sep="")
        jpeg(file = plotpath)
        #
        par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, mar = c(5.8, 4.4, 4.1, 2.1))
        #Test in red
        plot(par_discretization, statofcurvesoutput[[3]], type="l", lwd = 3, col="red", ylim=c(0, 3), xlab = "Filtration Parameter", ylab = "0-Dim Lifespan")
        # Control in blue
        lines(par_discretization, statofcurvesoutput[[4]], type="l", lwd = 3, col="blue")
        legend("topright", legend=c("Test", "Control"), col=c("red", "blue"), lty = 1:1, lwd = 3)
        title(paste("Lifespan 0 Curve", paste(chromosomeofstudy, armofstudy, "s", segmentofstudy, sep = ""), sep = " "))
        dev.off()
        
        # Retrieves the base pair range for the current chromosome arm segment
        currchromarmsegcond <- horlings_dict$Chrom==chromosomeofstudy & horlings_dict$Arm==armofstudy & horlings_dict$Segment==segmentofstudy
        currbpstart = horlings_dict[currchromarmsegcond, 7]
        currbpend = horlings_dict[currchromarmsegcond, 8]
        currbprange = paste(currbpstart, currbpend,  sep ="-")
        
        # Retrieves the cytoband range for the current chromosome arm segment
        currcytostart = horlings_dict[currchromarmsegcond, 9]
        currcytoend = horlings_dict[currchromarmsegcond, 10]
        currcytorange = paste(currcytostart, currcytoend,  sep ="-")
        
        pvaluetable[nrow(pvaluetable) +1, ] = list(typeName, "Rest", chromosomeofstudy, armofstudy, segmentofstudy, l1andl2pvals[1], l1andl2pvals[2], currbprange, currcytorange, 0, 0)
      }
      
      j = j + 1
      
    }
    
  }
  
  # Apply FDR to the p-values
  pvaluetable[, 10] = p.adjust(pvaluetable[, 6], method="fdr")
  pvaluetable[, 11] = p.adjust(pvaluetable[, 7], method="fdr")
  
  # Change the path based on where you want the file to end up
  path = paste("/Users/jkaslam/Desktop/P-Vals-L2/lifespan0-pvals/", typeName, "life0.csv", sep="")
  write.csv(pvaluetable, path)
}