# This script will classify a patient as significant for the center of mass 
# if the average for the arm (CM of the patient) lays outside de
# confidence interval from the control group of the phenotype

# The confidence interval is from a t distribution: 
# avg +/- t*sd/sqrt(n)   t=alfa/2 percentile from a t with n-1 d.f.

# INPUT
# size, mean and standard deviation for the control group from
# the file created by 9_mean_diff_perm_NoOut
# TODO: READ SIZE, MEAN AND SD FROM FILE INSTEAD AS ARGUMENT

# OUTPUT
# a new variable will be create in the phenotype file SET_phen.txt


# Here is the info for Horlings ER+, n_ctrl=28
#  16p mean_ctrl=0.02, sd=0.08
#  16q mean_ctrl= 0, sd=0.11

cmdata <- read.delim("/Users/jkaslam/ArsuagaLab/Results/horlings/CenterMass/horlings_mean_diff_perm_noOut_HER2_seed1sig.txt")

dataSet <- "horlings";
#dataSet <- "climent"
phenotype <- "HER2";
sig <- 0.05

# Working directory
begPath <- "~/Research";
#begPath <- "..";

# The following CGH_start is minus 1
CGH_start <- 7;


### INPUT/OUTPUT FILES####

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table("/Users/jkaslam/ArsuagaLab/TAaCGH/Research/Data/horlings/horlings_phen.txt", header=TRUE, sep="\t");
#phenData <- read.table("/Users/jkaslam/ArsuagaLab/TAaCGH/Research/Data/climent/climent_phen.txt", header=TRUE, sep="\t");

# Get the CGH data
dataPath = "./DATA/horlings_update_data_full.txt"
#dataPath = "./DATA/climent_segments_data_full.txt"
data <- read.table(dataPath, header=T, sep='\t', comment.char='');

cm_classifications = data.frame(matrix(nrow=nrow(cmdata), ncol=ncol(data) - 5 + 1))

for (sigarm in 1:nrow(cmdata)) {
  chrom <- cmdata$Chr[sigarm]
  arm <- cmdata$Arm[sigarm]
  n_ctrl <- cmdata$ControlObsNum[sigarm]
  mean_ctrl <- cmdata$mean_ctrl[sigarm]
  sd_ctrl <- cmdata$sd_ctrl[sigarm]
  # type: gain or del
  print(cmdata$Arm_conclusion[sigarm])
  if (startsWith(cmdata$Arm_conclusion[sigarm], "G")) {
    type <- "gain"
  }
  if (startsWith(cmdata$Arm_conclusion[sigarm], "L")) {
    type <- "del"
  } else{
    type = ""
  }
  #print(type)
  
  ######### BEGIN
  
  # Subset and get only the data for the arm
  dataIndices <- intersect(which(data$Chrom == chrom), which(data$Arm == arm));
  arm_data <- data[dataIndices,CGH_start:ncol(data)];
  
  # CM: center of mass by patient
  CM_i <- colMeans(arm_data)
  
  cmclassification <- rep(0, ncol(data) - 5)
  
  if (type=="gain") {
    upper <- mean_ctrl + qt(sig,df=n_ctrl-1, lower.tail=FALSE)*sd_ctrl/sqrt(n_ctrl)
    cmclassification[CM_i>upper] <- 1
  }else if (type=="del") {
    lower <- mean_ctrl - qt(sig,df=n_ctrl-1)*sd_ctrl/sqrt(n_ctrl)
    cmclassification[CM_i<lower] <- 1
    #print(CM_i<lower)
  } else { #print("incorrect input for argument type {gain, del}")}
    
  }
  cm_classifications[sigarm, ] = c(paste("arm", chrom, arm, sep=""), cmclassification)
}

colnames(cm_classifications) = c("Arm", colnames(data[,6:ncol(data)]))
rownames(cm_classifications) = cm_classifications[,1]
cm_classifications = cm_classifications[, 2:ncol(cm_classifications)]
normal_like = c("X187", "X219", "X256", "X370")
#cm_classifications = subset(cm_classifications, select = -c("X187", "X219", "X256", "X370"))
cm_classifications = cm_classifications[,!(colnames(cm_classifications)%in%normal_like)]

path = paste("./DATA/LogisticData/", phenotype, "armclass.csv", sep="")
write.csv(cm_classifications, path)