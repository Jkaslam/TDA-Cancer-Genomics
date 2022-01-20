phenotype = "Basal"
curvetypes = c("betti", "lifespan", "landscape2", "landscape3", "landscape4")
for (curve in curvetypes) {
  # Read in the model data
  horlings_phen = read.delim("./DATA/horlings_phen.txt")
  armfilename = paste("./DATA/LogisticData/", phenotype, "armclass.csv", sep="")
  if (file.exists(armfilename)) {
    arm_data = read.csv(armfilename, header = TRUE, row.names = 1)
  }
  phenotype_pos = t(horlings_phen[phenotype])
  rownames(phenotype_pos) = c(paste(phenotype, "pos", sep=""))
  colnames(phenotype_pos) = colnames(arm_data)
  segfilename = paste("./DATA/LogisticData/", phenotype, curve, "segclass.csv", sep="")
  if (file.exists(segfilename)) {
    seg_data = read.csv(segfilename, header = TRUE, row.names = 1)
    
  } 
  if (file.exists(armfilename) & file.exists(segfilename)) {
    model_data = t(rbind(arm_data, seg_data, phenotype_pos))
  } else {
    if (file.exists(armfilename)) {
      model_data = t(rbind(arm_data, phenotype_pos))
    }
    else {
      model_data = t(rbind(seg_data, phenotype_pos))
    }
  }

  model_data = as.data.frame(model_data)
  
  # # Read in the test data
  # climent_phen = read.delim("./DATA/climent_phen.txt")
  # climent_arm_data = read.csv(paste("./DATA/LogisticData/", phenotype, "climentarmclass.csv", sep=""), header = TRUE, row.names = 1)
  # climent_seg_data = read.csv(paste("./DATA/LogisticData/", phenotype, "landscape4climentsegclass.csv", sep=""), header = TRUE, row.names = 1)


  # full
  f = paste(paste(phenotype, "pos", sep=""), "~", paste(colnames(model_data)[1:ncol(model_data)-1], collapse=" + "))
  logit = glm(formula = f, data = model_data, family = binomial)
  #summary(logit)

  # step forward BIC
  nothing <- glm(paste(phenotype, "pos", " ~ 1", sep=""), data=model_data, family=binomial)
  full = glm(formula = f, data = model_data, family = binomial)
  logit = step(nothing,
       scope=list(lower=formula(nothing),upper=formula(full)),data = model_data, direction="forward", k=log(64))
  print(logit)

  # step forward AIC
  # nothing <- glm(paste(phenotype, "pos", " ~ 1", sep=""), data=model_data, family=binomial)
  # full = glm(formula = f, data = model_data, family = binomial)
  # step(nothing,
  #      scope=list(lower=formula(nothing),upper=formula(full)),data = model_data, direction="forward", k=2)
  # summary(full)

  # Predictions on the model-building set
  model_data$y_hat <- NA
  model_data$y_hat[ predict(logit, newdata=model_data, type = "response")>0.5] <- 1
  model_data$y_hat[ predict(logit, newdata=model_data, type = "response")<=0.5] <- 0

  phenopos = paste(phenotype, "pos", sep="")
  accuracy = length(which(model_data[,phenopos] == model_data$y_hat))/nrow(model_data)
  print(nrow(model_data))
  truepositives = length(which(model_data[,phenopos] == 1 & model_data[,phenopos] == model_data$y_hat))
  truenegatives = length(which(model_data[,phenopos] == 0 & model_data[,phenopos] == model_data$y_hat))
  falsepositives = length(which(model_data[,phenopos] == 0 & model_data$y_hat == 1))
  falsenegatives = length(which(model_data[,phenopos] == 1 & model_data$y_hat == 0))
  print(paste("True negatives", phenotype, curve, truenegatives, sep=" "))
  print(paste("True positives", phenotype, curve,  truepositives, sep=" "))
  print(paste("False negatives", phenotype, curve,  falsenegatives, sep=" "))
  print(paste("False positives", phenotype, curve, falsepositives, sep=" "))
  print(paste("Accuracy", phenotype, curve, accuracy, sep=" "))
}
# # Predictions on the test set
# climent_phenotype_pos = t(climent_phen[phenotype])
# rownames(climent_phenotype_pos) = c(paste(phenotype, "pos", sep=""))
# colnames(climent_phenotype_pos) = colnames(climent_arm_data)

# climent_model_data = t(rbind(climent_arm_data, climent_seg_data, climent_phenotype_pos))
# climent_model_data = as.data.frame(climent_model_data)
# climent_model_data$y_hat <- NA
# climent_model_data$y_hat[ predict(logit, newdata=climent_model_data, type = "response")>0.5] <- 1
# climent_model_data$y_hat[ predict(logit, newdata=climent_model_data, type = "response")<=0.5] <- 0
# 
# climent_accuracy = sum(climent_model_data$ERpos == climent_model_data$y_hat)/nrow(climent_model_data)
# climent_truepositives = length(which(climent_model_data$ERpos == 1 & climent_model_data$ERpos == climent_model_data$y_hat))
# climent_truenegatives = length(which(climent_model_data$ERpos == 0 & climent_model_data$ERpos == climent_model_data$y_hat))
# print(paste("Climent True negatives", climent_truenegatives, sep=" "))
# print(paste("Climent True positives", climent_truepositives, sep=" "))
# print(paste("Climent Accuracy", climent_accuracy, sep=" "))