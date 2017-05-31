############################
# TEDDY tools: Count Manipulation
# Author: Ricardo Ramirez
# 
# This collection of functions have the objective of transforming a count Matrix in different formats
# - This functions use filtered and normalised data (!!!)
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################

# You get the median or mean value of each feature for every group defined
# If you want to get the mean/median value of a feature at each time point for cases and controls:
# -Run this functions for each time-point
# -Create count and target objects for each group and use the function with time as group column
getGroupMedian = function(groupcolumn, stargets, sexprmat){
  #Inputs:
  #groupcolumn: attribute to be used for the separation
  #stargets: set of samples to be reduced
  #sexprmat: expression matrix related to selected annotations
  meanEX = c()
  colLabs = c()
  for(f in sort(unique(stargets[,groupcolumn]))){
    colLabs = c(colLabs,f)
    fsamples = sexprmat[,as.character(stargets$sample_mask_id[stargets[,groupcolumn]==f])]
    meanEX = cbind(meanEX, rowMedians(fsamples))
  }
  colLabs = as.character(colLabs)
  rownames(meanEX) = rownames(sexprmat)
  colnames(meanEX) = c(paste(groupcolumn,colLabs,sep = "_"))
  return(meanEX)
}

getGroupMean = function(groupcolumn, stargets, sexprmat){
  #Inputs:
  #groupcolumn: attribute to be used for the separation
  #stargets: set of samples to be reduced
  #sexprmat: expression matrix related to selected annotations
  meanEX = c()
  colLabs = c()
  for(f in unique(stargets[,groupcolumn])){
    colLabs = c(colLabs,f)
    fsamples = sexprmat[,as.character(stargets$sample_mask_id[stargets[,groupcolumn]==f])]
    meanEX = cbind(meanEX, rowMeans(fsamples))
  }
  meanEX = cbind(rownames(meanEX),meanEX)
  colnames(meanEX) = c("feature_name",paste(groupcolumn,colLabs,sep = "_"))
  return(meanEX)
}

# This function calculates the difference of the mean values of the time points used for teddy for the 4 possible periods
# timepoints = 12,9,6,3,0
# periods = 12-9/ 9-6 / 6-3 / 3-0
get_ChangeMatrix = function(CountMat,AnnotationMat){
  #
  # Calculates the mean expression change from one time point to its previous one, for each gene
  #
  
  
  meanMatrix = getGroupMean(groupcolumn = "time", stargets = AnnotationMat, sexprmat = CountMat)[,-1]
  rn = rownames(meanMatrix)
  meanMatrix = apply(meanMatrix,2,as.numeric)
  rownames(meanMatrix) = rn
  
  ChangeMatrix = matrix(0,nrow = length(rn),ncol = 4)
  rownames(ChangeMatrix) = rn
  colnames(ChangeMatrix) = c("3-0","6-3","9-6","12-9")
  
  ChangeMatrix[,"3-0"] = meanMatrix[,"time_3"] - meanMatrix[,"time_0"]
  ChangeMatrix[,"6-3"] = meanMatrix[,"time_6"] - meanMatrix[,"time_3"]
  ChangeMatrix[,"9-6"] = meanMatrix[,"time_9"] - meanMatrix[,"time_6"]
  ChangeMatrix[,"12-9"] =meanMatrix[,"time_12"] - meanMatrix[,"time_9"]
  
  return(ChangeMatrix)
  
}

# This function creates a matrix that merges time point measurements in a single vector
# so for each patient you'll have all of his data
meltCountMatrix = function(targetDF,countMat,selFeat,meltType = "both"){
  
  meltMat = c()
  
  if(meltType=="case"){
    targetDF = filter(targetDF,outcome==1)
  }else if(meltType=="control"){
    targetDF = filter(targetDF,outcome==0)
  }
  
  patients = unique(targetDF$mask_id)
  
  for(patient in patients){
    
    patient_row = c()
    patient_samples = filter(targetDF,mask_id == patient)
    
    for(t in unique(sort(targetDF$time))){
      
      if(t %in% patient_samples$time){
        
        sample_id = (patient_samples[patient_samples$time==t,"sample_mask_id"])[1]
        expr_vec = countMat[selFeat,as.character(sample_id)]
        names(expr_vec) = paste(selFeat,"_t",as.character(t),sep="")
        patient_row = c(patient_row, expr_vec)
        
      } else{
        
        expr_vec = rep(NA, length(selFeat))
        names(expr_vec) = paste(selFeat,"_t",as.character(t),sep="")
        patient_row = c(patient_row, expr_vec)
      }
      
    }
    
    meltMat = rbind(meltMat,patient_row)
    
  }
  
  rownames(meltMat) = patients
  return(meltMat)
}