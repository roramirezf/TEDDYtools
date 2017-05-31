############################
# TEDDY tools: Clustering 
# Author: Ricardo Ramirez
# 
# This collection of functions have the objective of clustering genes or patients, depending on their molecular profiles
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################

#
#
# Get best number of clusters
#
#

myHclust = function(testCounts, DIST = TRUE,Knum=8){
  #Uses 1-correlation, as the metric for HC
  
  if(DIST == TRUE){
    dcorrel = dist(testCounts)
  } else {
    # Uses correlation as distance
    dcorrel = as.dist(matrix(rep(1, nrow(testCounts)^2), nrow(testCounts), nrow(testCounts)) - cor(t(testCounts), use = "pairwise.complete.obs"))
  }
  
  # Get best number of clusters
  asw = numeric(Knum)
  for (k in 2:Knum){
    print(paste("trying",as.character(k),"clusters"))
    asw[[k]] = pam(dcorrel, k)$silinfo$avg.width
  }
  k.best = which.max(asw)
  Gclass = pam(x = dcorrel,k = k.best, cluster.only=TRUE)
  #Gclust = hclust(dcorrel)
  #Gclass = cutree(Gclust,k = k.best) 
  # Get lists of genes
  Glist = list()
  Nclust = unique(Gclass)
  for(i in Nclust){
    Glist[[paste("Clust",as.character(i),sep="_")]] = names(Gclass)[Gclass==i]
  }
  return(Glist)
  
}

get_BNumbClust_PAM = function(testCounts, DIST = TRUE,Knum=8){
  #Uses 1-correlation, as the metric for HC
  
  if(DIST == TRUE){
    dcorrel = dist(testCounts)
  } else {
    # Uses correlation as distance
    #dcorrel = as.dist(matrix(rep(1, nrow(testCounts)^2), nrow(testCounts), nrow(testCounts)) - cor(t(testCounts), use = "pairwise.complete.obs"))
    dcorrel = dist(testCounts)
  }
  
  # Get best number of clusters
  asw = numeric(Knum)
  for (k in 2:Knum){
    print(paste("trying",as.character(k),"clusters"))
    asw[[k]] = pam(dcorrel, k)$silinfo$avg.width
  }
  k.best = which.max(asw)
  Gclass = pam(x = dcorrel,k = k.best, cluster.only=TRUE)
  #Gclust = hclust(dcorrel)
  #Gclass = cutree(Gclust,k = k.best) 
  # Get lists of genes
  Glist = list()
  Nclust = unique(Gclass)
  for(i in Nclust){
    Glist[[paste("Clust",as.character(i),sep="_")]] = names(Gclass)[Gclass==i]
  }
  return(Glist)
  
}

get_BestNumberClusters = function(CountMatrix, method = "SOF"){
  #
  # For a complete (SOF) or incomplete (PAM) matrix
  #
  #
  if(method == "SOF"){ #Leandro's code
    
    mydata <- CountMatrix
    wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
    for (i in 2:16) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
    
    d_clust <- Mclust(as.matrix(mydata), G=1:16)
    Classification<-as.data.frame(d_clust$classification)
    rownames(Classification)<-rownames(mydata)
    colnames(Classification)<-"Classification"
    
    k_names = unique(Classification$Classification)
    
    k_list = list()
    for(k_name in k_names){
      k_label = paste("cluster",as.character(k_name),sep="_")
      k_list[[k_label]] = rownames(Classification)[Classification$Classification==k_name]
    }
    
    return(k_list)
    
  }else if(method == "PAM"){
    
    return(get_BNumbClust_PAM(testCounts = CountMatrix,DIST = FALSE,Knum = 16))
    
  }
}

#
#
# 
#
#

testSignatureDynamics = function(testCounts,testTargets,SignatureDictionary,folderPATH){
  # This function plots the dynamics of the ES of a list of gene sets in a specific data set (testCounts/testTargets)
  # -It calculates an Enrichment Score Matrix using the provided list (SignatureDictionary)
  # -Calculates a GLM per time point and exports the results (in plots)
  #
  #
  LM_globalFMAT = getFunctionMatrix_v2(eMat = testCounts,eTargets = testTargets,Signature_Dict = SignatureDictionary,weighted = TRUE)
  LM_globalFMAT = LM_globalFMAT$ESmat
  
  # Filter times with low number of samples
  
  completeD = names(table(testTargets$time))[table(testTargets$time) >=8]
  LM_globalTargets = filter(testTargets, time %in% completeD)
  LM_globalFMAT = LM_globalFMAT[,as.character(LM_globalTargets$sample_mask_id)]
  
  # Get results from T-test
  masigClusters = getSignatures(cMAT = LM_globalFMAT,targetMAT = LM_globalTargets,Features = NULL,time = sort(unique(LM_globalTargets$time)),g.pval = 0.05,s.pval=0.05,model_type = "glm")
  # Create matrix of results
  resultData = masigClusters$RESULTS
  
  pdf(file = paste(folderPATH,"Dynamics.pdf",sep = "_"),height = 10,width = 12)
  for(k in rownames(LM_globalFMAT)){
    plotES(ES_list = LM_globalFMAT[k,],sigTargets = LM_globalTargets,TRSH = 0.05,main_title = k)
  }
  dev.off()
  
  pdf(file = paste(folderPATH,"LMresults.pdf",sep = "_"),height = 8,width = 8)
  SummaryTables = summaryTestTable(resultData)
  dev.off()
  
  return(list("ESmat"=LM_globalFMAT,"GLM_Res"=resultData))
  
}

##########################################

testDynamics_Exhaustive = function(TestCounts_list, TestAnnotations_list, globalResultFolder, clusterList){
  # This functions is mostly used in the Age Comparison, where a set of clusters is Enriched in a test set
  # The medians of the enrichment score are plotted for every test group
  #
  #
  MedianValues = list()
  
  TestNames = names(TestCounts_list)
  
  Tset_Results = list()
  
  for(Tset in TestNames){
    
    print(paste("Testing",Tset))
    
    countT = TestCounts_list[[Tset]]
    annT = TestAnnotations_list[[Tset]]
    
    
    folderPATH = paste(globalResultFolder,Tset,"/",sep = "")
    system(paste("mkdir ",folderPATH))
    
    #1st part, Make an enrichment matrix for each testset
    #Compare if there's is a difference between cases and controls, using Enrichment Scores
    print(paste("Results at", folderPATH))
    dyn_RES = testSignatureDynamics(testCounts = countT,testTargets = annT,SignatureDictionary = clusterList, folderPATH = folderPATH)
    Tset_Results[[Tset]] = dyn_RES
    
    # Then, I calculated the median of the enrichment score of all the cases from the test set at all times
    testcases = filter(annT, outcome==1, sample_mask_id %in% colnames(dyn_RES$ESmat))
    medianResponse = getGroupMedian(groupcolumn = "time",stargets = testcases,sexprmat = dyn_RES$ESmat[,as.character(testcases$sample_mask_id)])
    
    MedianValues[[Tset]] = medianResponse
    
  }
  
  #Plot all medians together
  pdf(paste(globalResultFolder,"MedianComparison.pdf",sep=""),width = 12, height = 10)
  for(K in names(clusterList)){
    KmedianResps = list()
    for(Tset in TestNames){
      KmedianResps[[Tset]] = MedianValues[[Tset]][K,,drop=F]
    }
    plot_CompareDynamics(MedianResps = KmedianResps,main = K)
  }
  dev.off()
  
  #Plot all combinations of groups
  
  AgeCombinations = combn(TestNames,2)
  
  apply(AgeCombinations,2,function(x){
    Agecomb = paste(x,collapse = "_")
    G1 = x[1]
    G2 = x[2]
    filename = paste(globalResultFolder,Agecomb,".pdf",sep="")
    MedianResp = MedianValues[x]
    
    pdf(filename,width = 12, height = 10)
    
    for(K in names(clusterList)){
      MR = list()
      MR[[G1]] = MedianResp[[G1]][K,,drop=F]
      MR[[G2]] = MedianResp[[G2]][K,,drop=F]
      plot_CompareDynamics(MedianResps = MR,main = K)
    }
    
    dev.off()
  })
  
  return(Tset_Results)
}

write_clusters = function(globalResultFolder, clusterList){
  folderPATH = paste(globalResultFolder,"ClusterList","/",sep = "")
  system(paste("mkdir ",folderPATH))
  for(k in names(clusterList)){
    Elist = (na.omit(data.frame(mapIds(illuminaHumanv4.db, keys= clusterList[[k]], keytype = "SYMBOL", column = "ENTREZID",multiVals = "first"),
                                stringsAsFactors = F)))[,1]
    write.table(clusterList[[k]], file = paste(folderPATH,k,".txt",sep=""), sep = "\t",row.names = F,col.names = F,quote = F)
    write.table(Elist, file = paste(folderPATH,k,"_ENTREZ",".txt",sep=""), sep = "\t",row.names = F,col.names = F,quote = F)
  }
}























