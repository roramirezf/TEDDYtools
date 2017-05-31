############################
# TEDDY tools: Linear Models
# Author: Ricardo Ramirez
# 
# This collection of functions have the objective of making differential expression analyses
# - This functions use filtered and normalised data (!!!)
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################

#
# Time Specific Linear Models
#

# It only compares outcome!
runTEDDYlimma_basic = function(stargets, sEXPMAT){
  #
  # Functions fit a simple Sick vs Healthy contrast
  #
  
  geneList = list()
  
  Sclass = ifelse(stargets$outcome==0,"Healthy","Sick")
  stargets = cbind(stargets,Sclass)
  
  f = factor(stargets$Sclass,levels = sort(unique(stargets$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  fit <- lmFit(sEXPMAT, design)
  cont.IA = makeContrasts(
    Dif = Sick - Healthy,
    levels = design
  )
  fit2 = contrasts.fit(fit,cont.IA)
  fit2 = eBayes(fit2)
  
  lmres = topTable(fit2,coef = 1,adjust.method = "BH",number = Inf)
  
  return(lmres)
}

#Filtering by PVALUE...
filterLimmaGLSv3 = function(GLS,pTRSH=0.05,FC_col = "logFC"){
  
  filt_GLS = GLS[GLS$adj.P.Val<=pTRSH,]
  up_genes = rownames(filt_GLS[filt_GLS[,FC_col]>0,])
  down_genes = rownames(filt_GLS[filt_GLS[,FC_col]<0,])
  
  return(list("up"=up_genes,"down"=down_genes))
  
}

#If you have a list of results...
filterLimmaGLSv3_list = function(GLSlist,TRSH){
  
  geneList = list()
  
  for(contName in names(GLSlist)){
    GLS = GLSlist[[contName]]
    Tgenes = filterLimmaGLSv3(GLS = GLS,pTRSH = TRSH,FC_col = "logFC")
    geneList[[contName]] = Tgenes
  }
  
  return(geneList)
}

#Using a Generalized Linear Model (Not Limma)
fitGLM_v2 = function(cMAT,cTargets){
  # This function is used to fit a linear model at each time point for a given count matrix (it corrects with FDR)
  # Inputs:
  # cMAT=count matrix
  # cTargets=annotation table
  # Outputs:
  # - Table with results
  # - Matrix with the data use for the analysis
  
  lmResults = c()
  samples =  as.character(cTargets$sample_mask_id)
  rownames(cTargets) = samples
  outcome = cTargets[samples,c("outcome","time")]
  filtMAT = data.frame(cbind(outcome,t(cMAT[,samples])))
  filtMAT$outcome = factor(filtMAT$outcome)
  
  for(bf in colnames(filtMAT)[3:(ncol(filtMAT))]){
    cases_values = filtMAT[filtMAT$outcome == 1, bf]
    cases_values = sum(!is.na(cases_values))>1
    
    control_values = filtMAT[filtMAT$outcome == 0, bf]
    control_values = sum(!is.na(control_values))>1
    
    if(cases_values & control_values){
      
      naFREE_FMat = na.omit(filtMAT[,c(bf,"outcome")])
      
      t_lm = summary(lm(naFREE_FMat[,bf]~naFREE_FMat[,"outcome"]), data = naFREE_FMat)
      coef = coefficients(t_lm)[2,1]
      pval = coefficients(t_lm)[2,4]
      lmResults =  rbind(lmResults,c(bf,coef,pval))
      
    } 
  }
  
  lmResults = data.frame(lmResults,stringsAsFactors = F)
  colnames(lmResults) = c("Hallmark","Coefficient","P_val")
  lmResults$Coefficient = as.numeric(lmResults$Coefficient)
  lmResults$P_val = as.numeric(lmResults$P_val)
  adj_P_val = p.adjust(lmResults$P_val,method = "fdr")
  lmResults = cbind(lmResults,adj_P_val)
  
  return(lmResults)
  
}

#Filter general linear model
filterGLM = function(GLMres,THRS){
  
  SIG = GLMres[GLMres$adj_P_val<=THRS,]
  UP = SIG$Hallmark[SIG$Coefficient>0]
  DOWN = SIG$Hallmark[SIG$Coefficient<0]
  
  return(list("up"=UP,"down"=DOWN))
}

###############################################################3
# Method that generates exhaustively feature signatures 
###############################################################

getSignatures = function(cMAT,targetMAT,Features=NULL,time,g.pval,s.pval = NULL,model_type = "limma"){
  # For limma and GLM: It creates block of patients by Features + their combinations
  # Inputs:
  # - Count_matrix (cMAT) = count matrix containing the samples you wish to subset and model
  # - Annotation_file (targetMAT) = annotation of the count matrix
  # - Features = dictionary that contains the name of the features and the values to be used in the combinations
  # - Time = times to be used in the analysis
  # - globalTreshold = pval to be used
  #
  # Outputs:
  # - Signatures and direction of the change
  # - Results from the linear model
  # - Data used for each model and its annotation
  #
  
  if(is.null(Features)){ #If features is NULL, it creates time-specific models for the selected data
    # Objects that store the information of each test (Data,annotation,results,signatures)
    signatureDictionary = list()
    resultsDictionary = list()
    dataDictionary = list()
    annotDictionary = list()
    tag = "g"
    
    for(tm in time){ #Subselection divided by time
      tTargets = FilterTargets(targetTable = targetMAT,query_class= c("time"), query_values = c(tm),OR = FALSE) #Function that makes the filtering
      Ncases = sum(tTargets$outcome==1)
      Ncontrols = sum(tTargets$outcome==0)
      
      if(Ncases>2 & Ncontrols>2){
        tdata = cMAT[,as.character(tTargets$sample_mask_id)]
        dataDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tdata
        annotDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tTargets
      }
    }
    
    
    if(model_type == "limma"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
        signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = g.pval,FC_col = "logFC")
      }
      
    } else if(model_type == "glm"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
        signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = g.pval)
      }
      
    }
    
    return(list("DATA"=dataDictionary,"ANNOTATIONS"=annotDictionary,"RESULTS"=resultsDictionary,"SIGNATURES"=signatureDictionary))
    
  } else{ #If youy give a list of features for combinatorial analysis
    
    #1st, make a selection of all data and annotation
    
    signatureDictionary = list()
    resultsDictionary = list()
    dataDictionary = list()
    annotDictionary = list()
    tag = "g"
    
    for(tm in time){ #Subselection divided by time
      tTargets = FilterTargets(targetTable = targetMAT,query_class= c("time"), query_values = c(tm),OR = FALSE) #Function that makes the filtering
      Ncases = sum(tTargets$outcome==1)
      Ncontrols = sum(tTargets$outcome==0)
      
      if(Ncases>2 & Ncontrols>2){
        tdata = cMAT[,as.character(tTargets$sample_mask_id)]
        dataDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tdata
        annotDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tTargets
      }
    }
    
    #2nd, generate all the combinations
    
    Feat_names = names(Features)
    for(i in 1:length(Feat_names)){
      #Generate combinations of feauture classes, containing 1 to n number of features
      feature_class = combinations((length(Feat_names)),i,Feat_names)
      
      for(ir in 1:nrow(feature_class)){ #For each combination of features, obtain a list and get all the combinations of their values
        
        useful_vars = Features[feature_class[ir,]]
        useful_vars = c(useful_vars,list("time"=time))
        combination_df = expand.grid(useful_vars,stringsAsFactors = F)
        filtFeatures = colnames(combination_df)
        
        for(iq in 1:nrow(combination_df)){ #for every combination, filter patients
          
          query_values = as.character(combination_df[iq,])
          tTargets = FilterTargets(targetTable = targetMAT,query_class= filtFeatures, query_values = query_values,OR = FALSE) #Function that makes the filtering
          
          Ncases = sum(tTargets$outcome==1)
          Ncontrols = sum(tTargets$outcome==0)
          
          if(Ncases>2 & Ncontrols>2){
            
            tag = paste(paste(filtFeatures,query_values,sep="_"),collapse = ".")
            tdata = cMAT[,as.character(tTargets$sample_mask_id)]
            dataDictionary[[tag]] = tdata
            annotDictionary[[tag]] = tTargets
          }
        }
      }
    }
    
    # Apply_methods
    
    if(model_type == "limma"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        
        if(identical(grep("^g",dselection),integer(0))){
          resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = s.pval,FC_col = "logFC")
        } else{
          resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = g.pval,FC_col = "logFC")
        }
      }
      
    } else if(model_type == "glm"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        if(identical(grep("^g",dselection),integer(0))){
          resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = s.pval)
        } else{
          resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = g.pval)
        }
        
      }
      
    }
    
    return(list("DATA"=dataDictionary,"ANNOTATIONS"=annotDictionary,"RESULTS"=resultsDictionary,"SIGNATURES"=signatureDictionary))
    
  }
  
}

####################
# You can unlist the signatures to obtain a new set of features 
####################

unlistSignatures = function(Signature_Dictionary){
  
  All_signatures_dictionary_unlist = list()
  for(sign in names(Signature_Dictionary)){
    
    upT = paste(sign,"up",sep="_")
    downT = paste(sign,"down",sep="_")
    Flevel = Signature_Dictionary[[sign]]
    
    if(!is.null(Flevel)){
      if(!identical(Flevel[["up"]],character(0))){
        All_signatures_dictionary_unlist[[upT]] = Flevel[["up"]]
      }
      if(!identical(Flevel[["down"]],character(0))){
        All_signatures_dictionary_unlist[[downT]] = Flevel[["down"]]
      }
    }
  }
  
  return(All_signatures_dictionary_unlist)
  
}

############################
#
# Time Course Linear Models
#

runTimeCourse_test = function(countMAT, annotationDF, compCOL, Age0 = F){
  # It only supports 1-to-1 comparisons
  # A - B is always compared
  
  geneList = list()
  
  Gs = unique(annotationDF[,compCOL])
  ClassA = Gs[1]
  ClassB = Gs[2]
  
  Sclass = paste(ifelse(annotationDF[,compCOL]==ClassA,"ClassA","ClassB"),annotationDF$time,sep=".")
  annotationDF = cbind(annotationDF,Sclass)
  
  f = factor(annotationDF$Sclass,levels = sort(unique(annotationDF$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  countMAT = countMAT[,as.character(annotationDF$sample_mask_id)]
  
  fit <- lmFit(countMAT, design)
  
  if(Age0){
    cont.IA = makeContrasts(
      Dif03 = (ClassA.3 - ClassA.0) - (ClassB.3 - ClassB.0),
      Dif36 = (ClassA.6 - ClassA.3) - (ClassB.6 - ClassB.3),
      Dif69 = (ClassA.9 - ClassA.6) - (ClassB.9 - ClassB.6),
      #Dif129 = (ClassA.12 - ClassA.9) - (ClassB.12 - ClassB.9),
      levels = design
    )
    
    fit2 = contrasts.fit(fit,cont.IA)
    fit2 = eBayes(fit2)
    
    lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
    
    geneList[["global"]] = lmres
    
    for(i in 1:3){
      pr = paste("Period", as.character(i * -1))
      lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
      geneList[[pr]] = lmres
    }
    
    geneList[["ClassA"]] = ClassA
    geneList[["ClassB"]] = ClassB
    
  }else{
    cont.IA = makeContrasts(
      Dif03 = (ClassA.3 - ClassA.0) - (ClassB.3 - ClassB.0),
      Dif36 = (ClassA.6 - ClassA.3) - (ClassB.6 - ClassB.3),
      Dif69 = (ClassA.9 - ClassA.6) - (ClassB.9 - ClassB.6),
      Dif129 = (ClassA.12 - ClassA.9) - (ClassB.12 - ClassB.9),
      levels = design
    )
    
    fit2 = contrasts.fit(fit,cont.IA)
    fit2 = eBayes(fit2)
    
    lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
    
    geneList[["global"]] = lmres
    
    for(i in 1:4){
      pr = paste("Period", as.character(i * -1))
      lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
      geneList[[pr]] = lmres
    }
    
    geneList[["ClassA"]] = ClassA
    geneList[["ClassB"]] = ClassB
  }
  
  return(geneList)
}

#
#
#

runTimeCourse_test_multifactorial = function(countMAT_A, annotationDF_A, countMAT_B, annotationDF_B, compCOL, ClassA, ClassB,Age0 = F){
  # runTimeCourse_test but compares two different datasets, for example, age group
  # It calculate the difference of change between cases and controls, between a second comparison of groups
  
  geneList = list()
  
  GA_Sclass = paste("GA", ifelse(annotationDF_A[,compCOL]==ClassA,"ClassA","ClassB"),annotationDF_A$time,sep=".")
  GB_Sclass = paste("GB", ifelse(annotationDF_B[,compCOL]==ClassA,"ClassA","ClassB"),annotationDF_B$time,sep=".")
  
  annotationDF_A = cbind(annotationDF_A, Sclass = GA_Sclass)
  annotationDF_B = cbind(annotationDF_B, Sclass = GB_Sclass)
  
  annotationDF = rbind(annotationDF_A,annotationDF_B)
  
  f = factor(annotationDF$Sclass,levels = sort(unique(annotationDF$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  countMAT = cbind(countMAT_A,countMAT_B)
  countMAT = countMAT[,as.character(annotationDF$sample_mask_id)]
  
  fit <- lmFit(countMAT, design)
  
  if(Age0){
    cont.IA = makeContrasts(
      Dif03 = ((GA.ClassA.3 - GA.ClassA.0) - (GA.ClassB.3 - GA.ClassB.0)) - ((GB.ClassA.3 - GB.ClassA.0) - (GB.ClassB.3 - GB.ClassB.0)),
      Dif36 = ((GA.ClassA.6 - GA.ClassA.3) - (GA.ClassB.6 - GA.ClassB.3)) - ((GB.ClassA.6 - GB.ClassA.3) - (GB.ClassB.6 - GB.ClassB.3)),
      Dif69 = ((GA.ClassA.9 - GA.ClassA.6) - (GA.ClassB.9 - GA.ClassB.6)) - ((GB.ClassA.9 - GB.ClassA.6) - (GB.ClassB.9 - GB.ClassB.6)),
      #Dif129 = ((GA.ClassA.12 - GA.ClassA.9) - (GA.ClassB.12 - GA.ClassB.9)) - ((GB.ClassA.12 - GB.ClassA.9) - (GB.ClassB.12 - GB.ClassB.9)),
      levels = design
    )
    
    fit2 = contrasts.fit(fit,cont.IA)
    fit2 = eBayes(fit2)
    
    lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
    
    geneList[["global"]] = lmres
    
    for(i in 1:3){
      pr = paste("Period", as.character(i * -1))
      lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
      geneList[[pr]] = lmres
    }
    
    geneList[["ClassA"]] = ClassA
    geneList[["ClassB"]] = ClassB
    
  }else{
    cont.IA = makeContrasts(
      Dif03 = ((GA.ClassA.3 - GA.ClassA.0) - (GA.ClassB.3 - GA.ClassB.0)) - ((GB.ClassA.3 - GB.ClassA.0) - (GB.ClassB.3 - GB.ClassB.0)),
      Dif36 = ((GA.ClassA.6 - GA.ClassA.3) - (GA.ClassB.6 - GA.ClassB.3)) - ((GB.ClassA.6 - GB.ClassA.3) - (GB.ClassB.6 - GB.ClassB.3)),
      Dif69 = ((GA.ClassA.9 - GA.ClassA.6) - (GA.ClassB.9 - GA.ClassB.6)) - ((GB.ClassA.9 - GB.ClassA.6) - (GB.ClassB.9 - GB.ClassB.6)),
      Dif129 = ((GA.ClassA.12 - GA.ClassA.9) - (GA.ClassB.12 - GA.ClassB.9)) - ((GB.ClassA.12 - GB.ClassA.9) - (GB.ClassB.12 - GB.ClassB.9)),
      levels = design
    )
    
    fit2 = contrasts.fit(fit,cont.IA)
    fit2 = eBayes(fit2)
    
    lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
    
    geneList[["global"]] = lmres
    
    for(i in 1:4){
      pr = paste("Period", as.character(i * -1))
      lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
      geneList[[pr]] = lmres
    }
    
    geneList[["ClassA"]] = ClassA
    geneList[["ClassB"]] = ClassB
  }
  
  return(geneList)
}
