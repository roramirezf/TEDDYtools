############################
# TEDDY tools: Processing and Filtering
# Author: Ricardo Ramirez
# 
# This collection of functions have the objective of processing and filtering
# the raw counts based on the annotations provided by TEDDY
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################


#####################################################################################################
# TARGET LOADING
# The difference between these two versions is that the second one can be used in a complete data set
#

loadTargets = function(SELdisease,time_point="All",complete_case=NULL,targets){
  #
  # Function that provides a target table of the cases you are interested
  # call:
  # loadTargets(SELdisease = c("IA","T1D"), time_point = c("All", vector with specific times -integers-), complete_case = c("strict","flexible",NULL))
  # 
  # complete_case: It allows you to get only the samples with a corresponding control in a given time 
  # if NULL: this feature is deactivated
  # if flexible: it keeps all the patients that have at leats 1 case-control for the defined time points 
  # if strict: it keeps the patients that have 1 case-control in all specified time points
  #
  # Example:
  # test = loadTargets(SELdisease = "IA",time_point = "All",complete_case = "strict")
  #
  
  #Select disease
  selTargets = filter(targets,disease==SELdisease)
  
  #Select time point
  if(time_point=="All"){
    time_point= c(12,9,6,3,0)
    selTargets =  filter(selTargets,time %in% time_point)
  }else{
    selTargets =  filter(selTargets,time %in% time_point)
  }
  
  #Select data with complete information
  if(is.null(complete_case)){
    return(selTargets)
  }else{
    smids = c()
    for(case_control in unique(selTargets$CASE_IND)){
      strictFlag = TRUE
      cases = filter(selTargets,CASE_IND==case_control,outcome==1)
      controls = filter(selTargets,CASE_IND==case_control,outcome==0)
      ccsmids = c()
      
      for(t in time_point){
        
        caseFlag = t %in% cases$time
        controlFlag = t %in% controls$time
        time_flag = caseFlag & controlFlag
        
        if(time_flag){
          casetmp = filter(cases,time==t)
          controltmp = filter(controls,time==t)
          ccsmids = c(ccsmids,casetmp$sample_mask_id,controltmp$sample_mask_id)
        }else{
          if(complete_case == "strict"){
            strictFlag = FALSE
            break
          }
        }
      }
      
      if(strictFlag){
        smids = c(smids,ccsmids)
      }
      
    }
    
    selTargets =  filter(selTargets,sample_mask_id %in% smids)
    selTargets = arrange(selTargets,CASE_IND,time,outcome)
    if(sum(duplicated(selTargets$sample_mask_id))>0){
      bids = selTargets$sample_mask_id[duplicated(selTargets$sample_mask_id)]
      bcind = unique((filter(selTargets,sample_mask_id %in% bids))$CASE_IND)
      selTargets = filter(selTargets,!(CASE_IND %in% bcind))
    }
    
    return(selTargets)
  }
}

loadTargets_v2 = function(SELdisease,time_point="All",complete_case=NULL,targets){
  #
  # Function that provides a target table of the cases you are interested
  # call:
  # loadTargets(SELdisease = c("IA","T1D"), time_point = c("All", vector with specific times -integers-), complete_case = c("strict","flexible",NULL))
  # 
  # complete_case: It allows you to get only the samples with a corresponding control in a given time 
  # if NULL: this feature is deactivated
  # if flexible: it keeps all the patients that have at leats 1 case-control for the defined time points 
  # if strict: it keeps the patients that have 1 case-control in all specified time points
  #
  # Example:
  # test = loadTargets(SELdisease = "IA",time_point = "All",complete_case = "strict")
  #
  
  #Select disease
  selTargets = filter(targets,disease==SELdisease)
  
  #Select time point
  if(time_point[1]=="All"){
    time_point= unique(targets$time)
    selTargets =  filter(selTargets,time %in% time_point)
  }else{
    selTargets =  filter(selTargets,time %in% time_point)
  }
  
  #Select data with complete information
  if(is.null(complete_case)){
    return(selTargets)
  }else{
    smids = c()
    for(case_control in unique(selTargets$CASE_IND)){
      strictFlag = TRUE
      cases = filter(selTargets,CASE_IND==case_control,outcome==1)
      controls = filter(selTargets,CASE_IND==case_control,outcome==0)
      ccsmids = c()
      
      for(t in time_point){
        
        caseFlag = t %in% cases$time
        controlFlag = t %in% controls$time
        time_flag = caseFlag & controlFlag
        
        if(time_flag){
          casetmp = filter(cases,time==t)
          controltmp = filter(controls,time==t)
          ccsmids = c(ccsmids,casetmp$sample_mask_id,controltmp$sample_mask_id)
        }else{
          if(complete_case == "strict"){
            strictFlag = FALSE
            break
          }
        }
      }
      
      if(strictFlag){
        smids = c(smids,ccsmids)
      }
      
    }
    
    selTargets =  filter(selTargets,sample_mask_id %in% smids)
    selTargets = arrange(selTargets,CASE_IND,time,outcome)
    if(sum(duplicated(selTargets$sample_mask_id))>0){
      bids = selTargets$sample_mask_id[duplicated(selTargets$sample_mask_id)]
      bcind = unique((filter(selTargets,sample_mask_id %in% bids))$CASE_IND)
      selTargets = filter(selTargets,!(CASE_IND %in% bcind))
    }
    
    return(selTargets)
  }
}

#####################################################################################################
# TARGET Modification
# Functions used to add more information into target tables
#

addFeatures = function(newDFdata,TARGETS,selCol){
  #
  # Generates a new target table with a new column of data
  #
  # inputs:
  # newDFdata = New -data frame- to merge
  # TARGETS = Original target table (labels as in Teddy website)
  # selCol = name of the column used for merging
  #
  # output:
  # New target file
  #
  newTARGETS = left_join(TARGETS,newDFdata,by =selCol)
  return(newTARGETS)
}

#################
#
# Filter the targets by a MetaVariable, and the classes you like
#

FilterTargets = function(targetTable,query_class,query_values,OR = FALSE){
  #
  # A filter function that gives you the rows from targetTable that have all the features defined in query_class/query_values (OR=FALSE)
  # or at least one feature (OR=TRUE)
  #
  
  # Make individuals
  if(length(query_class) == length(query_values)){
    EvaluationMat = c()
    for(i in 1:length(query_class)){
      EvalVector = as.character(targetTable[,query_class[i]]) == query_values[i]
      EvaluationMat = cbind(EvaluationMat,EvalVector)
    }
  }else{
    print("Length of query class and query values are not the same")
    return()
  }
  
  if(OR==FALSE){
    
    EvalVector = rowSums(EvaluationMat)
    EvalVector = EvalVector == length(query_values)
    return(targetTable[EvalVector,])
    
  }else{
    EvalVector = rowSums(EvaluationMat)
    EvalVector = EvalVector > 0
    return(targetTable[EvalVector,])
  }
}

#####################################################################################################
# NORMALIZATION
# X within and group normalization
#

#X-within normalization
#It deals with 1-3, 1-1 design, and with incomplete Data (for example, Metabolomics)
normGEX_v2 = function(expMAT,targetTable,time = c(12,9,6,3,0)){
  #
  # Normalises gene expression by substracting for each time specific case-control case, the mean signal
  # X WITHIN
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(case_control in unique(targetTable$CASE_IND)){
    subsetTargets = filter(targetTable,CASE_IND==case_control)
    for(t in time){
      if(t %in% subsetTargets$time){
        
        case_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==1))$sample_mask_id)
        control_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==0))$sample_mask_id)
        
        # Verify how many controls you have, just to decide which version of the normalization is going to be used
        
        if(length(control_id)>1){
          # Get the means of the controls, and then the mean of means
          control_mean = rowMeans(expMAT[,control_id],na.rm = TRUE)
          mean_mat = cbind(expMAT[,case_id],control_mean)
          
          #sids_t =  as.character((filter(targetTable,CASE_IND==case_control,time==t))$sample_mask_id)
          subsetMat = expMAT[,c(case_id,control_id)]
          subsetMat = subsetMat - rowMeans(mean_mat,na.rm = TRUE)
          expMAT[,c(case_id,control_id)] = subsetMat
        } else if(length(control_id)==1){
          
          sids_t =  c(case_id,control_id)
          subsetMat = expMAT[,sids_t]
          subsetMat = subsetMat - rowMeans(subsetMat,na.rm = TRUE)
          expMAT[,sids_t] = subsetMat
          
        }
        
      }
    }
  }
  return(expMAT)
}
# First version
normGEX = function(expMAT,targetTable){
  #
  # Normalises gene expression by substracting for each time specific case-control case, the mean signal
  # Note: Only works to normalise the original 5 time points (12,9,6,3) 
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(case_control in unique(targetTable$CASE_IND)){
    subsetTargets = filter(targetTable,CASE_IND==case_control)
    for(t in c(12,9,6,3,0)){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(targetTable,CASE_IND==case_control,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}

#Group Normalization
#This function normalizes by group... In our case... Outcome, nonetheless is not used in any analysis
gnormGEX_v2 = function(expMAT,targetTable,timepoints){
  #
  # Normalizes gene expression by substracting for each outcome group, the mean signal
  # call:
  # normGEX(expMAT=expressionMatrix,targetTable=target table generated by loadTargets or a table with CASE_IND,time columns)
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(oc in unique(targetTable$outcome)){
    subsetTargets = filter(targetTable,outcome==oc)
    for(t in timepoints){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(subsetTargets,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}
#First version
gnormGEX = function(expMAT,targetTable){
  #
  # Normalizes gene expression by substracting for each outcome group, the mean signal
  # call:
  # normGEX(expMAT=expressionMatrix,targetTable=target table generated by loadTargets or a table with CASE_IND,time columns)
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(oc in unique(targetTable$outcome)){
    subsetTargets = filter(targetTable,outcome==oc)
    for(t in c(12,9,6,3,0)){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(subsetTargets,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}

#####################################################################################################
# Annotation of GEX matrices
# Different versions of annotations
#

annotateGEX = function(expMAT,annotation_set){
  #
  # Annotates a gene expression matrix and takes the mean of repeated symbols
  # inputs:
  # -expMAT = expression matrix
  # -annotation_set = your microarray database (AnnotationDBi)
  #
  # output:
  # -Annotated expression matrix
  #
  probe_names = rownames(expMAT)
  probe_symbol = na.omit(data.frame(mapIds(annotation_set, keys=probe_names, keytype = "PROBEID", column="SYMBOL",multiVals = "first"),stringsAsFactors = F))
  #FIlter exprmat with known probes
  expMAT = expMAT[rownames(probe_symbol),]
  rownames(expMAT) = probe_symbol[,1]
  expMAT = avereps(expMAT,ID=rownames(expMAT))
  return(expMAT)
}

annotateGEX_v2 = function(expMAT, annotation_set, annotation_column){
  #
  # Annotates a gene expression matrix and takes the mean of repeated symbols
  # inputs:
  # -expMAT = expression matrix
  # -annotation_set = your microarray database (AnnotationDBi)
  # -annotation_column = ENTREZID,SYMBOL
  # output:
  # -Annotated expression matrix
  #
  probe_names = rownames(expMAT)
  probe_symbol = na.omit(data.frame(mapIds(annotation_set, keys=probe_names, keytype = "PROBEID", column = annotation_column,multiVals = "first"),stringsAsFactors = F))
  #FIlter exprmat with known probes
  expMAT = expMAT[rownames(probe_symbol),]
  rownames(expMAT) = probe_symbol[,1]
  expMAT = avereps(expMAT,ID=rownames(expMAT))
  return(expMAT)
}

#####################################################################################################
# Recovery of Gene Expression Matrices
# Different versions
#

getExprmat = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Annotation of matrices
  expmat = annotateGEX(expMAT = expmat,annotation_set = illuminaHumanv4.db)
  
  return(expmat)
}

#Exactly as above, but allows you to use ENTREZ instead of SYMBOLS
getExprmat_v2 = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE, annotation_column = "SYMBOL"){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Annotation of matrices
  expmat = annotateGEX_v2(expMAT = expmat,annotation_set = illuminaHumanv4.db, annotation_column = annotation_column)
  
  return(expmat)
}

#The most robust version of getting an expression matrix
# time = the unique time points present in STARGETS *** YOU MUST REMEMBER THAT
getExprmat_v3 = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE, time,annotation_column = "SYMBOL"){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX_v2(expMAT=expmat,targetTable=STARGETS,time = time) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX_v2(expMAT=expmat,targetTable=STARGETS,timepoints = time) 
  }
  # Annotation of matrices
  expmat = annotateGEX_v2(expMAT = expmat,annotation_set = illuminaHumanv4.db, annotation_column = annotation_column)
  
  return(expmat)
}

#
# Functions that write Lists of dataframes/vectors
#

write_list = function(Olist, ResultFolder,ExtraDef=""){
  
  gFolder = paste("mkdir ",ResultFolder,sep = "")
  system(gFolder)
  
  for(E in names(Olist)){
    
    write.table(Olist[[E]],file = paste(ResultFolder,E,"_",ExtraDef,".txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
    
  }
  
}

write_genelists = function(Olist, ResultFolder,ExtraDef=""){
  
  gFolder = paste("mkdir ",ResultFolder,sep = "")
  system(gFolder)
  
  for(E in names(Olist)){
    
    glist = Olist[[E]]
    
    glist_Annotation = AnnotationDbi::select(illuminaHumanv4.db, keys=glist,columns = c("SYMBOL","GENENAME","ENTREZID"),keytype="SYMBOL")
    
    glist_Annotation = glist_Annotation[!duplicated(glist_Annotation$SYMBOL),]
    
    glist_Annotation = arrange(glist_Annotation,SYMBOL)
    
    write.table(glist_Annotation,file = paste(ResultFolder,E,"_",ExtraDef,".txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
    
  }
  
}


