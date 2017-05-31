############################
# TEDDY tools: Transrcriptional Signatures Intersect
# Author: Ricardo Ramirez
# 
# This collection of functions have the objective of comparing sets of transcriptional signatures
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Functions to plot intersection
similix = function(vectorA, vectorB, n = 14302, significancePercentage = .50){
  I = length(intersect(vectorA,vectorB))
  similix = I/n
  
  referenceV = min(c(length(vectorA)),length(vectorB))
  
  if(I >= round(referenceV*significancePercentage)){
    sig = 1
  }else{
    sig = 0
  }
  
  return(c(similix,sig))
}

jaccardIx = function(set1,set2){
  intersect_val = length(intersect(set1,set2))
  union_val = length(union(set1,set2))
  return(intersect_val/union_val)
}


getIntersect_matrix = function(SignatureDictionary, significancePercentage = .50){
  # Define matrices that are used for plotting (Directed, non-directed and their corresponding significance matrix)
  # NOTE: All NaN are converted to 0 and the similarity indexes are not rezised
  # Inputs: 
  # -SIgnatureDictionary: List with signatures to be compared
  # -significancePercentage: how much an intersection should cover from the smallest signature in a given comparison, to be considered significant
  # Outputs:
  # -list containing 4 matrices, 
  
  sigNames = names(SignatureDictionary)
  N = length(unique(unlist(SignatureDictionary)))
  
  sigmat = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(sigmat) = sigNames
  colnames(sigmat) = sigNames
  
  sigmat_sig = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(sigmat_sig) = sigNames
  colnames(sigmat_sig) = sigNames
  
  dir_sigmat = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(dir_sigmat) = sigNames
  colnames(dir_sigmat) = sigNames
  
  dir_sigmat_sig = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(dir_sigmat_sig) = sigNames
  colnames(dir_sigmat_sig) = sigNames
  
  # Fill the matrices
  
  for(i in 1:(length(sigNames)-1)){
    
    #Fill identities
    
    i_sig = SignatureDictionary[[sigNames[i]]]
    
    i_sig_up = i_sig$up
    i_sig_down = i_sig$down
    i_sig = unique(unlist(i_sig))
    
    simil_score = similix(i_sig, i_sig, n = N, significancePercentage = significancePercentage)
    
    sigmat[sigNames[i],sigNames[i]] = simil_score[1]
    sigmat_sig[sigNames[i],sigNames[i]] = simil_score[2]
    
    dir_sigmat[sigNames[i],sigNames[i]] = simil_score[1]
    dir_sigmat_sig[sigNames[i],sigNames[i]] = simil_score[2]
    
    for(j in (i+1):length(sigNames)){
      
      #Fill comparisons
      
      j_sig = SignatureDictionary[[sigNames[j]]]
      
      j_sig_up = j_sig$up
      j_sig_down = j_sig$down
      j_sig = unlist(j_sig)
      
      jix = similix(i_sig,j_sig, n = N, significancePercentage = significancePercentage)
      up_jix = similix(i_sig_up,j_sig_up, n = N, significancePercentage = significancePercentage)
      down_jix = similix(i_sig_down,j_sig_down, n = N, significancePercentage = significancePercentage)
      
      #Fill undirected
      sigmat[sigNames[i],sigNames[j]] = jix[1]
      sigmat[sigNames[j],sigNames[i]] = jix[1]
      sigmat_sig[sigNames[i],sigNames[j]] = jix[2]
      sigmat_sig[sigNames[j],sigNames[i]] = jix[2]
      
      #Fill directed
      dir_sigmat[sigNames[i],sigNames[j]] = up_jix[1]
      dir_sigmat_sig[sigNames[i],sigNames[j]] = up_jix[2]
      
      dir_sigmat[sigNames[j],sigNames[i]] = down_jix[1]
      dir_sigmat_sig[sigNames[j],sigNames[i]] = down_jix[2]
    }
  }
  
  i_sig = SignatureDictionary[[sigNames[j]]]
  
  i_sig_up = i_sig$up
  i_sig_down = i_sig$down
  i_sig = unique(unlist(i_sig))
  
  simil_score = similix(i_sig, i_sig, n = N, significancePercentage = significancePercentage)
  
  sigmat[sigNames[j],sigNames[j]] = simil_score[1]
  sigmat_sig[sigNames[j],sigNames[j]] = simil_score[2]
  
  dir_sigmat[sigNames[j],sigNames[j]] = simil_score[1]
  dir_sigmat_sig[sigNames[j],sigNames[j]] = simil_score[2]
  
  # Delete all NAN
  
  sigmat[is.nan(sigmat)] = 0
  sigmat_sig[is.nan(sigmat_sig)] = 0
  dir_sigmat[is.nan(dir_sigmat)] = 0
  dir_sigmat_sig[is.nan(dir_sigmat_sig)] = 0
  
  return(list("undirected" = sigmat,"undirected_significance" = sigmat_sig,"directed" = dir_sigmat,"directed_significance" = dir_sigmat_sig))
}


plotIntersection = function(Signature_mats){
  # Normalises similarity matrix and creates plots of the matrices, directed and not directed 
  # INPUT: 
  # - Signature_mats: getIntersect_matrix output
  #
  directed = Signature_mats$directed
  directed_sig = Signature_mats$directed_significance
  
  undirected = Signature_mats$undirected
  undirected_sig = Signature_mats$undirected_significance
  
  # First Plot directed mats
  
  up_regulated = get_upper_tri(directed)
  regulation="upregulated"
  up_df = cbind(melt(up_regulated,na.rm=TRUE),regulation)
  
  up_regulated_sig = get_upper_tri(directed_sig)
  up_df_sig = melt(up_regulated_sig,na.rm=TRUE)
  
  down_regulated = get_lower_tri(directed)
  regulation="downregulated"
  down_df = cbind(melt(down_regulated,na.rm=TRUE),regulation)
  
  down_regulated_sig = get_lower_tri(directed_sig)
  down_df_sig = melt(down_regulated_sig,na.rm=TRUE)
  
  significance = c(up_df_sig$value,down_df_sig$value)
  mat_df = rbind(up_df,down_df)
  mat_df = cbind(mat_df, significance)
  
  #Manipulate data
  colnames(mat_df) = c("SigA","SigB","Similarity_Ix","Regulation","Significance")
  mat_df$Regulation = as.character(mat_df$Regulation)
  mat_df$Regulation[mat_df$SigA == mat_df$SigB] = "identity"
  mat_df = unique(mat_df)
  
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==0]=NA
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==Inf]=NA
  mat_df$Similarity_Ix[is.nan(mat_df$Similarity_Ix)]=NA
  mat_df$Similarity_Ix = mat_df$Similarity_Ix/max(mat_df$Similarity_Ix,na.rm = T)
  
  mat_df$Regulation = as.factor(mat_df$Regulation)
  mat_df$Significance = factor(mat_df$Significance,levels=c(1,0))
  
  p = ggplot(data = mat_df, aes(x=SigA, y=SigB, fill=Regulation)) + geom_tile(colour="black") + geom_point(aes(size = Similarity_Ix,colour=Significance)) + scale_colour_manual(values=c("gold","black")) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  
  # Plot undirected mats
  
  up_regulated = get_upper_tri(undirected)
  up_df = melt(up_regulated,na.rm=TRUE)
  
  up_regulated_sig = get_upper_tri(undirected_sig)
  up_df_sig = melt(up_regulated_sig,na.rm=TRUE)
  significance = up_df_sig$value
  
  mat_df = cbind(up_df, significance)
  
  #Manipulate data
  colnames(mat_df) = c("SigA","SigB","Similarity_Ix","Significance")
  
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==0]=NA
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==Inf]=NA
  mat_df$Similarity_Ix[is.nan(mat_df$Similarity_Ix)]=NA
  mat_df$Similarity_Ix = mat_df$Similarity_Ix/max(mat_df$Similarity_Ix,na.rm = T)
  
  mat_df$Significance = factor(mat_df$Significance,levels=c(1,0))
  
  p = ggplot(data = mat_df, aes(x=SigA, y=SigB)) + geom_tile(colour="black",fill="lightgrey") + geom_point(aes(size = Similarity_Ix,colour=Significance)) + scale_colour_manual(values=c("red", "black")) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  
}

#
# Plot similarity of selected tags
# 

SignatureIntersect = function(SignatureDictionary, collapse = FALSE, tags = c("Time"), significancePercentage = .50){
  # For a given subset of signatures, create directed and not directed similarity matrices and plot them.
  # You can create tags to create subgroups of signatures and collapse them, or simply just to reduce the set of signatures in a signature dictionary
  #
  # Inputs: 
  # -SignatureDictionary: A list of signatures
  # -collapse: TRUE if you like to create and evaluate TAG signatures
  # -tags: strings or regular expressions to subset the signatures. It is relative to names(SignatureDictionary)
  # -significancePercentage: how much an intersection should cover from the smallest signature in a given comparison, to be considered significant
  #
  # Outputs:
  # -list with similarity matrices and their respective significance matrices
  # -plot of the matrices
  #
  #
  #
  selectedSignatures = c()
  
  if(length(tags)<2){
    collapse == FALSE
    print("Collapse has been set to FALSE: Only 1 tag was provided")
  }
  
  if(collapse == TRUE){
    collapseDict = list()
  }
  
  for(tag in tags){
    selectedNames = names(SignatureDictionary)[grep(tag,names(SignatureDictionary))]
    selectedSignatures = c(selectedSignatures, selectedNames)
    
    if(collapse == TRUE){
      
      collapseDict[[tag]] = list("up"=c(),"down"=c())
      
      up_genes = c()
      down_genes = c()
      
      for(n in selectedNames){
        up_genes = c(up_genes,SignatureDictionary[[n]]$up)
        down_genes = c(up_genes,SignatureDictionary[[n]]$down)
      }
      
      up_genes = unique(up_genes)
      down_genes = unique(down_genes)
      
      #Not allow to have genes in both signatures
      #up_genes = setdiff(up_genes, down_genes)
      #down_genes = setdiff(down_genes, up_genes)
      
      collapseDict[[tag]]$up = up_genes
      collapseDict[[tag]]$down = down_genes
      
    }
  }
  
  if(collapse == TRUE){
    SignatureDictionary = collapseDict
  }else{
    SignatureDictionary = SignatureDictionary[selectedSignatures]
  }
  
  #Signatures matrix
  
  Signature_mats = getIntersect_matrix(SignatureDictionary = SignatureDictionary, significancePercentage = significancePercentage)
  
  #Plot Signatures
  plotIntersection(Signature_mats)
  
  return(Signature_mats)
  
}