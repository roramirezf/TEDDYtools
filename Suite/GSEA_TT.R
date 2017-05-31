############################
# TEDDY tools: GSEA
# Author: Ricardo Ramirez and Broad
# 
# This collection of functions have the objective of generating function matrices using GSEA
#  -
#
# Some functions may have more than 1 version because they were adapted for general pospuses
############################

getFunctionMatrix_v2 = function(eMat, eTargets, Signature_Dict, weighted=TRUE){
  Fmat = c()
  Ledge = list()
  for(geneSet_names in names(Signature_Dict)){
    print(geneSet_names)
    geneSet = Signature_Dict[[geneSet_names]]
    ES_results = GSEA_ES_v3(eMat = eMat,gene.set = geneSet,eTargets = eTargets,weighted = weighted)
    ES_vector = ES_results$ES_vector
    Fmat = rbind(Fmat, ES_vector)
    Ledge[[geneSet_names]] = ES_results$LeadingEdge
  }
  rownames(Fmat) = names(Signature_Dict)
  return(list(ESmat = Fmat, LeadingEdge = Ledge))
}

GSEA_ES_v3 = function(eMat, eTargets, gene.set, weighted = TRUE){
  # For a given set of samples (expression matrix) calculates the enrichment score of a gene.set
  # If you decide to weight the list, Signal2Noise is calculated. (Broad implementation: difference of means scaled by the standard deviation)
  #
  ES_list = c()
  LeadingEdge_list = list()
  
  if(weighted == FALSE){
    
    for(i in 1:ncol(eMat)){
      sampleName = colnames(eMat)[i]
      gene.list = names(sort(eMat[,i],decreasing = T))
      #ADD NULL robustness
      int.gene.set = intersect(gene.list,gene.set)
      if(length(int.gene.set>0)){
        ES_results = GSEA.EnrichmentScore(gene.list = gene.list,gene.set = int.gene.set,weighted.score.type = 0)
        ES = ES_results$ES
        ES_list = c(ES_list,ES)
        peak = ES_results$arg.ES
        if(ES>0){
          coordinates = as.logical(ES_results$indicator[1:peak])
          leadingedge = gene.list[1:peak]
          leadingedge =  leadingedge[coordinates]
          LeadingEdge_list[[sampleName]] = leadingedge
        }else{
          coordinates = as.logical(ES_results$indicator[peak:length(gene.list)])
          leadingedge = gene.list[peak:length(gene.list)]
          leadingedge =  leadingedge[coordinates]
          LeadingEdge_list[[sampleName]] = leadingedge
        }
      } else{
        ES_list = c(ES_list,NA)
        LeadingEdge_list[[sampleName]] = NA
      }
      
    }
    names(ES_list) = colnames(eMat)
    return(list(ES_vector = ES_list, LeadingEdge = LeadingEdge_list))
    
    
  }
  
  if(weighted == TRUE){
    
    O = GSEA.GeneRanking(A=eMat, class.labels=eTargets$outcome, gene.labels=rownames(eMat), nperm=2, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F)
    obs.correl.matrix = O$obs.s2n.matrix
    obs.s2n = apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
    obs.s2n = sort(obs.s2n, decreasing=T)   
    
    for(i in 1:ncol(eMat)){
      sampleName = colnames(eMat)[i]
      gene.list = names(sort(eMat[,i],decreasing = T))
      signal.noise.scores = obs.s2n[gene.list]
      int.gene.set = intersect(gene.list,gene.set)
      ES_results = GSEA.EnrichmentScore(gene.list = gene.list,gene.set = int.gene.set,weighted.score.type = 1,correl.vector = signal.noise.scores)
      ES = ES_results$ES
      ES_list = c(ES_list,ES)
      peak = ES_results$arg.ES
      if(ES>0){
        coordinates = as.logical(ES_results$indicator[1:peak])
        leadingedge = gene.list[1:peak]
        leadingedge =  leadingedge[coordinates]
        LeadingEdge_list[[sampleName]] = leadingedge
      }else{
        coordinates = as.logical(ES_results$indicator[peak:length(gene.list)])
        leadingedge = gene.list[peak:length(gene.list)]
        leadingedge =  leadingedge[coordinates]
        LeadingEdge_list[[sampleName]] = leadingedge
      }
    }
    names(ES_list) = colnames(eMat)
    return(list(ES_vector = ES_list, LeadingEdge = LeadingEdge_list))
    
  }
  
}

#########################################################################

#Broad Institute Functions - Adapted

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag  <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F) { 
  
  # This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap  
  # subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation 
  # in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced. 
  # It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
  # the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for 
  # all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second 
  # matrix for the null distribution will still have the values for the random permutations 
  # (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
  # It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain 
  # smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity 
  # checks before trusting the code.
  #
  # Inputs:
  #   A: Matrix of gene expression values (rows are genes, columns are samples) 
  #   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
  #   gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of the expression matrix 
  #   nperm: Number of random permutations/bootstraps to perform 
  #   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
  #   sigma.correction: Correction to the signal to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package) 
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F) 
  #   reverse.sign: Reverse direction of gene list (default = F)
  #
  # Outputs:
  #   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios (rows are genes, columns are permutations or bootstrap subsamplings
  #   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
  #   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
  #   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  A <- A + 0.00000001
  
  N <- length(A[,1])
  Ns <- length(A[1,])
  
  subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  
  order.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
  s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  
  obs.gene.labels <- vector(length = N, mode="character")
  obs.gene.descs <- vector(length = N, mode="character")
  obs.gene.symbols <- vector(length = N, mode="character")
  
  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)
  
  gc()
  
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
  
  for (r in 1:nperm) {
    class1.subset <- sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
    class2.subset <- sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
    class1.subset.size <- length(class1.subset)
    class2.subset.size <- length(class2.subset)
    subset.class1 <- rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] <- 1
      }
    }
    subset.class2 <- rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] <- 1
      }
    }
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size/Ns
    fraction.class2 <- class2.size/Ns
    
    if (permutation.type == 0) { # random (unbalanced) permutation
      full.subset <- c(class1.subset, class2.subset)
      label1.subset <- sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        m1 <- sum(!is.na(match(label1.subset, i)))
        m2 <- sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] <- m1
        reshuffled.class.labels2[i, r] <- m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
      
    } else if (permutation.type == 1) { # proportional (balanced) permutation
      
      class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
      class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 <- sum(!is.na(match(class1.label1.subset, i)))
          m2 <- sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          m1 <- sum(!is.na(match(class2.label1.subset, i)))
          m2 <- sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    }
  }
  
  # compute S2N for the random permutation matrix
  
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  s2n.matrix <- M1/S1
  
  if (reverse.sign == T) {
    s2n.matrix <- - s2n.matrix
  }
  gc()
  
  for (r in 1:nperm) {
    order.matrix[, r] <- order(s2n.matrix[, r], decreasing=T)            
  }
  
  # compute S2N for the "observed" permutation matrix
  
  P <- class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  } 
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  obs.s2n.matrix <- M1/S1
  gc()
  
  if (reverse.sign == T) {
    obs.s2n.matrix <- - obs.s2n.matrix
  }
  
  for (r in 1:nperm) {
    obs.order.matrix[,r] <- order(obs.s2n.matrix[,r], decreasing=T)            
  }
  
  return(list(s2n.matrix = s2n.matrix, 
              obs.s2n.matrix = obs.s2n.matrix, 
              order.matrix = order.matrix,
              obs.order.matrix = obs.order.matrix))
}

###############################################
#
#
# REGULAR ENRICHMENT METHODS
#
#
###############################################

#
#     A!        ( A )
# ---------- =  ( a )
#  a!(A-a)!
#
getCombinations = function(A, a){
  return( factorial(A)/(factorial(a) * factorial(A-a)) )
}

# Code the Hypergeometric test
# Over-representation
# The probability of randomly drawing k or more successes from the population in n total draws
#
#     (K)(N-K)
#     (k)(n-k)
#   ------------ =  probability of randomly drawing k successes from the population
#        (N)
#        (n)
#
#
# Where:
# N = number of annotated genes in a specific data base 
# K = number of genes in pathway Alpha
# And:
# n = number of genes in the genelist
# k = number of genes in the genelist that are in pathway Alpha
#
HypergeometricProb = function(N,K,n,k){
  
  numerator = choose(K,k) * choose(N-K,n-k)
  denominator = choose(N,n)
  
  return(numerator/denominator)
  
}

HypergeometricTest = function(overRepres = TRUE, N,K,n,k){
  
  if(overRepres){
    
    if(n>K){
      maxk = K
    }else{
      maxk = n
    }
    
    p=0
    for(i in k:maxk){
      p = p + HypergeometricProb(N,K,n,k=i)
    }
    
    return(p)
    
  }else{
    
    p=0
    
    for(i in 0:k){
      p = p + HypergeometricProb(N,K,n,k=i)
    }
    
    return(p)
  }
}


#
# GSE analysis using hypergeometric tests
#

GSE_analysis = function(geneList,Annotation_DB){
  
  # TO DO
  # Filter AnnotationDB
  # test = lapply(Annotation_DB,length)
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    #ResultsDF[gset,"p_value"] = HypergeometricTest(overRepres = TRUE,N = N,K = K,n = n,k = k)
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  return(ResultsDF)
  
}


#################################################################
#
# Piano Analysis
#

runMultiGSA = function(GSC, GLS){
  pval = data.frame(GLS$adj.P.Val)
  rownames(pval) = rownames(GLS)
  Fchange = data.frame(GLS$logFC)
  rownames(Fchange) = rownames(GLS)
  tval = data.frame(GLS$t)
  rownames(tval) = rownames(GLS)
  
  gsaRes1 = runGSA(tval,gsc=GSC,geneSetStat = "mean",gsSizeLim=c(5,300))
  gsaRes2 = runGSA(tval,gsc=GSC,geneSetStat = "median",gsSizeLim=c(5,300))
  gsaRes3 = runGSA(tval,gsc=GSC,geneSetStat = "sum",gsSizeLim=c(5,300))
  gsaRes4 = runGSA(tval,gsc=GSC,geneSetStat = "maxmean",gsSizeLim=c(5,300))
  
  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4)
  names(resList) <- c("mean","median","sum","maxmean")
  
  return(resList)
}


