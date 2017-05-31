############################
# TEDDY tools: VisualManager
# Author: Ricardo Ramirez
# 
# This collection of plotting functions can be used in most of TEDDY's data
#
# Some functions may have more than 1 version because they were adapted for general porpuses
############################

#
#
# FEATURE DYNAMICS
#
#

# This function allows you to compare feature's dynamics of cases/controls (or the variable coded as "outcome")
# A linear model is performed at each timepoint for each feature to evaluate the differences of the classes defined in "outcome"

genedynamics = function(EXPRMAT, GENES, STARGETS, TRSH = 0.05, SIGNIFICANCE = FALSE,CASE_CONTROL=FALSE){
  # Input:
  # EXPRMAT - count matrix
  # GENES - features you want to test
  # STARGETS - target matrix
  # TRSH - Non-corrected p-value threshold for significance
  # SIGNIFICANCE - TRUE, if you only want to plot the genes with at least one point significant
  # CASE_CONTROL - TRUE, if your targets have a CASE_IND column and you are interested in case-control groups (to be put in the plot, not recommended)
  #
  #Load libraries
  library(gridExtra)
  library(ggplot2)
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    #For each gene at each time adjust a linear model and save results
    lm_results = c()
    for(GENE in GENES){
      for(t in unique(genedframe$time)){
        t_GEX = genedframe[genedframe$time==t,]
        t_lm = summary(lm(t_GEX[,GENE]~t_GEX[,"outcome"]), data = t_GEX)
        coef = coefficients(t_lm)[2,1]
        pval = coefficients(t_lm)[2,4]
        lm_results =  rbind(lm_results,c(t,GENE,coef,pval))
      }
    }
    #Manipulate pval results
    colnames(lm_results) = c("TIME","GENE","COEF","PVAL")
    lm_results = data.frame(lm_results,stringsAsFactors = F)
    lm_results[,2] = as.factor(lm_results[,2]) 
    lm_results[,3] = as.numeric(lm_results[,3])
    lm_results[,4] = p.adjust(p=as.numeric(lm_results[,4]),method = "fdr")
    
    for(GENE in GENES){
      print(GENE)
      significant = FALSE
      gene_results = lm_results[lm_results$GENE==GENE,]
      
      if(SIGNIFICANCE==TRUE){
        if(sum(gene_results$PVAL < TRSH) > 0){
          significant = TRUE
        }
      }else{
        significant = TRUE
      }
      
      gene_results$PVAL = -log(gene_results$PVAL)
      
      if(significant){
        
        vp = ggplot(data=gene_results, aes(x=COEF, y=PVAL, label=TIME),guide=FALSE) + geom_point(colour="white", fill="lightblue", shape=21) +
          geom_text(size=5,colour="black") + scale_size_continuous(range = c(8, 18)) + scale_x_continuous(name="Coefficient") +
          scale_y_continuous(name="-log(pval)") + geom_hline(yintercept=-log(TRSH), linetype="dashed") + theme(plot.margin = unit(c(1,3.3,1,1), "cm"),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        
        if(CASE_CONTROL){
          GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = CASE_IND, group = mask_id, linetype = outcome)) + 
            labs(y = GENE) + geom_point() + geom_line() + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ylab("Expression") + theme(plot.margin = unit(c(1,0.6,1,0.6), "cm"), axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        }else{
          GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = outcome, group = mask_id)) + 
            labs(y = GENE) + geom_point() + geom_line() + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ylab("Expression") + theme(plot.margin = unit(c(1,0.6,1,0.6), "cm"),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        }
        
        fgenedframe = genedframe
        fgenedframe$time = factor(genedframe$time,levels=c(12,9,6,3,0))
        gp = ggplot(fgenedframe, aes(x = time, y = fgenedframe[,GENE])) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + 
          stat_summary(fun.y=mean, geom="line", aes(group=outcome, col=outcome),size=1.5)  + ylab("Expression") +
          theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm"))
        
        
        grid.arrange(GEX_plot,gp,vp, nrow=3,top = GENE)
        
      }
    }
    
    return(lm_results)
  }
  print("We failed to find genes in expression matrix") 
  return(NULL)
}

# This function is similar to genedynamics, nonetheless, a significance test is not performed

plotgenetrajectory = function(EXPRMAT, GENES, STARGETS, CASE_CONTROL=FALSE){
  
  #Load libraries
  library(gridExtra)
  library(ggplot2)
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    
    for(GENE in GENES){
      
      if(CASE_CONTROL){
        GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = CASE_IND, group = mask_id, linetype = outcome)) + 
          labs(y = GENE) + geom_point() + geom_line() + theme(plot.margin = unit(c(1,1,1,1), "cm"))
        print(GEX_plot)
      }else{
        GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = outcome, group = mask_id)) + 
          labs(y = GENE) + geom_point() + geom_line() + theme(plot.margin = unit(c(1,1,1,1), "cm"))
        print(GEX_plot)
      }
      
    }
    
  }
  else{
    print("We failed to find genes in expression matrix") 
    return(NULL)
  }
}

# This function compares feature's dynamics of two classes defined in "outcome", 
# and performs a time specific linear model to assess the significance of the difference   
plotES = function(ES_list,sigTargets,TRSH=0.05,main_title){
  # This function compares the Enrichment Scores of patients and controls, with a T-test
  # Input:
  # -ES_list: A vector with the expression values of the feature you want to analyse
  # -sigTargets =  Annotation of the count matrix where ES_list was obtained
  # -TRSH =  Significance threshold of the T-test
  # -main_title = plot Title
  # Output:
  # -Plot of ES dynamics
  # -Boxplot of groups across time
  plotdf = data.frame(cbind(ES_list,sigTargets$time,ifelse(sigTargets$outcome==0,"Healthy","Sick")),stringsAsFactors = F)
  colnames(plotdf) = c("ES","TIME","OUTCOME")
  plotdf$ES = as.numeric(plotdf$ES) 
  plotdf$TIME = as.numeric(plotdf$TIME) * -1
  
  new_plotdf = c()
  
  for(t in unique(plotdf$TIME)){
    
    tsp_df = filter(plotdf,TIME==t)
    tsp_df$OUTCOME = as.factor(tsp_df$OUTCOME)
    ttestRes = t.test(tsp_df$ES~tsp_df$OUTCOME)
    if(ttestRes$p.value<=TRSH){
      pval = "significant"
      tsp_df = cbind(tsp_df,pval)
    }else{
      pval = "not-significant"
      tsp_df = cbind(tsp_df,pval)
    }
    
    new_plotdf = rbind(new_plotdf,tsp_df)
  }
  
  new_plotdf$pval = factor(new_plotdf$pval,levels = c("significant","not-significant"))
  new_plotdf$TIME = as.numeric(new_plotdf$TIME)
  
  ESplot = ggplot(new_plotdf, aes(x=TIME, y=ES, color=OUTCOME, shape = pval)) + geom_point(size=2.5) + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ggtitle(main_title) + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
  print(ESplot)
  
  new_plotdf$TIME = as.factor(new_plotdf$TIME)
  gp = ggplot(new_plotdf, aes(x = TIME, y = ES)) + geom_boxplot(aes(fill = OUTCOME,linetype = pval), alpha = 0.5) + ylab("Enrichment Score") + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm")) + ggtitle(main_title)
  print(gp)
  
}

# This function is similar to plotES, nonetheless, a significance test is not performed (boxplot option)

plotgenetrajectorybp = function(EXPRMAT,STARGETS,GENES){
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    genedframe$time = factor(genedframe$time,levels=c(12,9,6,3,0))
    
    for(GENE in GENES){
      
      gp = ggplot(genedframe, aes(x = time, y = genedframe[,GENE])) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + 
        stat_summary(fun.y=median, geom="line", aes(group=outcome, col=outcome),size=1.5)  + labs(title = GENE) + ylab("Expression") +
        theme(axis.title.x = element_text(size = 14,face="bold"), axis.title.y = element_text(size = 14,face="bold"),axis.text.x  = element_text(size=14),axis.text.y  = element_text(size=14))
      print(gp)
    }
    
  }
  else{
    print("We failed to find genes in expression matrix") 
    return(NULL)
  }
}


#
#
# HEATMAPS
#
#

plotHMps = function(genes,expmat,targetdf,time.point){
  #
  # Plots time-specific heatmaps of a count matrix
  #
  # inputs:
  # -genes: a list of genes to be clustered and plotted
  # -expmat: a count matrix that contains the subset of data wanted to be plotted
  # -
  color = ifelse(targetdf$outcome == 0,"lightcoral","lightblue")
  targetdf = cbind(targetdf,color)
  
  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("darkred", "lightgrey", "darkblue"))(n = 299)
  
  # Divide the matrix in time points
  if(is.null(time.point)){
    for(t in sort(unique(targetdf$time))){
      ttargets = filter(targetdf,time==t)
      tmat = expmat[genes,as.character(ttargets$sample_mask_id)]
      
      heatmap.2(tmat,
                main=paste("Time ",as.character(t*-1),sep=""),
                Rowv = TRUE,
                Colv = TRUE,
                distfun = dist,
                hclustfun = hclust,
                trace = "none",
                sepwidth=c(0.0005,0.0005),
                dendrogram = "column",
                col=my_palette,       # use on color palette defined earlier
                ColSideColors = as.character(ttargets$color)
      )
      par(xpd=TRUE)
      legend("bottomleft",      # location of the legend on the heatmap plot
             legend = c("control","case"), # category labels
             col = c("lightcoral","lightblue"),  # color key
             lty= 1,             # line style
             lwd = 10,            # line width
             inset=c(-0.1,0)
      )
      
    }
  }else{
    for(t in time.point){
      #t = time.point
      ttargets = filter(targetdf,time==t)
      tmat = expmat[genes,as.character(ttargets$sample_mask_id)]
      
      heatmap.2(tmat,
                main=paste("Time ",as.character(t*-1),sep=""),
                Rowv = TRUE,
                Colv = TRUE,
                distfun = dist,
                hclustfun = hclust,
                trace = "none",
                dendrogram = "column",
                col=my_palette,       # use on color palette defined earlier
                ColSideColors = as.character(ttargets$color)
      )
      par(xpd=TRUE)
      legend("bottomleft",      # location of the legend on the heatmap plot
             legend = c("Control","Case"), # category labels
             col = c("lightcoral","lightblue"),  # color key
             lty= 1,             # line style
             lwd = 10,           # line width
             inset=c(-0.1,0)
      )
      
    }
  }
}


#
#
# MDS plots
#
#

# TEDDY's data is incomplete, so it is impossible to calculate the distance 
# among all patients. What I do, is to heuristically determine the biggest
# distance matrix that can be calculated with the data given and perform an MDS with that

# This functions use patient matrices (rows= patients, columns = features)

#This only uses different colors for the outcome column
plot_MDS_Tv = function(meltMat,targetDF,main){
  
  #dist_pat = as.matrix(daisy(meltMat,metric = "gower"))
  dist_pat = as.matrix(dist(meltMat))
  dist_pat_temp = dist_pat
  test_coef = as.numeric(names(table(rowSums(is.na(dist_pat)))))
  
  maxSamples = c()
  for(coef in test_coef){
    dist_pat = dist_pat_temp
    dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef]
    dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1]
    maxSamples = c(maxSamples,dim(dist_pat)[1])
  }
  
  coef = test_coef[which(maxSamples == max(maxSamples))]
  dist_pat = dist_pat_temp
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef]
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1]
  
  mds = cmdscale(as.dist(dist_pat))
  
  pat_annotation = as.matrix(unique(targetDF[,c("mask_id","outcome")]))
  rownames(pat_annotation) = pat_annotation[,1] 
  pat_annotation = pat_annotation[,-1]
  pat_annotation = ifelse(pat_annotation==1,"case","control")
  pat_color = ifelse(pat_annotation=="case","red","blue")
  
  
  plot(mds, type = 'p',pch=21,main=main,col = pat_color[rownames(mds)],bg = pat_color[rownames(mds)])
  #text(mds[, 1], mds[, 2], pat_annotation[rownames(mds)],col = pat_color[rownames(mds)])
  
}

#This can show you different classes 
plot_MDS_Tv_v2 = function(meltMat,targetDF,main,GroupColumn){
  
  #dist_pat = as.matrix(daisy(meltMat,metric = "gower"))
  dist_pat = as.matrix(dist(meltMat))
  dist_pat_temp = dist_pat
  test_coef = as.numeric(names(table(rowSums(is.na(dist_pat)))))
  
  maxSamples = c()
  for(coef in test_coef){
    if(coef==0){
      dist_pat = dist_pat_temp
      dist_pat = dist_pat[rowSums(is.na(dist_pat))==0,rowSums(is.na(dist_pat))==0,drop=F]
      maxSamples = c(maxSamples,dim(dist_pat)[1])
      
    }else{
      dist_pat = dist_pat_temp
      dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef,drop=F]
      dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1,drop=F]
      maxSamples = c(maxSamples,dim(dist_pat)[1])
    }
    
  }
  
  coef = test_coef[which(maxSamples == max(maxSamples))][1]
  dist_pat = dist_pat_temp
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef]
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1]
  
  mds = cmdscale(as.dist(dist_pat))
  
  pat_annotation = as.matrix(unique(targetDF[,c("mask_id",GroupColumn)]))
  rownames(pat_annotation) = pat_annotation[,1] 
  pat_annotation = pat_annotation[,-1]
  
  clist = c("blue","red","green","purple","orange","yellow","black")
  GFact = unique(targetDF[,GroupColumn])
  clist = clist[1:length(GFact)]
  pat_color = pat_annotation
  
  for(i in 1:length(GFact)){
    pat_color[pat_color==GFact[i]] = clist[i]
  }
  
  plot(mds, type = 'n',pch=21,main=main,col = pat_color[rownames(mds)],bg = pat_color[rownames(mds)])
  text(mds[, 1], mds[, 2], pat_annotation[rownames(mds)],col = pat_color[rownames(mds)])
}

#
#
# LINEAR MODELS
#
#

#resultData is a GLM result from getSignatures
summaryTestTable = function(resultData){
  #Function that plots the coefficients and pvalues of a GLM Result list object
  PVAL_mat = matrix(1,nrow = nrow((resultData[[1]])), ncol = length(resultData))
  COEF_mat = matrix(1,nrow = nrow((resultData[[1]])), ncol = length(resultData))
  
  rownames(PVAL_mat) = rownames(COEF_mat) = resultData[[1]]$Hallmark
  colnames(PVAL_mat) = colnames(COEF_mat) = names(resultData)
  
  for(Tres in names(resultData)){
    
    COEF_mat[,Tres] = resultData[[Tres]]$Coefficient
    PVAL_mat[,Tres] = resultData[[Tres]]$adj_P_val
    
  }
  
  COEF_mat = round(COEF_mat,2)
  PVAL_mat = round(PVAL_mat,2)
  
  my_palette <- colorRampPalette(c("lightcoral", "white", "lightgreen"))(n = 20)
  
  PVALplot = heatmap.2(PVAL_mat,Rowv = FALSE,Colv = FALSE,dendrogram = 'none',cellnote = PVAL_mat,trace = "none",notecol = "black",col=my_palette,lwid=  c(0.05,.5),lhei= c(0.05,.30),margins(30,15),key = F,srtCol=35)
  COEFplot = heatmap.2(COEF_mat,Rowv = FALSE,Colv = FALSE,dendrogram = 'none',cellnote = COEF_mat,trace = "none",notecol = "black",col=my_palette,lwid=  c(0.05,.5),lhei= c(0.05,.30),margins(30,15),key = F,srtCol=35)
  
  print(PVALplot)
  print(COEFplot)
  
  return(list("pval" = PVAL_mat,"coef" = COEF_mat))
}

#
#
# Cluster
#
#

# This function requires a list of median responses (or any other quantification) of a given feature
#
#
#       time 12 \ time 9 \ ...
#GroupA
#GroupB
#GroupC
#

plot_CompareDynamics = function(MedianResps, main){
  # Plot the median/mean responses of any count
  plotDF = c()
  for(TG in names(MedianResps)){
    plotDF = rbind(plotDF, cbind(melt(MedianResps[[TG]]),TestGroup = TG))
  }
  colnames(plotDF) = c("ClusterName","Time","MedianES","TestGroup")
  plotDF = cbind(plotDF,Tclass = paste(plotDF$ClusterName,plotDF$TestGroup,sep = "_"))
  
  plotDF$Time =  as.numeric(gsub("time_","",plotDF$Time))
  plotDF$ClusterName = factor(plotDF$ClusterName)
  plotDF$MedianES = as.numeric(plotDF$MedianES)
  plotDF$TestGroup = factor(plotDF$TestGroup)
  plotDF$Tclass = factor(plotDF$Tclass)
  GEX_plot <- ggplot(plotDF, aes(x = plotDF[,"Time"] * -1, y = plotDF[,"MedianES"], group = Tclass, color = TestGroup)) + 
    scale_x_continuous(name="Time",breaks=sort(unique(plotDF[,"Time"] * -1))) + ylab("Median Expression") + labs(title = main) +
    theme(axis.title.x = element_text(size = 14,face="bold"), axis.title.y = element_text(size = 14,face="bold"),axis.text.x  = element_text(size=14),axis.text.y  = element_text(size=14)) +
    geom_point() + geom_line()
  print(GEX_plot)
  #return(GEX_plot)
}

# This function compares the dynamics of any feature for any groups you want to compare
plot_CompareDynamics_Boxplots = function(CountMatrix,TargetDF,comp_col,selFeat){
  # CountMatrix : rows=Individuals, columns = Features(genes,metabolites, etc)
  # TargetDF : Target data.frame (it must contain a column called "time" and a column specifying your groups)
  # comp_col : The name of the column that contains your groups
  # selFeat : Features you want to plot (metabolite names, genes, etc) (In a vector!!)
 
  #
  selFeatCounts = CountMatrix[selFeat,,drop=F]
  plotDF_str = cbind(TargetDF,t(selFeatCounts))
  
  plotDF_str[,comp_col] = factor(plotDF_str[,comp_col])
  plotDF_str$time = plotDF_str$time * -1
  plotDF_str$time = factor(plotDF_str$time)
  
  for(Feat in selFeat){
    
    gp = ggplot(plotDF_str, aes(x = time, y = plotDF_str[,Feat])) + geom_boxplot(aes(fill = agegroup), alpha = 0.5) + 
      ylab("Expression") + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm")) + 
      ggtitle(Feat)
    
    print(gp)
    
  }
  
}





