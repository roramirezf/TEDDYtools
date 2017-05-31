############################
# TEDDY tools: Requirements
# Author: Ricardo Ramirez
# 
# Loads all the functions, objects and collection of functions of TEDDYtools
#
#
############################

parent = getwd()

# The location of TEDDYtoolsV2
setwd("/home/rramirez/Dropbox (UFL)/myTEDDY_V2/TEDDYtoolsV2/")

# Global Variables : RAW DATA

# Gene Expression

load("GlobalData/GEXtargets.ro") #Target Matrix
load("GlobalData/rawTEDDY_GEX.ro") #Count Matrix (Original Table)
load("GlobalData/rawTEDDY_GEX_inversevst.ro") #New Matrix

# Metabolomics

load("GlobalData/MetabolicTargets.ro") #Target Matrix
load("GlobalData/Metabolic_Counts_Raw.ro") #Complete Matrix
load("GlobalData/GCTOFRAW_counts.ro") #GCTOF
load("GlobalData/negLipRAW_counts.ro") #NegLip
load("GlobalData/posLipRAW_counts.ro") #PosLip

# GeneSets

load("GlobalData/TEDDY_geneSets.ro")

# External Libraries
library(dplyr)
library(data.table)
library(illuminaHumanv4.db)
library(GO.db)
library(splines)
library(limma)
library(piano)
library(matrixStats)
library(maSigPro)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(cluster)
library(gtools)
library(stringdist)
library(RDAVIDWebService)
library(fpc)
library(mclust)
library(scales)

# Load TEDDYtools Suite

print("Loading Processing and filtering tools")
source("Suite/Processing_Filtering_TT.R")
print("Loading Count Manipulation tools")
source("Suite/CountManipulation_TT.R")
print("Loading GSEA tools")
source("Suite/GSEA_TT.R")
print("Loading Linear Modelling tools")
source("Suite/LinearModels_DEA_TT.R")
print("Loading Visual Manager tools")
source("Suite/VisualManager_TT.R")
print("Loading Transcriptional Signature Comparison tools")
source("Suite/TranscriptionalSignComp_TT.R")
print("Loading Clustering tools")
source("Suite/Clustering_TT.R")


#Return to parent directory
setwd(parent)
