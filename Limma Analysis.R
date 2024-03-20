# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load required packages
library(edgeR)
library(DESeq2)
library(dendextend)
library(heatmap3)
library(limma)
library(dplyr)
library(tidyverse)
library(scales)
library(sva)
library(ggrepel)
library(patchwork)

############################################################
### Load Data and Metadata
############################################################
expression_data <- read.csv("data1/GSE152004_norm_data.csv", header=T, row.names=1)
metadata<-read.csv("data1/GSE152004_metadata.csv", row.names = 1)

############################################################
### Filtering Dataset Based on Coefficient of Variation
############################################################
# compute  coefficient of variation(CV) of each gene
expression_data$CV<-apply(expression_data,1, function(x) sd(x) / mean(x) * 100) 
# Select genes based on CV greater than 4 % 
expression_data<-expression_data[expression_data$CV > 4,]
expression_data <- subset(expression_data, select = -CV)

#Normalize with DESeq using variance stabilizingg transfomation
expression_data[is.na(expression_data)] <- 0
dds <- DESeqDataSetFromMatrix(countData=round(expression_data), colData = metadata , design = ~Type)
dds <- DESeq(dds)
a <- varianceStabilizingTransformation(dds)
b <- getVarianceStabilizedData(dds)
b1 <- rowVars(b)
q00_wpn <- quantile( rowVars(b)) 
expression_data <- b[ b1 > q00_wpn, ]

# Use Limma package for DEG analysis for normalized GEO dataset 
contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('asthmatic', 'healthy control'), "Diff")
sample<-metadata$Type
design.mat<-model.matrix(~0+sample)
colnames(design.mat)<-levels(sample)
fit<-lmFit(expression_data, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" )
DEG1 <- as.data.frame(deg1)
write.csv(DEG1, "data2/GSE67472_Lima_Results.csv")

