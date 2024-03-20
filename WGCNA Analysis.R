# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load required packages
library(DESeq2)
library(edgeR)
library(WGCNA)
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
### Data processing and selection of GSE152004 Data for DEGs.
############################################################

### Load Raw count GSE152004 RNA-seq matrix
GSE152004_count<- read.csv("data/GSE152004_count_data.csv", header=T, row.names=1)
head(GSE152004_count)
GSE152004_metadata<-read.csv("data/GSE152004_metadata.csv", row.names = 1)

#Normalize Counts with DESeq using variance stabilizingg transfomation
dds <- DESeqDataSetFromMatrix(countData=GSE152004_count, colData = GSE152004_metadata , design = ~Type)
dds <- DESeq(dds)
a <- varianceStabilizingTransformation(dds)
b <- getVarianceStabilizedData(dds)
b1 <- rowVars(b)
summary(b1)
q00_wpn <- quantile( rowVars(b)) 
expr_normalized_Data <- b[ b1 > q00_wpn, ]
xx<-expr_normalized_Data
dim(xx)
write.csv(xx, "data/GSE152004_Normalized_data_matrix.csv")

############################################################
### Filtering GSE152004 Dataset Based on Coefficient of Variation
############################################################
# Load normalized gene expression data of GSE152004
NormalizedGSE152004<-read.csv("data/GSE152004_Normalized_data_matrix.csv", 
                              row.names = 1)
# compute  coefficient of variation(CV) of each gene
NormalizedGSE152004$CV<-apply(NormalizedGSE152004,1, function(x) sd(x) / mean(x) * 100) 

# Select genes based on CV greater than 4 % 
NormalizedGSE152004<-NormalizedGSE152004[NormalizedGSE152004$CV > 4,]
NormalizedGSE152004 <- subset(NormalizedGSE152004, select = -CV)
write.csv(NormalizedGSE152004, "data/FilteredNormalizedGSE152004.csv")


############################################################
### #WGCNA analysis on GSE152004 Dataset 
############################################################
dat<-read.csv("data/FilteredNormalizedGSE152004.csv", row.names = 1) 
input_mat = t(dat)
# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


picked_power = 5
temp_cor <- cor       
cor <- WGCNA::cor         
netwk <- blockwiseModules(input_mat,
                          power = picked_power,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 50,
                          maxBlockSize = 5000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3)

module_eigengenes <- netwk$MEs
netwk$MEs
# Print out a preview
head(module_eigengenes)
dim(module_eigengenes)
write.csv(module_eigengenes,"data/GSE152004_module_eigengeneMerge.csv")

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df
write.csv(module_df, "data/GSE152004_modules.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

module_order


MEs0
write.csv(MEs0, "data/GSE152004_MEs0_GSEMerge.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
MEList = moduleEigengenes(input_mat, colors = mergedColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs0, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

#**************************************************
# Define numbers of genes and samples
nGenes = ncol(input_mat);
nSamples = nrow(input_mat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
dim(MEs)
MEs$MEgreen

# Read clinical traits of asthma and control subjects
clinical_traits = read.csv("data/GSE152004_metadata.csv", row.names = 1)
clinical_traits$Type <- as.numeric(factor(clinical_traits$Type))
clinical_traits <- subset(clinical_traits, select = Type)
rownames(clinical_traits)
dim(clinical_traits)

# sample names should be consistent in eigen genes and traits !!!!
table(rownames(MEs) == rownames(clinical_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, clinical_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
colnames(moduleTraitCor)[colnames(moduleTraitCor) == 'Type'] <- 'Correlation'
colnames(moduleTraitPvalue)[colnames(moduleTraitPvalue) == 'Type'] <- 'Pvalue'
r_and_Pvalue_data = cbind(moduleTraitCor,moduleTraitPvalue)
head(r_and_Pvalue_data)
write.csv(r_and_Pvalue_data, "data/GSE152004_r_and_Pvalue_data.csv")

xx2<-moduleTraitPvalue
moduleTraitCor
xx3<-moduleTraitCor

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


############################################################
### Differential Gene expression Analysis on GSE152004 Dataset 
############################################################
# For count data, we used DESeq2 package for DEGs analysis
GSE152004_count_data <- read.csv("data/GSE152004_count_data.csv", 
                                 header=T, row.names=1)
GSE152004_metadata <-read.csv("data/GSE152004_metadata.csv", row.names = 1)

#Create a Datasheet from count matrix and labels
dds <- DESeqDataSetFromMatrix(GSE152004_count_data ,GSE152004_metadata , design = ~Type)

# healthy control->set reference control subjects
dds$Type <- relevel(dds$Type, ref = "healthy control") 
dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
write.csv(res, "data/GSE152004_Deseq2_Results.csv")





############################################################
### Filtering GSE67472 Dataset Based on Coefficient of Variation
############################################################
# Load normalized gene expression data of GSE67472
NormalizedGSE67472<-read.csv("data/GSE67472_normalized_data.csv", row.names = 1)
# compute  coefficient of variation(CV) of each gene
NormalizedGSE67472$CV<-apply(NormalizedGSE67472,1, function(x) sd(x) / mean(x) * 100) 

# Select genes based on CV greater than 4 % 
NormalizedGSE67472<-NormalizedGSE67472[NormalizedGSE67472$CV > 4,]
NormalizedGSE67472 <- subset(NormalizedGSE67472, select = -CV)
write.csv(NormalizedGSE67472, "data/FilteredNormalizedGSE67472.csv")


############################################################
### #WGCNA analysis on GSE67472 Dataset 
############################################################
dat<-read.csv("data/FilteredNormalizedGSE67472.csv", row.names = 1) 
input_mat = t(dat)
# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


picked_power = 7
temp_cor <- cor       
cor <- WGCNA::cor         
netwk <- blockwiseModules(input_mat,
                          power = picked_power,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 50,
                          maxBlockSize = 5000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3)

module_eigengenes <- netwk$MEs
netwk$MEs
# Print out a preview
head(module_eigengenes)
dim(module_eigengenes)
write.csv(module_eigengenes,"data/GSE67472_module_eigengeneMerge.csv")

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df
write.csv(module_df, "data/GSE67472_modules.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

module_order


MEs0
write.csv(MEs0, "data/GSE67472_MEs0_GSEMerge.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
MEList = moduleEigengenes(input_mat, colors = mergedColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs0, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

#**************************************************
# Define numbers of genes and samples
nGenes = ncol(input_mat);
nSamples = nrow(input_mat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
dim(MEs)
MEs$MEgreen

# Read clinical traits of asthma and control subjects
clinical_traits = read.csv("data/GSE67472_metadata.csv", row.names = 1)
clinical_traits$Type <- as.numeric(factor(clinical_traits$Type))
clinical_traits <- subset(clinical_traits, select = Type)
rownames(clinical_traits)
dim(clinical_traits)

# sample names should be consistent in eigen genes and traits !!!!
table(rownames(MEs) == rownames(clinical_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, clinical_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
colnames(moduleTraitCor)[colnames(moduleTraitCor) == 'Type'] <- 'Correlation'
colnames(moduleTraitPvalue)[colnames(moduleTraitPvalue) == 'Type'] <- 'Pvalue'
r_and_Pvalue_data = cbind(moduleTraitCor,moduleTraitPvalue)
head(r_and_Pvalue_data)
write.csv(r_and_Pvalue_data, "data/GSE67472_r_and_Pvalue_data.csv")

xx2<-moduleTraitPvalue
moduleTraitCor
xx3<-moduleTraitCor

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Limma package was used for DEG analysis for normalized GEO dataset 
GSE67472_normalized_data<-read.csv("data/GSE67472_normalized_data.csv", row.names = 1)
GSE67472_metadata<-read.csv("data/GSE67472_metadata.csv")
GSE67472_metadata<-GSE67472_metadata[,-1]
dim(GSE67472_metadata)

contrast.mat<-matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat)<-list(c('asthmatic', 'healthy control'), "Diff")
sample<-GSE67472_metadata$Type
design.mat<-model.matrix(~0+sample)
colnames(design.mat)<-levels(sample)
design.mat
dim(design.mat)
fit<-lmFit(GSE67472_normalized_data, design.mat)
fit2<-contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
deg<-topTable(fit3) # top ranked DEGs
deg1 <- topTable(fit3, n=Inf, coef=1,adjust.method="BH" )
DEG1 <- as.data.frame(deg1)
View(DEG1)
write.csv(DEG1, "data/GSE67472_Lima_Results.csv")

# # Validation datasets having batch effect were corrected using SVA package 
# datx<-read.csv("data with batch_effect.csv", row.names = 1)
# b<-read.csv("Design_matrix.csv")
# b$batch<-as.factor(b$batch)
# dim(b)
# dim(datx)
# batch = b$batch
# batch
# # parametric adjustment
# combat_gdata = ComBat(dat= datx, batch= batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# combat_gdata
# write.csv(combat_gdata, "Batch_adjusted_by_sva.csv")
# # non-parametric adjustment, mean-only version
# combat_edata2 = ComBat(dat=datx, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TR)
# combat_edata2 
# write.csv(combat_gdata2, "Batch_adjusted_by_sva.csv")

