# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load required packages
library(DESeq2)
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

# define function
perform_WGCNA_analysis <- function(count_data_path, metadata_path, normalize_with_deseq2, modules_df_path, r_and_Pvalue_data_path) {
  # Load Data and Metadata
  expression_data <- read.csv(count_data_path, header = TRUE, row.names = 1)
  metadata <- read.csv(metadata_path, row.names = 1)
  row.names(metadata) <- colnames(expression_data)
  
  # expression_data[is.na(expression_data)] <- 0
  
  if (normalize_with_deseq2) {
    dds <- DESeqDataSetFromMatrix(countData = round(expression_data), colData = metadata, design = ~Type)
    dds <- DESeq(dds)
    a <- varianceStabilizingTransformation(dds)
    b <- getVarianceStabilizedData(dds)
    b1 <- rowVars(b)
    q00_wpn <- quantile(rowVars(b))
    expression_data <- b[b1 > q00_wpn, ]
  }
  
  
  # Compute coefficient of variation (CV) of each gene
  expression_data <- data.frame(expression_data)
  expression_data$CV <- apply(expression_data, 1, function(x) sd(x) / mean(x) * 100)
  
  # Select genes based on CV greater than 4%
  expression_data <- expression_data[expression_data$CV > 4, ]
  expression_data <- subset(expression_data, select = -CV)
  
  # Perform WGCNA
  cor <- WGCNA::cor
  input_mat <- t(expression_data)
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)
  picked_power <- sft$powerEstimate
  netwk <- blockwiseModules(input_mat,
                            power = picked_power,
                            networkType = "signed",
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minModuleSize = 50,
                            maxBlockSize = 5000,
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "ER",
                            numericLabels = TRUE,
                            verbose = 3)
  
  # Write module data to file
  module_df <- data.frame(
    gene_id = names(netwk$colors),
    colors = labels2colors(netwk$colors)
  )
  write.csv(module_df, modules_df_path)
  
  # Get Module Eigengenes per cluster
  nGenes = ncol(input_mat);
  nSamples = nrow(input_mat);
  mergedColors = labels2colors(netwk$colors)
  MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
  MEList = moduleEigengenes(input_mat, colors = mergedColors)
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)
  
  # Read clinical traits of asthma and control subjects
  metadata$Type <- as.numeric(factor(metadata$Type))
  metadata <- subset(metadata, select = Type)
  
  moduleTraitCor <- cor(MEs, metadata$Type, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  colnames(moduleTraitCor) <- c("Correlation")
  colnames(moduleTraitPvalue) <- c("Pvalue")
  r_and_Pvalue_data <- cbind(moduleTraitCor, moduleTraitPvalue)
  write.csv(r_and_Pvalue_data, r_and_Pvalue_data_path)
}


# Define a list of parameter combinations
parameter_combinations <- list(
  list(
    count_data_path = "data1/GSE152004_norm_data.csv",
    metadata_path = "data1/GSE152004_metadata.csv",
    normalize_with_deseq2 = TRUE,
    modules_df_path = "data2/GSE152004_modules.csv",
    r_and_Pvalue_data_path = "data2/GSE152004_r_and_Pvalue_data.csv"
  ),
  list(
    count_data_path = "data1/GSE201955_norm_data.csv",
    metadata_path = "data1/GSE201955_metadata.csv",
    normalize_with_deseq2 = FALSE,
    modules_df_path = "data2/GSE201955_modules.csv",
    r_and_Pvalue_data_path = "data2/GSE201955_r_and_Pvalue_data.csv"
  ),
  list(
    count_data_path = "data1/GSE67472_norm_data.csv",
    metadata_path = "data1/GSE67472_metadata.csv",
    normalize_with_deseq2 = FALSE,
    modules_df_path = "data2/GSE67472_modules.csv",
    r_and_Pvalue_data_path = "data2/GSE67472_r_and_Pvalue_data.csv"
  ),
  list(
    count_data_path = "data1/GSE58434_norm_data.csv",
    metadata_path = "data1/GSE58434_metadata.csv",
    normalize_with_deseq2 = TRUE,
    modules_df_path = "data2/GSE58434_modules.csv",
    r_and_Pvalue_data_path = "data2/GSE58434_r_and_Pvalue_data.csv"
  ),
  list(
    count_data_path = "data1/GSE69683_norm_data.csv",
    metadata_path = "data1/GSE69683_metadata.csv",
    normalize_with_deseq2 = FALSE,
    modules_df_path = "data2/GSE69683_modules.csv",
    r_and_Pvalue_data_path = "data2/GSE69683_r_and_Pvalue_data.csv"
  )
  # Add more parameter combinations as needed
)

# Iterate over parameter combinations and run differential gene expression analysis
for (params in parameter_combinations) {
  perform_WGCNA_analysis(
    count_data_path = params$count_data_path,
    metadata_path = params$metadata_path,
    normalize_with_deseq2 = params$normalize_with_deseq2,
    modules_df_path = params$modules_df_path,
    r_and_Pvalue_data_path = params$r_and_Pvalue_data_path
  )
}



