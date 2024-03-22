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

# define function
perform_DGE_analysis <- function(input_expression_data_path, 
                                 metadata_path, 
                                 output_path, 
                                 normalize_with_deseq2) {
  # Load Data and Metadata
  expression_data <- read.csv(input_expression_data_path, header = TRUE, row.names = 1)
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
  
  # Use Limma package for DEG analysis for normalized GEO dataset 
  contrast.mat <- matrix(c(1, -1), ncol = 1)
  dimnames(contrast.mat) <- list(c('asthmatic', 'healthy control'), "Diff")
  sample <- metadata$Type
  design.mat <- model.matrix(~0 + sample)
  colnames(design.mat) <- levels(sample)
  fit <- lmFit(expression_data, design.mat)
  fit2 <- contrasts.fit(fit, contrast.mat)
  fit3 <- eBayes(fit2)
  deg1 <- topTable(fit3, n = Inf, coef = 1, adjust.method = "BH")
  DEG1 <- as.data.frame(deg1)
  
  # Write results to output_path
  write.csv(DEG1, output_path)
}

# Define a list of parameter combinations
parameter_combinations <- list(
  list(
    input_expression_data_path = "data1/GSE201955_norm_data.csv",
    metadata_path = "data1/GSE201955_metadata.csv",
    output_path = "data2/GSE201955_Limma_Results.csv",
    normalize_with_deseq2 = FALSE
  ),
  list(
    input_expression_data_path = "data1/GSE152004_norm_data.csv",
    metadata_path = "data1/GSE152004_metadata.csv",
    output_path = "data2/GSE152004_Limma_Results.csv",
    normalize_with_deseq2 = TRUE
  ),
  list(
    input_expression_data_path = "data1/GSE67472_norm_data.csv",
    metadata_path = "data1/GSE67472_metadata.csv",
    output_path = "data2/GSE67472_Limma_Results.csv",
    normalize_with_deseq2 = FALSE
  ),
  list(
    input_expression_data_path = "data1/GSE58434_norm_data.csv",
    metadata_path = "data1/GSE58434_metadata.csv",
    output_path = "data2/GSE58434_Limma_Results.csv",
    normalize_with_deseq2 = TRUE
  ),
  list(
    input_expression_data_path = "data1/GSE69683_norm_data.csv",
    metadata_path = "data1/GSE69683_metadata.csv",
    output_path = "data2/GSE69683_Limma_Results.csv",
    normalize_with_deseq2 = FALSE
  )
  # Add more parameter combinations as needed
)

# Iterate over parameter combinations and run differential gene expression analysis
for (params in parameter_combinations) {
  perform_DGE_analysis(
    input_expression_data_path = params$input_expression_data_path,
    metadata_path = params$metadata_path,
    output_path = params$output_path,
    normalize_with_deseq2 = params$normalize_with_deseq2
  )
}
