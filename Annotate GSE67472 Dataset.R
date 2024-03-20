# load required packages
library(Biobase)
library(GEOquery)
library(annotate)
library(hgu133plus2.db)
library(pd.hg.u133.plus.2)
library(jsonlite)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################################################
########## Annotate GSE67472
######################################################################
# Load data GSE67472
data<-read.csv("data/GSE67472_normalized_data.csv", row.names = 1) 
gene_ids = rownames(data)
aliasgenids = select(hgu133plus2.db, gene_ids, c("SYMBOL")) 
aliasgenids <- na.omit(aliasgenids)
# Convert the table to a named list
data_list <- as.list(setNames(aliasgenids$SYMBOL, aliasgenids$PROBEID))
# Convert the list to JSON format
json_data <- toJSON(data_list, auto_unbox = TRUE)
# Save the JSON data to a file
write(json_data, file = "GSE67472_annotation.json")