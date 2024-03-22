#!/bin/bash

# Define the URLs of the files to download
urls=(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/suppl/GSE152004%5F695%5Fexpr%5Fnorm%2Etxt%2Egz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/matrix/GSE152004_series_matrix.txt.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201955/suppl/GSE201955%5FRNAseq%5F118%5Fprocesseddata%2Etxt%2Egz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201955/matrix/GSE201955-GPL20301_series_matrix.txt.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58434/suppl/GSE58434%5FAll%5FSample%5FFPKM%5FMatrix%2Etxt%2Egz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58434/matrix/GSE58434_series_matrix.txt.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67472/matrix/GSE67472_series_matrix.txt.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69683/matrix/GSE69683_series_matrix.txt.gz"
)

# Loop through the URLs and download each file
for url in "${urls[@]}"; do
    wget "$url"
done
