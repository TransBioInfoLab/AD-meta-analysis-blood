
# Introduction 

This repository contains the analysis R code for Integrative meta-analysis of epigenome-wide association studies in Alzheimer’s disease. Some of the analysis were made using SAS were not included in this repository.

## Article

Integrative meta-analysis of epigenome-wide association studies identifies genomic and epigenomics differences in the brain and the blood in Alzheimer’s disease

## Authors

- Tiago C Silva
- Juan I. Young
- Lanyu Zhang
- Lissette Gomez
- Michael A. Schmidt
- Achintya Varma
- Xi Chen
- Eden R. Martin
- Lily Wang


# Folders in this repository

DRAFT-TABLES_FIGURES_4-17-2021: Has tables used in the code
code: main folder with code analysis
analysis_results: Folder with some of the analysis results

# Code

You should start processing the data in `code/ADNI` and `code/AIBL`.


## Install required packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.13",ask = FALSE) # Install last version of Bioconductor

list.of.packages <- c(
  "readr",
  "readxl",
  "plyr",
  "dplyr",
  "tidyr",
  "GenomicRanges",
  "SummarizedExperiment",
  "myGene",
  "ggpubr",
  "stats",
  "gt",
  "fgsea",
  "minfi",
  "bacon",
  "wateRmelon",
  "DMRcate",
  "lumi",
  "RPMM",
  "sm",
  "doParallel",
  "EpiDISH",
  "GEOquery",
  "ggrepel",
  "DT",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "RVenn",
  "GWASTools",
  "meta",
  "metap"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr")
```

# Platform information

```
version  R version 4.1.0 (2021-05-18)
 os       macOS Big Sur 11.4          
 system   x86_64, darwin17.0          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       America/New_York            
 date     2021-07-12      
```
