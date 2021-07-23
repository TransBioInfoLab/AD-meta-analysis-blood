
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


# Code




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
  "mygene",
  "ggpubr",
  "stats",
  "gt",
  "fgsea",
  "minfi",
  "bacon",
  "wateRmelon",
  "missMethyl",
  "DMRcate",
  "lumi",
  "RPMM",
  "sm",
  "MethReg",
  "writexl",
  "readxl",
  "gridExtra",
  "doParallel",
  "ReMapEnrich",
  "EpiDISH",
  "GEOquery",
  "ggrepel",
  "DT",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "RVenn",
  "GWASTools",
  "meta",
  "metap",
  "ExperimentHub",
  "lubridate"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr")
```

For ADNIMERGE, download it from https://ida.loni.usc.edu/: Merged ADNI 1/GO/2 Packages for R

```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```

## Analysis
### Study cohorts, Preprocessing of DNA methylation data, Single Cohort analysis


| File                 | Dataset | Link |
|----------------------|-------------|-------------|
| code/ADNI/ADNI_SAS.Rmd        |   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/ADNI/ADNI_SAS.Rmd) |
| code/ADNI/ADNI_DMR_analysis_SAS.Rmd        |   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/ADNI/ADNI_DMR_analysis_SAS.Rmd) |
| code/AIBL/AIBL.Rmd           |   AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/AIBL/AIBL.Rmd) |
| code/AIBL/AIBL_DMR.Rmd           |   AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/AIBL/AIBL_DMR.Rmd) |
| code/Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals.R          |   ADNI    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals.R) |
| code/Clinical/clinical_info.Rmd          |   ADNI & AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Clinical/clinical_info.Rmd) |


### Blood samples meta-analysis

| File                 | Link |
|----------------------|-------------|
| code/meta-analysis/meta-analysis-two-cohorts_glm.Rmd       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/meta-analysis/meta-analysis-two-cohorts_glm.Rmd) |

### Cross-tissue meta-analysis	

| File                 | Link |
|----------------------|-------------|
| code/cross_tissue_meta_analysis/cross_tissue_meta_analysis.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis.Rmd) |
| code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_DMR.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_DMR.Rmd) |


### Functional annotation of significant methylation differences

| File                 | Link |
|----------------------|-------------|
| code/annotations/create_great_annotation.R     |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/annotations/create_great_annotation.R) |
| code/annotations/annotate_enhancer.R    |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/annotations/annotate_enhancer.R) |

### Correlations between methylation levels of significant CpGs and DMRs in AD with expressions of nearby genes

| File                 | Link |
|----------------------|-------------|
| code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_DMR.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_DMR.R) |
| code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg.R) |
| code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_DMR.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_DMR.R) |
| code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg.R) |


### MethReg integrative analysis


| File                 | Link |
|----------------------|-------------|
| code/MethReg/Blood_MethReg_DMR_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_cpg.R) |
| code/MethReg/Blood_MethReg_DMR_median.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_median.R) |
| code/MethReg/Blood_MethReg_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_cpg.R) |
| code/MethReg/Brain_MethReg_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Brain_MethReg_cpg.R) |



### Integrative analysis of DNA methylation differences in the brain and blood with transcriptome-wide gene expressions

| File                 | Link |
|----------------------|-------------|
| code/TWAS_pathway_analysis/TWAS_approach_blood.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/TWAS_approach_blood.R) |
| code/TWAS_pathway_analysis/TWAS_approach_brain.R       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/TWAS_approach_brain.R) |
| code/TWAS_pathway_analysis/jaccard_idx.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/jaccard_idx.R) |


### Correlation and overlap with genetic susceptibility loci

| File                 | Link |
|----------------------|-------------|
| code/Overlap_with_AD_associated_genetics_loci/Overlap_with_AD_associated_genetics_loci.R       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Overlap_with_AD_associated_genetics_loci/Overlap_with_AD_associated_genetics_loci.R) |


### Correlation of AD-associated CpGs and DMRs methylation levels in blood and brain samples	

| File                 | Link |
|----------------------|-------------|
| code/MethReg/Blood_MethReg_DMR_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_cpg.R) |




# Platform information

```r
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
