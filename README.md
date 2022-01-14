
# Cross-tissue meta-analysis of blood and brain epigenome-wide association studies in Alzheimer’s disease 
Tiago C. Silva, Juan I. Young, Lanyu Zhang, Lissette Gomez, Michael A. Schmidt, Achintya Varma, Xi Chen, Eden R. Martin, Lily Wang

### Description

This github repository includes scripts used for the analyses in the above manuscript. 

In this work, we performed a meta-analysis of two large independent blood-based epigenome-wide association studies, the ADNI and AIBL studies, and identified 5 CpGs mapped to the SPIDR, CDH6 genes, and intergenic regions that were significantly associated with AD diagnosis. Furthermore, to identify blood-based DNA methylation markers that also change with underlying neuropathology in the brain, we next performed a cross-tissue meta-analysis by combining these blood DNA methylation datasets with four additional DNA methylation datasets, which included a total of  1030 brain prefrontal cortex samples. Our findings provide a useful resource for future biomarker studies in AD.  

### 1. Study cohorts, Preprocessing of DNA methylation data, Single Cohort analysis


| File                 | Dataset | Link |
|----------------------|-------------|-------------|
| code/ADNI/ADNI_SAS.Rmd        |   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/ADNI/ADNI_SAS.Rmd) |
| code/ADNI/ADNI_DMR_analysis_SAS.Rmd        |   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/ADNI/ADNI_DMR_analysis_SAS.Rmd) |
|code/ADNI/GLMM_models_ADNIdata_all.sas|   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/ADNI/GLMM_models_ADNIdata_all.sas)  |    
| code/AIBL/AIBL.Rmd           |   AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/AIBL/AIBL.Rmd) |
| code/AIBL/AIBL_DMR.Rmd           |   AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/AIBL/AIBL_DMR.Rmd) |
| code/Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals.R          |   ADNI    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals.R) |
| code/Clinical/clinical_info.Rmd          |   ADNI & AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Clinical/clinical_info.Rmd) |


### 2. Blood samples meta-analysis

| File                 | Link |
|----------------------|-------------|
| code/meta-analysis/meta-analysis-two-cohorts_glm.Rmd       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/meta-analysis/meta-analysis-two-cohorts_glm.Rmd) |
| code/meta-analysis/meta-analysis-two-cohorts_glm_DMR.Rmd      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/meta-analysis/meta-analysis-two-cohorts_glm_DMR.Rmd) |

 
### 3. Cross-tissue meta-analysis	

| File                 | Link |
|----------------------|-------------|
| code/cross_tissue_meta_analysis/cross_tissue_meta_analysis.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis.Rmd) |
| code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_DMR.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_DMR.Rmd) |


### 4. Functional annotation of significant methylation differences

| File                 | Link |
|----------------------|-------------|
| code/annotations/create_great_annotation.R     |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/annotations/create_great_annotation.R) |
| code/annotations/annotate_enhancer.R    |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/annotations/annotate_enhancer.R) |

### 5. Correlations between methylation levels of significant CpGs and DMRs in AD with expressions of nearby genes

| File                 | Link |
|----------------------|-------------|
| code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_DMR.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_DMR.R) |
| code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg.R) |
| code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_DMR.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_DMR.R) |
| code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg.R) |


### 6. MethReg integrative analysis


| File                 | Link |
|----------------------|-------------|
| code/MethReg/Blood_MethReg_DMR_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_cpg.R) |
| code/MethReg/Blood_MethReg_DMR_median.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_median.R) |
| code/MethReg/Blood_MethReg_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_cpg.R) |
| code/MethReg/Brain_MethReg_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Brain_MethReg_cpg.R) |



### 7. Integrative analysis of DNA methylation differences in the brain and blood with transcriptome-wide gene expressions

| File                 | Link |
|----------------------|-------------|
| code/TWAS_pathway_analysis/TWAS_approach_blood.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/TWAS_approach_blood.R) |
| code/TWAS_pathway_analysis/TWAS_approach_brain.R       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/TWAS_approach_brain.R) |
| code/TWAS_pathway_analysis/jaccard_idx.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/TWAS_pathway_analysis/jaccard_idx.R) |


### 8. Correlation and overlap with genetic susceptibility loci

| File                 | Link |
|----------------------|-------------|
| code/Overlap_with_AD_associated_genetics_loci/Overlap_with_AD_associated_genetics_loci.R       |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/Overlap_with_AD_associated_genetics_loci/Overlap_with_AD_associated_genetics_loci.R) |


### 9. Correlation of AD-associated CpGs and DMRs methylation levels in blood and brain samples	

| File                 | Link |
|----------------------|-------------|
| code/MethReg/Blood_MethReg_DMR_cpg.R        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood/blob/main/code/MethReg/Blood_MethReg_DMR_cpg.R) |

# For reproducible research

The following R packages are required: 

```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.13",ask = FALSE) # Install last version of Bioconductor

list.of.packages <- c("readr", "readxl", "plyr", "dplyr", "tidyr", "GenomicRanges", "SummarizedExperiment", "mygene", "ggpubr", "stats", "gt", "fgsea", "minfi", "bacon", "wateRmelon", "missMethyl", "DMRcate", "lumi", "RPMM","sm", "MethReg", "writexl", "readxl", "gridExtra", "doParallel", "ReMapEnrich", "EpiDISH", "GEOquery", "ggrepel",
"DT", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "RVenn", "GWASTools", "meta", "metap", "ExperimentHub", "lubridate")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr")
```

For ADNIMERGE, download it from https://ida.loni.usc.edu/: Merged ADNI 1/GO/2 Packages for R

```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```
The platform information are: 

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

# Acknowledgement
Data used in preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative
(ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design
and implementation of ADNI and/or provided data but did not participate in analysis or writing of this report. A complete listing of ADNI investigators can be found at:
http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

# References

1. Vasanthakumar, A. et al. Harnessing peripheral DNA methylation differences in the Alzheimer's Disease Neuroimaging Initiative (ADNI) to reveal novel biomarkers of disease. Clin Epigenetics 12, 84 (2020).

2. Ellis, K.A. et al. Enabling a multidisciplinary approach to the study of ageing and Alzheimer's disease: an update from the Australian Imaging Biomarkers and Lifestyle (AIBL) study. Int Rev Psychiatry 25, 699-710 (2013).
