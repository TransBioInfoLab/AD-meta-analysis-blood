#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
# Article
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
#  Integrative meta-analysis of epigenome-wide association studies identifies 
#  genomic and epigenomics differences in the brain and the blood in Alzheimer’s disease
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
# Authors
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
# - Tiago C Silva
# - Juan I. Young
# - Lanyu Zhang
# - Lissette Gomez
# - Michael A. Schmidt
# - Achintya Varma
# - Xi Chen
# - Eden R. Martin
# - Lily Wang
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
library(dplyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
# Before you start
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-
#  London data processing was performed in our previous study
#  Code to acquire the data is available at:
# https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/single_cohort_analysis/London.Rmd
# https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/single_cohort_analysis/London_blood.Rmd
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-

#-------------------------------------------------------------------------------
# Data retrieval
#-------------------------------------------------------------------------------
main.dir <- "../coMethDMR_metaAnalysis/code_validation/Meta_analysis_code/"
data.brain.beta <- file.path(main.dir,"DATASETS/LONDON/step5_pca_filtering/")
data.brain.pheno <- file.path(main.dir,"DATASETS/LONDON/step6_neuron_comp/")
data.blood.beta <- file.path(main.dir,"DATASETS/LONDON_blood/step5_pca_filtering/")
data.blood.pheno <- file.path(main.dir,"DATASETS/LONDON_blood/step6_neuron_comp/")
data.dmr <- file.path(main.dir,"meta_analysis_region_results/step4_dmr_vs_cpgs/")
data.cpg <- file.path(main.dir,"meta_analysis_single_cpg_results/")
data.final <- file.path(main.dir,"London_blood_brain_correlation_results/")
data.final.beta <- file.path(main.dir,"London_blood_brain_correlation_results/using_betas/")
data.final.resid <- file.path(main.dir,"London_blood_brain_correlation_results/using_residuals/")
data.BECon <- file.path(main.dir,"DATASETS/LONDON_blood/step10_blood_brain_correlation/")

brain_beta <- readRDS(
  paste0(data.brain.beta, "London_PFC_QNBMIQ_PCfiltered_withStageExclude.RDS")
) # 450793    111

brain_pheno <- readRDS(
  paste0(data.brain.pheno, "pheno107_PFC_withNeuronProp_withStageExclude_df.RDS")
)

blood_beta <- readRDS(
  paste0(data.blood.beta, "London_QNBMIQ_PCfiltered_withStatusExclude.RDS")
) #  450793     77

blood_pheno <- readRDS(
  paste0(data.blood.pheno, "pheno_BLOOD_withBloodProp_withStatusExclude_df.rds")
)

# Limit samples in both datasets

### Renames variables
colnames(brain_pheno)[c(1, 3:ncol(brain_pheno))] <- paste0(
  "brain_", colnames(brain_pheno)[c(1, 3:ncol(brain_pheno))]
)
colnames(blood_pheno)[c(1, 3:ncol(blood_pheno))] <- paste0(
  "blood_", colnames(blood_pheno)[c(1, 3:ncol(blood_pheno))]
)

### Merge datasets
pheno_final <- merge(
  brain_pheno, blood_pheno,
  by = "subject.id"
) #dim: 69 23

### Limit beta matrices to samples in pheno_final
brain_beta_final <- brain_beta[, pheno_final$brain_sample]
blood_beta_final <- blood_beta[, pheno_final$blood_sample]

# Calculate blood and brain correlation (without taking residuals)

### Call in datasets with sig DMRs and CpGs
main_dmrs <- read.csv(
  paste0(data.dmr, "meta_analysis_sig_no_crossHyb_smoking_ov_comb_p_with_sig_single_cpgs.csv")
)

main_cpgs <- read.csv(
  paste0(data.cpg, "meta_analysis_single_cpg_sig_no_crossHyb_smoking_df.csv")
)

## for sig. regions

### Get probes from regions
probes.cluster.all <- coMethDMR:::getPredefinedCluster(
  arrayType = "450k",
  clusterType = "regions"
)

idx <- gsub("450k_Gene_3_200.|450k_InterGene_3_200.","",names(probes.cluster.all)) %in% main_dmrs$inputRegion
main_dmrs_cpgs <- probes.cluster.all[idx] %>% unlist %>% as.character() %>% unique

### Limit blood_beta and brain_beta to the probes above
brain_beta_regions <- brain_beta_final
blood_beta_regions <- blood_beta_final

identical(dim(brain_beta_regions), dim(blood_beta_regions))
identical(row.names(brain_beta_regions), row.names(blood_beta_regions))

doParallel::registerDoParallel(cores = 10)
blood_brain_cor <- plyr::adply(
  seq_len(nrow(brain_beta_regions)),
  .margins = 1,
  .fun =  function(row){
    
    spearman_cor <- cor.test(
      brain_beta_regions[row,],
      blood_beta_regions[row,],
      method = "spearman"
    )
    
    data.frame(
      cpg = row.names(brain_beta_regions)[row],
      spearman_cor = spearman_cor$estimate,
      pVal = spearman_cor$p.value,
      stringsAsFactors = FALSE
    )
    
  },.id = NULL, .progress = "time",.parallel = TRUE)

blood_brain_cor$fdr <- p.adjust(blood_brain_cor$pVal, method = "fdr")
blood_brain_cor

write.csv(
  blood_brain_cor,
  paste0("analysis_results/London_blood/London_blood_brain_beta_correlation.csv"),
  row.names = FALSE
)

# Calculate blood and brain correlation (with residuals)

## Take residuals

### for brain beta matrix

### Compute M values
mvalue_mat <- log2( brain_beta_final /(1 - brain_beta_final))

### Reorder samples based on pheno_df
mvalue_mat <- mvalue_mat[, pheno_final$brain_sample]

identical(colnames(mvalue_mat), pheno_final$brain_sample)

### Take residuals
lmF <- function(mval){
  fitE <- lm(
    as.numeric(mval) ~ brain_age.brain + brain_sex + brain_prop.neuron + as.character(brain_slide), #add batch if rosmap
    data = pheno_final,
    na.action = na.exclude
  )
  residuals (fitE)
}

library(doParallel)
registerDoParallel(detectCores()/2)
resid <- plyr::adply(mvalue_mat,1,.fun = lmF,.progress = "time",.parallel = FALSE)
rownames(resid) <- resid[,1]
resid[,1] <- NULL
colnames(resid) <- colnames(mvalue_mat)

saveRDS(
  resid,
  paste0(data.final.resid, "London_PFC_QNBMIQ_PCfiltered_mvalResiduals.RDS")
)

### for blood beta matrix

### Compute M values
mvalue_mat <- log2(blood_beta_final / (1 - blood_beta_final))

### Reorder samples based on pheno_df
mvalue_mat <- mvalue_mat[, pheno_final$blood_sample]

identical(colnames(mvalue_mat),  pheno_final$blood_sample)

lmF <- function(mval){
  fitE <- lm(
    as.numeric(mval) ~ blood_age.blood + blood_sex + blood_slide +
      blood_B + blood_NK + blood_CD4T + blood_CD8T + blood_Mono + blood_Neutro + blood_Eosino,
    data = pheno_final,
    na.action = na.exclude
  )
  residuals (fitE)
}

resid <- plyr::adply(mvalue_mat,1,.fun = lmF,.progress = "time",.parallel = TRUE)
rownames(resid) <- resid[,1]
resid[,1] <- NULL
colnames(resid) <- colnames(mvalue_mat)

saveRDS(
  resid,
  paste0(data.final.resid, "LONDON_blood_QNBMIQ_PCfiltered_mvalResiduals.RDS")
)

## Call in datasets

### Call in brain and blood residual matrices
brain_beta_final <- as.matrix(
  readRDS(
    paste0(data.final.resid, "London_PFC_QNBMIQ_PCfiltered_mvalResiduals.RDS")
  )
)
blood_beta_final <- as.matrix(
  readRDS(
    paste0(data.final.resid, "LONDON_blood_QNBMIQ_PCfiltered_mvalResiduals.RDS")
  )
)

### Limit blood_beta and brain_beta to the probes above
brain_beta_regions <- brain_beta_final
blood_beta_regions <- blood_beta_final

identical(dim(brain_beta_regions), dim(blood_beta_regions))
identical(row.names(brain_beta_regions), row.names(blood_beta_regions))

blood_brain_cor <- plyr::adply(
  .data = seq_len(nrow(brain_beta_regions)),
  .margins = 1,
  .fun = function(row){
    spearman_cor <- cor.test(
      brain_beta_regions[row,],
      blood_beta_regions[row,],
      method = "spearman"
    )
    
    data.frame(
      cpg = row.names(brain_beta_regions)[row],
      spearman_cor = spearman_cor$estimate,
      pVal = spearman_cor$p.value,
      stringsAsFactors = FALSE
    )
  },.progress = "time",.parallel = TRUE,.id = NULL)

blood_brain_cor$fdr <- p.adjust(blood_brain_cor$pVal, method = "fdr")

write.csv(
  blood_brain_cor,
  paste0("analysis_results/London_blood/London_blood_brain_residuals_correlation.csv"),
  row.names = FALSE
)


# Merge final results

### Call in datasets
beta <- read.csv(
  "analysis_results/London_blood/London_blood_brain_beta_correlation.csv"
)
resid <- read.csv(
  "analysis_results/London_blood/London_blood_brain_residuals_correlation.csv"
)

### Rename variables
colnames(beta)[2:4] <- paste0("beta_", colnames(beta)[2:4])
colnames(resid)[2:4] <- paste0("residual_", colnames(resid)[2:4])

### Merge datasets
cor <- merge(
  beta, 
  resid,
  by = "cpg"
)

# Merge results with results from BECon

### Call in BECon results
becon.files <- dir("analysis_results/London_blood/",pattern = "becon_part",full.names = TRUE)
becon <- plyr::adply(becon.files,.margins = 1,.fun = function(f){read.csv(f)})


### Select and rename variables
becon <- becon[
  ,c("CpG.ID", "Cor.Blood.BA7", "Cor.Blood..BA10", "Cor.Blood..BA20", "Mean.Cor.All.Brain")
]
colnames(becon) <- c(
  "cpg", "BECon_cor_BA7", "BECon_cor_BA10", "BECon_cor_BA20", "BECon_cor_mean"
)


### Merge BECon results with our results
final <- merge(
  cor, 
  becon,
  by = "cpg",
  all.x = TRUE
)

write.csv(
  final,
  paste0("analysis_results/London_blood/London_blood_brain_correlation_cpgs.csv"),
  row.names = FALSE
)


# Merge results 


library(dplyr)
tab1 <- read.csv( paste0("analysis_results/London_blood/London_blood_brain_correlation_cpgs.csv"))

tab2 <- "analysis_results/RNA_vs_DNAm//Blood_and_brain_merged_all_blood_cpgs.xlsx" %>% 
  readxl::read_xlsx( sheet = 2) %>%
  dplyr::rename(cpg = probeID)

# File in code/others
# Does not have all probes due to hm450 and EPIC difference
tab3 <- readxl::read_xlsx( paste0("tables//all-CpGs_and_DMRs-crossTissue_brain_blood_with_NAs.xlsx"),sheet = 1)

merged <- tab2 %>% dplyr::full_join(tab3) %>% dplyr::full_join(tab1)
merged <- merged[!is.na(merged$cpgs_from),]
writexl::write_xlsx(merged,"tables//merged_table_RNA_DNAm_cross_tissue_meta_analysis_london_brain_blood_cor.xlsx")


