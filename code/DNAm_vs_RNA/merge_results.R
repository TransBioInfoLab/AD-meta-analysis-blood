brain.dmr.promoter <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_DMR_Target_vs_DNAm.xlsx", sheet = 1)
brain.dmr.promoter$analysis <- "promoter"
brain.dmr.distal <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_DMR_Target_vs_DNAm.xlsx",sheet = 2)
brain.dmr.distal$analysis <- "distal"
brain.dmr.window <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_DMR_Target_vs_DNAm.xlsx",sheet = 3)
brain.dmr.window$analysis <- "window"
brain.dmr <- rbind(brain.dmr.promoter, brain.dmr.distal)
colnames(brain.dmr)[grep("met.residual",colnames(brain.dmr))] <- paste0("Brain_",gsub("RLM_met.residual_","",colnames(brain.dmr)[grep("met.residual",colnames(brain.dmr))]))
brain.dmr <- brain.dmr %>% dplyr::select(-contains("Braak"))
brain.dmr.sig <- brain.dmr %>% dplyr::filter(Brain_fdr < 0.05)

brain.cpg.promoter <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_cpg_Target_vs_DNAm.xlsx", sheet = 1)
brain.cpg.promoter$analysis <- "promoter"
brain.cpg.distal <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_cpg_Target_vs_DNAm.xlsx",sheet = 2)
brain.cpg.distal$analysis <- "distal"
brain.cpg.window <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Brain_cpg_Target_vs_DNAm.xlsx",sheet = 3)
brain.cpg.window$analysis <- "distal"
brain.cpg <- rbind(brain.cpg.promoter, brain.cpg.distal)
colnames(brain.cpg)[grep("met.residual",colnames(brain.cpg))] <- paste0("Brain_",gsub("RLM_met.residual_","",colnames(brain.cpg)[grep("met.residual",colnames(brain.cpg))]))
brain.cpg <- brain.cpg %>% dplyr::select(-contains("Braak"))
brain.cpg.sig <- brain.cpg %>% dplyr::filter(Brain_fdr < 0.05)

blood.dmr.promoter <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_DMR_Target_vs_DNAm.xlsx", sheet = 1)
blood.dmr.promoter$analysis <- "promoter"
blood.dmr.distal <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_DMR_Target_vs_DNAm.xlsx",sheet = 2)
blood.dmr.distal$analysis <- "distal"
blood.dmr.window <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_DMR_Target_vs_DNAm.xlsx",sheet = 3)
blood.dmr.window$analysis <- "window"
blood.dmr <- rbind(blood.dmr.promoter, blood.dmr.distal)
colnames(blood.dmr)[grep("met.residual",colnames(blood.dmr))] <- paste0("Blood_",gsub("RLM_met.residual_","",colnames(blood.dmr)[grep("met.residual",colnames(blood.dmr))]))
blood.dmr <- blood.dmr %>% dplyr::select(-contains("AD"))
blood.dmr.sig <- blood.dmr  %>% dplyr::filter(Blood_fdr < 0.05)

blood.cpg.promoter <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm//Blood_cpg_Target_vs_DNAm.xlsx", sheet = 1)
blood.cpg.promoter$analysis <- "promoter"
blood.cpg.distal <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_cpg_Target_vs_DNAm.xlsx",sheet = 2)
blood.cpg.distal$analysis <- "distal"
blood.cpg.window <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_cpg_Target_vs_DNAm.xlsx",sheet = 3)
blood.cpg.window$analysis <- "window"
blood.cpg <- rbind(blood.cpg.promoter, blood.cpg.distal)
colnames(blood.cpg)[grep("met.residual",colnames(blood.cpg))] <- paste0("Blood_",gsub("RLM_met.residual_","",colnames(blood.cpg)[grep("met.residual",colnames(blood.cpg))]))
blood.cpg <- blood.cpg %>% dplyr::select(-contains("AD"))
blood.cpg.sig <- blood.cpg  %>% dplyr::filter(Blood_fdr < 0.05)

cpgs <- inner_join(blood.cpg,brain.cpg) %>% dplyr::relocate(c("analysis","probeID","cpgs_from",1))
dmrs <- inner_join(blood.dmr,brain.dmr) %>% dplyr::relocate(c("analysis","regions_from",1))
writexl::write_xlsx(
  list(
    "CpGs (FDR < 0.05)" = cpgs[cpgs$Brain_fdr < 0.05 | cpgs$Blood_fdr < 0.05,],
    "CpGs" = cpgs,
    "DMR (FDR < 0.05)" = dmrs[dmrs$Brain_fdr < 0.05 | dmrs$Blood_fdr < 0.05,], 
    "DMR" = dmrs
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_and_brain_merged.xlsx"
)


cpgs <- left_join(blood.cpg,brain.cpg) %>% dplyr::relocate(c("analysis","probeID","cpgs_from",1))
dmrs <- left_join(blood.dmr,brain.dmr) %>% dplyr::relocate(c("analysis","regions_from",1))
writexl::write_xlsx(
  list(
    "CpGs (FDR < 0.05)" = cpgs[which(cpgs$Brain_fdr < 0.05 | cpgs$Blood_fdr < 0.05),],
    "CpGs" = cpgs,
    "DMR (FDR < 0.05)" = dmrs[which(dmrs$Brain_fdr < 0.05 | dmrs$Blood_fdr < 0.05),], 
    "DMR" = dmrs
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/Blood_and_brain_merged_all_blood_cpgs.xlsx"
)
