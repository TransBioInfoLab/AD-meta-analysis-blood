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
library(SummarizedExperiment)
library(coMethDMR)
library(readr)
library(readxl)

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted-V2.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

cpgs.prioritized <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 3 prioritized-CpGs-crossTissue_brain_blood.xlsx",skip = 3
)
cpgs.prioritized  <- cpgs.prioritized$CpG %>% na.omit() %>% as.character
length(cpgs.prioritized)

cpgs.all <- c(
  cpgs.prioritized,
  cpgs.ad.cn
) %>% unique

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 2 DMRs-Combp-AD_vs_CN_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

combp.probes <- GetCpGsInRegion(
  combp_AD_vs_CN$DMR,
  arrayType = "EPIC"
) 

#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 3 Top 10 prioritized-CpGs_and_DMRs-crossTissue_brain_blood-V2.xlsx",sheet = 1, skip = 17)
prioritized.dmrs <- prioritized.dmrs$DMRs %>% na.omit %>% as.character
length(prioritized.dmrs) # 10

prioritized.dmrs.probes <- GetCpGsInRegion(
  prioritized.dmrs,
  arrayType = "EPIC"
) # 19 DMRs

df <- rbind(
  data.frame(
    "cpg" = cpgs.prioritized,
    "cpg is from" = "prioritized cpg"
  ),
  data.frame(
    "cpg" = cpgs.ad.cn,
    "cpg is from" = "ad vs. cn"
  ),
  data.frame(
    "cpg" = combp.probes %>% strsplit(split = ";") %>% unlist %>% unique,
    "cpg is from" = "comb-p"
  ),
  data.frame(
    "cpg" = prioritized.dmrs.probes %>% strsplit(split = ";") %>% unlist %>% unique,
    "cpg is from" = "prioritized dmr"
  )
)
blood.resuls <- readr::read_csv("analysis_results/London_blood/London_blood_brain_correlation_cpgs.csv")
load("datasets/Aux/great_EPIC_array_annotation.rda")

ret <- dplyr::left_join(df, blood.resuls) %>% dplyr::left_join(great[,c("cpg","GREAT_annotation")])

g1 <- c("prioritized cpg", "prioritized dmr")
g2 <- c("ad vs. cn", "comb-p")         

ret$beta_fdr[which(ret$cpg.is.from %in% g1)] <- p.adjust(  ret$beta_pVal[which(ret$cpg.is.from %in% g1)],method = "fdr")
ret$residual_fdr[which(ret$cpg.is.from %in% g1)] <- p.adjust(  ret$residual_pVal[which(ret$cpg.is.from %in% g1)],method = "fdr")

ret$beta_fdr[which(ret$cpg.is.from %in% g2)] <- p.adjust(  ret$beta_pVal[which(ret$cpg.is.from %in% g2)],method = "fdr")
ret$residual_fdr[which(ret$cpg.is.from %in% g2)] <- p.adjust(  ret$residual_pVal[which(ret$cpg.is.from %in% g2)],method = "fdr")

# Update FDR

writexl::write_xlsx(
  list(
    "Intersection of results" = ret[!is.na(ret$beta_spearman_cor),],
    "All input" = ret
  ), 
  path = "tables/Correlation_in_brain_and_blood_samples.xlsx"
)

