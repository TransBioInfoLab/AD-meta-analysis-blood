library(dplyr)
library(SummarizedExperiment)
#-----------------------------------------------------------------------------
# MethReg analysis without TF
# target gene ~ CpG 
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
library(readr)
AD_vs_CN <- readxl::read_xlsx(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

cpgs.prioritized <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0
)
cpgs.prioritized  <- cpgs.prioritized[[5]] %>% na.omit() %>% as.character
length(cpgs.prioritized)

cpgs.all <- c(
  cpgs.prioritized,
  cpgs.ad.cn
) %>% unique


#-----------------------------------------------------------------------------
# target gene ~ CpG 
#-----------------------------------------------------------------------------
devtools::load_all("~/Documents/packages/coMethDMR/")

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",sheet = 2)
prioritized.dmrs <- prioritized.dmrs[[4]] %>% na.omit %>% as.character
length(prioritized.dmrs) # 10

prioritized.dmrs.probes <- GetCpGsInAllRegion(
  prioritized.dmrs,
  arrayType = "EPIC"
) # 19 DMRs
prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(
  prioritized.dmrs.probes,
  FUN = function(x) paste(x,collapse = ";")
)

regions <- rbind(combp_AD_vs_CN[,c("DMR","Probes")],prioritized.dmrs[,c("DMR","Probes")])
regions.cpgs <- strsplit(regions$Probes,split = ";") %>% unlist %>% unique

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
    "cpg" = combp_AD_vs_CN$Probes %>% strsplit(split = ";") %>% unlist %>% unique,
    "cpg is from" = "comb-p"
  ),
  data.frame(
    "cpg" = prioritized.dmrs$Probes %>% strsplit(split = ";") %>% unlist %>% unique,
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

