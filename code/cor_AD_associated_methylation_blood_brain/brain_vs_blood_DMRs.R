dmr.brain <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/DRAFT_REVISION_NatComm_10-12-2020/Supplementary Data 1-24.xlsx",2,skip = 3)

library(dplyr)
library(SummarizedExperiment)

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

# no overlap between brain meta-analysis DMRs and blood DMRs
subsetByOverlaps(
  combp_AD_vs_CN$DMR %>% make_granges_from_names(), 
  dmr.brain$DMR %>% make_granges_from_names()
)
