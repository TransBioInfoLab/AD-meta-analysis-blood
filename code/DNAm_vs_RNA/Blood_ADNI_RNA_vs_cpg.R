library(dplyr)
library(SummarizedExperiment)

setwd("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples")
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
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0
)
cpgs.prioritized  <- cpgs.prioritized[[5]] %>% na.omit() %>% as.character
length(cpgs.prioritized) # 97


# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)
combp.cpgs <- combp_AD_vs_CN$Probes %>% sapply(FUN = function(x){stringr::str_split(x,";")}) %>% unlist

cpgs.all <- c(
  combp.cpgs,
  cpgs.prioritized,
  cpgs.ad.cn
) %>% unique

load("~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/datasets/Aux/ADNI_matched_rna_dnam_residuals.rda")

#-----------------------------------------------------------------------------
# get residuals we will only be using CN and Dementia; MCI will be removed
#-----------------------------------------------------------------------------

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

library(MethReg)

#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------
auxfunction <- function(row){
  
  rna.target <- residuals.matched.exp[which(rownames(residuals.matched.exp) == row$target), , drop = FALSE]
  met <- dnam[which(rownames(dnam) == as.character(row$regionID)), , drop = FALSE]
  
  data <- data.frame(
    "met.residual" = met %>% as.numeric(),
    "rna.target.residual" = rna.target %>% as.numeric(),
    "AD_Status" = factor(metadata.dnam$DX,levels = c("Dementia","CN"))
  )
  
  rlm <- MASS::rlm(
    rna.target.residual ~ met.residual + AD_Status, 
    data = data,
    psi = MASS::psi.bisquare,
    maxit = 100
  ) %>% summary %>% coef %>% data.frame
  
  degrees.freedom.value <- nrow(data) - 3
  rlm$pval <- 2 * (1 - pt( abs(rlm$t.value), df = degrees.freedom.value) )
  
  quant.pval <- rlm[-1,4,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")
  
  quant.estimate <- rlm[-1,1,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")
  return(cbind(quant.pval, quant.estimate))
}


add_cpgs_from_and_do_fdr <- function(results){
  idx.ad.vs.cn <- which(results$probeID %in% cpgs.ad.cn)
  idx.prioritized <- which(results$probeID %in% cpgs.prioritized)
  idx.combp <- which(results$probeID %in% combp.cpgs)
  
  results <- results[c(idx.ad.vs.cn,idx.prioritized,idx.combp),]
  results$cpgs_from <- c(
    rep("50 CpGs from AD vs CN meta-analysis", length(idx.ad.vs.cn)),
    rep("97 prioritized CpGs", length(idx.prioritized)),
    rep("CpGs in 9 DRMs from Combp", length(idx.combp))
  )
  
  
  results$RLM_met.residual_fdr <- c(
    p.adjust(results$RLM_met.residual_pvalue[which( results$cpgs_from %in% "50 CpGs from AD vs CN meta-analysis")],method = "fdr"),
    p.adjust(results$RLM_met.residual_pvalue[which( results$cpgs_from %in% "97 prioritized CpGs")],method = "fdr"),
    p.adjust(results$RLM_met.residual_pvalue[which( results$cpgs_from %in% "CpGs in 9 DRMs from Combp")],method = "fdr")
  )
  return(results)
}

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
dnam <- MethReg:::map_probes_to_regions(
  dnam = residuals.matched.met,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")
promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

# 12 cpgs in cpgs.ad.cn are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[cpgs.ad.cn]))

# 40 cpgs in cpgs.prioritized are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[cpgs.prioritized]))

nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(promoter.gene.dnam.pair)
results.promoter.analysis <- plyr::adply(promoter.gene.dnam.pair,.margins = 1,.fun = auxfunction)

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr

results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05]  %in%  make_names_from_granges(EPIC.hg19[cpgs.prioritized])
results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05]  %in% make_names_from_granges(EPIC.hg19[cpgs.ad.cn])
results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

# 38 cpgs in cpgs.ad.cn are promoter
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[cpgs.ad.cn]))

# 54 cpgs in cpgs.ad.cn are promoter
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[cpgs.prioritized]))


nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(distal.gene.dnam.pair)

results.distal.analysis <- plyr::adply(distal.gene.dnam.pair,.margins = 1,.fun = auxfunction)

results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr

results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05]  %in% make_names_from_granges(EPIC.hg19[cpgs.prioritized])
results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05]  %in% make_names_from_granges(EPIC.hg19[cpgs.ad.cn])
results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# window analysis
#-------------------------------------------------------------------------------
window.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "window",
  genome = "hg19",
  window.size = 500 * 10^3
)
window.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(window.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]

window.gene.dnam.pair <- window.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
results.window.analysis <- plyr::adply(window.gene.dnam.pair,.margins = 1,.fun = auxfunction)

results.window.analysis <- results.window.analysis %>% add_cpgs_from_and_do_fdr

results.window.analysis$regionID[results.window.analysis$RLM_met.residual_fdr < 0.05]  %in% make_names_from_granges(EPIC.hg19[cpgs.prioritized])
results.window.analysis$regionID[results.window.analysis$RLM_met.residual_fdr < 0.05]  %in% make_names_from_granges(EPIC.hg19[cpgs.ad.cn])
results.window.analysis[results.window.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------

path.RNA_vs_DNAm <- "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/"
writexl::write_xlsx(
  list(
    "Promoter" = results.promoter.analysis,
    "Distal_10_up_10_down" = results.distal.analysis,
    "Window_500kb" = results.window.analysis
  ),
  path = file.path(path.RNA_vs_DNAm,"Blood_cpg_Target_vs_DNAm.xlsx")
)
