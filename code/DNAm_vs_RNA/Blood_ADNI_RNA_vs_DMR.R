#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Integrative meta-analysis of epigenome-wide association studies
# identifies genomic and
# epigenomics differences in the brain and the blood in Alzheimer’s disease
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Authors: 
# - Tiago C. silva
# - Lily Wang
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Date: 12 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#-----------------------------------------------------------------------------
# MethReg analysis without TF
# target gene ~ CpG 
#-----------------------------------------------------------------------------
library(coMethDMR)
path.mathReg <- "AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm"
dir.create(path.mathReg,recursive = TRUE,showWarnings = FALSE)

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - CpGs within significant DMRs (with length > 3cpgs) identified by combp

combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 2 DMRs-Combp-AD_vs_CN_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

combp_AD_vs_CN.probes <- plyr::alply(combp_AD_vs_CN$DMR,.margins = 1,.fun = function(dmr) {
  GetCpGsInRegion(
    dmr,
    arrayType = "EPIC"
  )}) # 19 DMRs
combp_AD_vs_CN$Probes <- sapply(
  combp_AD_vs_CN.probes,
  FUN = function(x) paste(x,collapse = ";")
)


#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 3 Top 10 prioritized-CpGs_and_DMRs-crossTissue_brain_blood-V2.xlsx",sheet = 1, skip = 17)
prioritized.dmrs <- prioritized.dmrs$DMRs %>% na.omit %>% as.character
length(prioritized.dmrs) # 10


prioritized.dmrs.probes <- plyr::alply(prioritized.dmrs,.margins = 1,.fun = function(dmr) {
  GetCpGsInRegion(
  dmr,
  arrayType = "EPIC"
)}) # 19 DMRs
prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(
  prioritized.dmrs.probes,
  FUN = function(x) paste(x,collapse = ";")
)

regions <- rbind(combp_AD_vs_CN[,c("DMR","Probes")],prioritized.dmrs[,c("DMR","Probes")])

load("datasets/Aux/ADNI_matched_rna_dnam_residuals_DMR.rda")

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------
auxfunction <- function(row){
  
  rna.target <- residuals.matched.exp[which(rownames(residuals.matched.exp) == row$target), , drop = FALSE]
  met <- residuals.matched.met[which(rownames(residuals.matched.met) == as.character(row$regionID)), , drop = FALSE]
  
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

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------

promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(residuals.matched.met) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

table(promoter.gene.dnam.pair$regionID %>% unique   %in% prioritized.dmrs$DMR)
table(promoter.gene.dnam.pair$regionID  %>% unique %in% combp_AD_vs_CN$DMR)

promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
ret <- plyr::adply(promoter.gene.dnam.pair,.margins = 1,.fun = auxfunction)
results.promoter.analysis <- ret

results.promoter.analysis$regions_from <- "Regions_prioritized"
results.promoter.analysis$regions_from[results.promoter.analysis$regionID %in% combp_AD_vs_CN$DMR] <- "comb-p analysis"
results.promoter.analysis$regions_from[results.promoter.analysis$regionID %in% base::intersect(prioritized.dmrs$DMR,combp_AD_vs_CN$DMR)] <- "comb-p analysis/regions_prioritized"

results.promoter.analysis$RLM_met.residual_fdr <- NA
results.promoter.analysis$RLM_met.residual_fdr[results.promoter.analysis$regionID %in% combp_AD_vs_CN$DMR] <- p.adjust(results.promoter.analysis$RLM_met.residual_pvalue[results.promoter.analysis$regionID %in% combp_AD_vs_CN$DMR],method = "fdr")
results.promoter.analysis$RLM_met.residual_fdr[results.promoter.analysis$regionID %in% prioritized.dmrs$DMR] <- p.adjust(results.promoter.analysis$RLM_met.residual_pvalue[results.promoter.analysis$regionID %in% prioritized.dmrs$DMR],method = "fdr")

results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05]  %in% prioritized.dmrs$DMR
results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05]  %in% combp_AD_vs_CN$DMR
results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(residuals.matched.met) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  genome = "hg19"
)

table(distal.gene.dnam.pair$regionID %>% unique %in% prioritized.dmrs$DMR)
table(distal.gene.dnam.pair$regionID %>% unique %in% combp_AD_vs_CN$DMR)


distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
ret <- plyr::adply(distal.gene.dnam.pair,.margins = 1,.fun = auxfunction)
results.distal.analysis <- ret

results.distal.analysis$regions_from <- "Regions_prioritized"
results.distal.analysis$regions_from[results.distal.analysis$regionID %in% combp_AD_vs_CN$DMR] <- "comb-p analysis"
results.distal.analysis$regions_from[results.distal.analysis$regionID %in% base::intersect(prioritized.dmrs$DMR,combp_AD_vs_CN$DMR)] <- "comb-p analysis/regions_prioritized"
results.distal.analysis$RLM_met.residual_fdr <- NA
results.distal.analysis$RLM_met.residual_fdr[results.distal.analysis$regionID %in% combp_AD_vs_CN$DMR] <- p.adjust(results.distal.analysis$RLM_met.residual_pvalue[results.distal.analysis$regionID %in% combp_AD_vs_CN$DMR],method = "fdr")
results.distal.analysis$RLM_met.residual_fdr[results.distal.analysis$regionID %in% prioritized.dmrs$DMR] <- p.adjust(results.distal.analysis$RLM_met.residual_pvalue[results.distal.analysis$regionID %in% prioritized.dmrs$DMR],method = "fdr")

results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05]  %in% prioritized.dmrs$DMR
results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05]  %in% combp_AD_vs_CN$DMR
results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    "Promoter" = results.promoter.analysis,
    "Distal_10_up_10_down" = results.distal.analysis
  ),
  path = file.path(path.mathReg,"Blood_DMR_Target_vs_DNAm.xlsx")
)
