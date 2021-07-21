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
library(dplyr)
library(SummarizedExperiment)
library(coMethDMR)
#-----------------------------------------------------------------------------
# Analysis: target gene ~ CpG  uisng ROSMAP data
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
library(readr)
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


combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 2 DMRs-Combp-AD_vs_CN_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

combp.cpgs <- GetCpGsInRegion(
  combp_AD_vs_CN$DMR,
  arrayType = "EPIC"
) 

cpgs.all <- c(
  combp.cpgs,
  cpgs.prioritized,
  cpgs.ad.cn
) %>% unique


#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP DNAm data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# Data not included
load("../coMethDMR_metaAnalysis/DNAm_RNA/data/matched_data.rda")
dim(matched.dnam)
dim(matched.exp)
matched.exp <- matched.exp[rowSums(matched.exp) > 0,]


# 1) remove confounding effects in DNAm data: 
resid_met <- coMethDMR:::GetResiduals(
  dnam = matched.dnam[rownames(matched.dnam) %in% cpgs.all,],
  betaToM = TRUE, #converts to Mvalues for fitting linear model 
  pheno_df = matched.phenotype,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","msex","age_death"), 
  nCores_int = 1,
  progressbar = TRUE  
)

resid_dnam <- MethReg:::map_probes_to_regions(
  dnam = resid_met,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# Get linked genes
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")

promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]

distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]


window.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "window",
  genome = "hg19",
  window.size = 500 * 10^3,
  rm.promoter.regions.from.distal.linking = FALSE
)
window.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(window.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]


#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP Gene expression data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
matched.exp.log2 <- log2(matched.exp + 1) # + 1 is required otherwise -INF
markers <-
  t(matched.exp.log2[c(
    "ENSG00000111674",
    "ENSG00000129226",
    "ENSG00000131095",
    "ENSG00000205927",
    "ENSG00000174059"
  ), ])
colnames(markers) <- c(
  "markers_ENO2",
  "markers_OLIG2",
  "markers_CD34",
  "markers_CD68",
  "markers_GFAP"
)

matched.exp.log2 <- matched.exp.log2[ rownames(matched.exp.log2) %in% 
                                        c(distal.gene.dnam.pair$target,
                                          promoter.gene.dnam.pair$target,
                                          window.gene.dnam.pair$target),]

matched.phenotype$rnaseq_id  <- map$rnaseq_id[match(matched.phenotype$Sample,map$mwas_id)]
residuals.matched.exp <- plyr::adply(
  .data = matched.exp.log2,
  .margins = 1, 
  function(row){
    val <- t(row)
    colnames(val) <- "val"
    dat <- cbind(
      val, 
      matched.phenotype,
      markers
    )
    dat$val <- as.numeric(dat$val)
    fitE <- lm(
      "val ~ age_death + msex + markers_ENO2 + markers_OLIG2 + markers_CD34 + markers_CD68 + markers_GFAP", 
      data = dat, 
      na.action = na.exclude
    )
    residuals(fitE)
  }, .progress = "time",
  .parallel = FALSE)
rownames(residuals.matched.exp) <- rownames(matched.exp.log2)

#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------
auxfunction <- function(row){
  
  rna.target <- residuals.matched.exp[which(rownames(residuals.matched.exp) == row$target), , drop = FALSE]
  met <- resid_dnam[which(rownames(resid_dnam) == as.character(row$regionID)), , drop = FALSE]
  
  data <- data.frame(
    "met.residual" = met %>% as.numeric(),
    "rna.target.residual" = rna.target %>% as.numeric(),
    Braak_stage = matched.phenotype$braaksc %>% as.numeric()
  )
  
  rlm <- MASS::rlm(
    rna.target.residual ~ met.residual + Braak_stage, 
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

promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
results.promoter.analysis <- plyr::adply(promoter.gene.dnam.pair,.margins = 1,.fun = auxfunction)

# Where did the cpg come from ? 

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr

results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05]  %in%  make_names_from_granges(EPIC.hg19[cpgs.prioritized])
results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
nrow(distal.gene.dnam.pair) #  1360
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
nrow(distal.gene.dnam.pair) # 1044

results.distal.analysis <- plyr::adply(distal.gene.dnam.pair,.margins = 1,.fun = auxfunction)

results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr

results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05]  %in%  make_names_from_granges(EPIC.hg19[cpgs.prioritized])
results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Save resuls
#-------------------------------------------------------------------------------
path.RNA_vs_DNAm <- "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm/"
writexl::write_xlsx(
  list(
    "Promoter" = results.promoter.analysis,
    "Distal_10_up_10_down" = results.distal.analysis
  ),
  path = file.path(path.RNA_vs_DNAm,"Brain_cpg_Target_vs_DNAm.xlsx")
)


