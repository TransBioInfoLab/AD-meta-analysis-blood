library(dplyr)
library(SummarizedExperiment)

path.mathReg <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/analysis_results/RNA_vs_DNAm"
dir.create(path.mathReg,recursive = TRUE,showWarnings = FALSE)

#-----------------------------------------------------------------------------
# Analysis: target gene ~ CpG  uisng ROSMAP data
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
#  in fdr significant DMRs prioritized in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",sheet = 2)
prioritized.dmrs <- prioritized.dmrs[[4]] %>% na.omit %>% as.character
length(prioritized.dmrs)

devtools::load_all("~/Documents/packages/coMethDMR/")
prioritized.dmrs.probes <- GetCpGsInAllRegion(
  prioritized.dmrs,
  arrayType = "EPIC"
) # 19 DMRs
prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(prioritized.dmrs.probes,FUN = function(x) paste(x,collapse = ";"))

regions <- rbind(prioritized.dmrs[,c("DMR","Probes")])

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP DNAm data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
load("~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/DNAm_RNA/data/matched_data.rda")
dim(matched.dnam)
dim(matched.exp)
matched.exp <- matched.exp[rowSums(matched.exp) > 0,]


matched.dnam.median <- plyr::ldply(
  regions$Probes,
  .fun = function(x){
    cpgs <- stringr::str_split(x,";") %>% unlist; 
    colMedians(matched.dnam[rownames(matched.dnam) %in% cpgs,],na.rm = TRUE)
  },.id = NULL
)
colnames(matched.dnam.median) <- colnames(matched.dnam)
rownames(matched.dnam.median) <- regions$DMR
matched.dnam.median <- na.omit(matched.dnam.median)


# 1) remove confounding effects in DNAm data: 
resid_dnam <- GetResiduals(
  dnam = matched.dnam.median,
  betaToM = TRUE, #converts to Mvalues for fitting linear model 
  pheno_df = matched.phenotype,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","msex","age_death"), 
  nCores_int = 1,
  progressbar = TRUE  
)

all(colnames(resid_dnam) == matched.phenotype$Sample[match(gsub("_[0-9]$","",colnames(matched.exp.log2)),matched.phenotype$rnaseq_id)])

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# Get linked genes
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=

promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)



window.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(resid_dnam) %>% MethReg::make_granges_from_names(),
  method = "window",
  genome = "hg19",
  window.size = 500 * 10^3
)

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

matched.exp.log2 <- matched.exp.log2[ rownames(matched.exp.log2) %in% c(window.gene.dnam.pair$target, distal.gene.dnam.pair$target,promoter.gene.dnam.pair$target),]

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

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------

promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
results.promoter.analysis <- plyr::adply(promoter.gene.dnam.pair,.margins = 1,.fun = auxfunction)

# Where did the DMRs come from ? 
results.promoter.analysis$regions_from <- "Regions_prioritized"
results.promoter.analysis$RLM_met.residual_fdr <- p.adjust(results.promoter.analysis$RLM_met.residual_pvalue,method = "fdr")

results.promoter.analysis$regionID[results.promoter.analysis$RLM_met.residual_fdr < 0.05] %in% prioritized.dmrs$DMR
results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
nrow(distal.gene.dnam.pair) #  1360
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
nrow(distal.gene.dnam.pair) # 1044

results.distal.analysis <- plyr::adply(distal.gene.dnam.pair,.margins = 1,.fun = auxfunction)

# Where did the region come from ? 
results.distal.analysis$regions_from <- "Regions_prioritized"
results.distal.analysis$RLM_met.residual_fdr <- p.adjust(results.distal.analysis$RLM_met.residual_pvalue,method = "fdr")

results.distal.analysis$regionID[results.distal.analysis$RLM_met.residual_fdr < 0.05] %in% prioritized.dmrs$DMR
results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# window analysis
#-------------------------------------------------------------------------------
window.gene.dnam.pair <- window.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
results.window.analysis <- plyr::adply(window.gene.dnam.pair,.margins = 1,.fun = auxfunction)

# Where did the region come from ?
results.window.analysis$regions_from <- "Regions_prioritized"
results.window.analysis$RLM_met.residual_fdr <- p.adjust(results.window.analysis$RLM_met.residual_pvalue,method = "fdr")

results.window.analysis$regionID[results.window.analysis$RLM_met.residual_fdr < 0.05] %in% prioritized.dmrs$DMR
results.window.analysis[results.window.analysis$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------
writexl::write_xlsx(
  list(
    "Promoter" = results.promoter.analysis,
    "Distal_10_up_10_down" = results.distal.analysis,
    "Window_500kb" = results.window.analysis
  ),
  path = file.path(path.mathReg,"Brain_DMR_Target_vs_DNAm.xlsx")
)