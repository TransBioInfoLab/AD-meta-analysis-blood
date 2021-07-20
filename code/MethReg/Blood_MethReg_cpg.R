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
library(readr)
library(MethReg)
library(readxl)
library(ReMapEnrich)
library(writexl)

#-----------------------------------------------------------------------------
# MethReg analysis
# target gene ~ TF_activity (dorothea) + CpG + CpG * TF
#-----------------------------------------------------------------------------
path.mathReg <- "analysis_results/methReg/Blood"
path.mathReg.plot <- file.path(path.mathReg, "plots/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

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

# need to add TF activity
load("datasets/Aux/ADNI_matched_rna_dnam_residuals.rda")

#-------------------------------------------------------------------------------
# Analysis
#-------------------------------------------------------------------------------
# Get triplets using remap
dir.base <- "."
dir.data.aux <- file.path(dir.base,"datasets/Aux/") 
remapCatalog2018hg19 <- downloadRemapCatalog(dir.data.aux, assembly = "hg19")
remapCatalog <- bedToGranges(remapCatalog2018hg19)

#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------

update_met_IQR <- function(results, dnam){
  iqr <- calculate_IQR(dnam)
  results$met.IQR <- iqr$IQR[match(results$regionID, iqr$ID)]
  results <- results %>% dplyr::relocate(dplyr::contains("IQR"), .after = last_col())
  return(results)
}


add_annot_cpgs <- function(results){
  
  results <- cbind(
    results,
    AD_vs_CN[
      match(results$probeID,AD_vs_CN$cpg),
      c("Islands.UCSC.Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Group","GREAT_annotation","E073_15_coreMarks_segments_state","sig.in.brain")
    ]
  )
  results
}


add_percent_zero_q1_q4 <- function(results, dnam, exp){
  
  aux <- plyr::adply(
    unique(results[,c("probeID","target")]),
    .margins = 1,
    .fun = function(row) {
      rna.target <- exp[rownames(exp) == row$target, , drop = FALSE]
      met <- dnam[rownames(dnam) == as.character(row$probeID), ]
      data <- data.frame(
        rna.target = rna.target %>% as.numeric,
        met = met %>% as.numeric
      )
      quant.met <-  quantile(data$met,na.rm = TRUE)
      low.cutoff <- quant.met[2]
      upper.cutoff <- quant.met[4]
      
      data.high.low <- data %>% filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
      data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)
      pct.zeros.in.quant.samples <- sum(
        data.high.low$rna.target == 0,
        na.rm = TRUE) / nrow(data.high.low)
      
      data.frame("% of 0 target genes (Q1 and Q4)" = paste0(round(pct.zeros.in.quant.samples * 100,digits = 2)," %"))
    }
  )
  results$`% of 0 residual target genes (Q1 and Q4)` <- results$`% of 0 target genes (Q1 and Q4)`
  results$`% of 0 target genes (Q1 and Q4)` <- aux$X..of.0.target.genes..Q1.and.Q4.[
    match(paste0(results$probeID,results$target),paste0(aux$probeID,aux$target))
  ]
  return(results)
}

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")

triplet.promoter.ewas <- create_triplet_distance_based(
  region = EPIC.hg19[rownames(residuals.matched.met),],
  genome = "hg19",
  target.method =  "genes.promoter.overlap",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
triplet.promoter.ewas <- triplet.promoter.ewas[!is.na(triplet.promoter.ewas$TF),]
nrow(triplet.promoter.ewas) # 4642 triplets

triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es))
triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))

nrow(triplet.promoter.ewas) # 1995 triplets
triplet.promoter.ewas$probeID <- names(EPIC.hg19)[match(triplet.promoter.ewas$regionID,make_names_from_granges(EPIC.hg19))]


triplet.promoter.ewas$TF %>% unique %>% length # 239
triplet.promoter.ewas$regionID %>% unique %>% length # 36
triplet.promoter.ewas$target %>% unique %>% length # 31

file.promoter <- file.path(path.mathReg, "promoter/ADNI_and_remap_promoter_analysis_using_TF_es_dorothea_all_triplet.csv")
dir.create(dirname(file.promoter),recursive = TRUE,showWarnings = FALSE)
dnam <- MethReg:::map_probes_to_regions(residuals.matched.met,genome = "hg19",arrayType = "EPIC")

cores <- 1
results.promoter.analysis <- 
  triplet.promoter.ewas %>% 
  cor_tf_target_gene(
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores
  ) %>% cor_dnam_target_gene(
    dnam = dnam,
    exp = residuals.matched.exp,
    cores = cores,
    filter.results = FALSE, 
    min.cor.estimate = 0.2,
    min.cor.pval = 0.05
  ) %>%  interaction_model(
    dnam = dnam,
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = FALSE,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE
  ) %>%  stratified_model(
    dnam = dnam,
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores
  )

results.promoter.analysis <- results.promoter.analysis %>% add_annot_cpgs() %>%  
  add_percent_zero_q1_q4(dnam = residuals.matched.met, exp = residuals.matched.exp) %>%
  update_met_IQR(dnam = residuals.matched.met)

results.promoter.analysis$RLM_TF_fdr <- p.adjust(results.promoter.analysis$RLM_TF_pvalue,method = "fdr")
results.promoter.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.promoter.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
results.promoter.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.promoter.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")

readr:::write_csv(
  x = results.promoter.analysis,
  file = file.promoter
)
results.promoter.analysis.sig.fdr.int <- results.promoter.analysis %>%
  dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)

readr:::write_csv(
  x = results.promoter.analysis.sig.fdr.int,
  file = gsub("all_triplet","sig_fdr_int_triplet",file.promoter)
)

results.promoter.analysis.sig.int <- results.promoter.analysis %>%
  dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)

readr:::write_csv(
  x = results.promoter.analysis.sig.int,
  file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.promoter)
)

results.promoter.analysis.sig <- results.promoter.analysis %>%
  filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))

readr:::write_csv(
  x = results.promoter.analysis.sig,
  file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.promoter)
)

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.distal.ewas <- create_triplet_distance_based(
  region = EPIC.hg19[rownames(se.selected.cpgs),],
  genome = "hg19",
  target.method =  "nearby.genes",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500,
  target.rm.promoter.regions.from.distal.linking = TRUE
) 
triplet.distal.ewas <- triplet.distal.ewas[!is.na(triplet.distal.ewas$TF),]
dim(triplet.distal.ewas) # 57730 triplets
triplet.distal.ewas <- triplet.distal.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es))
dim(triplet.distal.ewas) # 39360  triplets
triplet.distal.ewas$probeID <- names(EPIC.hg19)[match(triplet.distal.ewas$regionID,make_names_from_granges(EPIC.hg19))]


triplet.distal.ewas$TF %>% unique %>% length # 265
triplet.distal.ewas$regionID %>% unique %>% length # 92
triplet.distal.ewas$target %>% unique %>% length # 862


file.distal <- file.path(path.mathReg, "distal/ADNI_and_remap_distal_analysis_using_TF_es_dorothea_all_triplet.csv")
dir.create(dirname(file.distal),recursive = TRUE)
cores <- 4
results.distal.analysis <- 
  triplet.distal.ewas %>% 
  cor_tf_target_gene(
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores
  ) %>% cor_dnam_target_gene(
    dnam = dnam,
    exp = residuals.matched.exp,
    cores = cores,
    filter.results = FALSE, 
    min.cor.estimate = 0.2,
    min.cor.pval = 0.05
  ) %>%  interaction_model(
    dnam = dnam,
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = FALSE,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE
  ) %>%  stratified_model(
    dnam = dnam,
    exp = residuals.matched.exp,
    tf.activity.es = rnaseq.tf.es,
    cores = cores
  )

#results.distal.analysis <- results.distal.analysis[results.distal.analysis$probeID %in% cpgs.iqr.higher.0.03,]
results.distal.analysis <- results.distal.analysis %>% add_annot_cpgs() %>%  
  add_percent_zero_q1_q4(dnam = residuals.matched.met, exp = residuals.matched.exp) %>%
  update_met_IQR(dnam = residuals.matched.met)

results.distal.analysis$RLM_TF_fdr <- p.adjust(results.distal.analysis$RLM_TF_pvalue,method = "fdr")
results.distal.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.distal.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
results.distal.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.distal.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")

readr:::write_csv(
  x = results.distal.analysis,
  file = file.distal
)

results.distal.analysis.sig.fdr.int <- results.distal.analysis %>%
  dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)
results.distal.analysis.sig.fdr.int <- results.distal.analysis.sig.fdr.int[order(results.distal.analysis.sig.fdr.int$`RLM_DNAmGroup:TF_fdr`),]


readr:::write_csv(
  x = results.distal.analysis.sig.fdr.int,
  file = gsub("all_triplet","sig_fdr_int_triplet",file.distal)
)


results.distal.analysis.sig.int <- results.distal.analysis %>%
  dplyr::filter(results.distal.analysis$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)

readr:::write_csv(
  x = results.distal.analysis.sig.int,
  file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.distal)
)

results.distal.analysis.sig <- results.distal.analysis %>%
  filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))

readr:::write_csv(
  x = results.distal.analysis.sig,
  file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.distal)
)

#-------------------------------------------------------------------------------------------
# Save table
#------------------------------------------------------------------------------------------
writexl::write_xlsx(
  list(
    "Methreg_promoter.analysis_all" = results.promoter.analysis,
    "Methreg_distal.analysis_all" = results.distal.analysis
  ),
  path = file.path(path.mathReg,"cpg/MethReg_distal_promoter.xlsx")
)

#-------------------------------------------------------------------------------------------
# Plot triplets
#------------------------------------------------------------------------------------------
triplets <- plyr::rbind.fill(results.promoter.analysis,results.distal.analysis)
triplets <- triplets[which(triplets$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05),]

plots <- plot_interaction_model(
  triplet.results = triplets, 
  dnam = dnam, 
  label.dnam = "residuals",
  label.exp = "residuals",
  exp = residuals.matched.exp,
  tf.activity.es = rnaseq.tf.es,
  genome = "hg19"
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.mathReg, paste0("plots/Distal_promoter_RLM_DNAmGroup_TF_triplet_stage_wise_adj_pvalue_less_than_005.pdf")),
  plot = plots.one.page,
  width = 11,
  height = 13
)  

