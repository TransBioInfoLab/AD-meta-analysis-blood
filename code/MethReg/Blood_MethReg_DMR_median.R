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
library(coMethDMR)
library(MethReg)
library(ReMapEnrich)

path.mathReg <- "AD-meta-analysis-blood-samples/analysis_results/methReg_DMR_median"
dir.create(path.mathReg, recursive = TRUE, showWarnings = FALSE)
#-----------------------------------------------------------------------------
# MethReg analysis
# target gene ~ TF_activity (dorothea) + median CpG at DMR + median CpG at DMR * TF
#-----------------------------------------------------------------------------

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
  )
}) # 19 DMRs

combp_AD_vs_CN$Probes <- sapply(
  combp_AD_vs_CN.probes,
  FUN = function(x) paste(x,collapse = ";")
)

#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(
  path = "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 3 Top 10 prioritized-CpGs_and_DMRs-crossTissue_brain_blood-V2.xlsx",
  sheet = 1, 
  skip = 17
)
prioritized.dmrs <- prioritized.dmrs$DMRs %>% na.omit %>% as.character
length(prioritized.dmrs) # 10


prioritized.dmrs.probes <- plyr::alply(prioritized.dmrs,.margins = 1,.fun = function(dmr) {
  GetCpGsInRegion(
    dmr,
    arrayType = "EPIC"
  )
}) # 19 DMRs

prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(
  prioritized.dmrs.probes,
  FUN = function(x) paste(x,collapse = ";")
)

regions <- rbind(combp_AD_vs_CN[,c("DMR","Probes")],prioritized.dmrs[,c("DMR","Probes")])

load("datasets/Aux/ADNI_matched_rna_dnam_residuals_DMR.rda")

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
    combp_AD_vs_CN[
      match(results$regionID,combp_AD_vs_CN$DMR),
      c(
        "Relation to Island",
        "UCSC RefGene Name",
        "UCSC RefGene Group",
        "GREAT Annotation",
        "Chromatin State"
      )
    ]
  )
  results
}


add_percent_zero_q1_q4 <- function(results, dnam, exp){
  
  aux <- plyr::adply(
    unique(results[,c("regionID","target")]),
    .margins = 1,
    .fun = function(row) {
      rna.target <- exp[rownames(exp) == row$target, , drop = FALSE]
      met <- dnam[rownames(dnam) == as.character(row$regionID), ]
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
triplet.promoter.ewas <- create_triplet_distance_based(
  region = rownames(residuals.matched.met),
  genome = "hg19",
  target.method =  "genes.promoter.overlap",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
triplet.promoter.ewas <- triplet.promoter.ewas[!is.na(triplet.promoter.ewas$TF),]
nrow(triplet.promoter.ewas) # 1076 triplets

triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es))
triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))

nrow(triplet.promoter.ewas) # 587 triplets

triplet.promoter.ewas$TF %>% unique %>% length # 211
triplet.promoter.ewas$regionID %>% unique %>% length # 9
triplet.promoter.ewas$target %>% unique %>% length # 8

file.promoter <- file.path(path.mathReg, "promoter/ADNI_and_remap_promoter_analysis_using_TF_es_dorothea_all_triplet.csv")
dir.create(dirname(file.promoter),showWarnings = FALSE,recursive = TRUE)
dnam <- residuals.matched.met

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
  region = rownames(se.selected.median),
  genome = "hg19",
  target.method =  "nearby.genes",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500,
  target.rm.promoter.regions.from.distal.linking = TRUE
) 
triplet.distal.ewas <- triplet.distal.ewas[!is.na(triplet.distal.ewas$TF),]
dim(triplet.distal.ewas) # 6280 triplets
triplet.distal.ewas <- triplet.distal.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es))
dim(triplet.distal.ewas) # 4300  triplets

triplet.distal.ewas$TF %>% unique %>% length # 166
triplet.distal.ewas$regionID %>% unique %>% length # 6
triplet.distal.ewas$target %>% unique %>% length # 60


file.distal <- file.path(path.mathReg, "distal/ADNI_and_remap_distal_analysis_using_TF_es_dorothea_all_triplet.csv")
dir.create(dirname(file.distal),showWarnings = FALSE,recursive = TRUE)
cores <- 1
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


#-------------------------------------------------------------------------------
# Regulon analysis, triplets using remap
#-------------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
triplet.regulon.ewas <- create_triplet_regulon_based(
  region = rownames(se.selected.median),
  tf.target = regulons.dorothea,
  genome = "hg19",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
triplet.regulon.ewas <- triplet.regulon.ewas[!is.na(triplet.regulon.ewas$TF),]
dim(triplet.regulon.ewas) # 164 triplets
triplet.regulon.ewas <- triplet.regulon.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es))
dim(triplet.regulon.ewas) # 164   triplets
triplet.regulon.ewas$probeID <- names(EPIC.hg19)[match(triplet.regulon.ewas$regionID,make_names_from_granges(EPIC.hg19))]

triplet.regulon.ewas$TF %>% unique %>% length # 192
triplet.regulon.ewas$regionID %>% unique %>% length # 53
triplet.regulon.ewas$target %>% unique %>% length # 342

file.regulon <- file.path(path.mathReg, "regulon/ADNI_and_remap_regulon_analysis_using_TF_es_dorothea_all_triplet.csv")
dir.create(dirname(file.regulon),showWarnings = FALSE,recursive = TRUE)
cores <- 1
results.regulon.analysis <- 
  triplet.regulon.ewas %>% 
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

# results.regulon.analysis <- results.regulon.analysis[results.regulon.analysis$probeID %in% cpgs.iqr.higher.0.03,]
results.regulon.analysis <- results.regulon.analysis %>% 
  add_annot_cpgs() %>%  
  add_percent_zero_q1_q4(dnam = residuals.matched.met, exp = residuals.matched.exp) %>%
  update_met_IQR(dnam = residuals.matched.met)


results.regulon.analysis$RLM_TF_fdr <- p.adjust(results.regulon.analysis$RLM_TF_pvalue,method = "fdr")
results.regulon.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.regulon.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
results.regulon.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.regulon.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")


readr:::write_csv(
  x = results.regulon.analysis,
  file = file.regulon
)
results.regulon.analysis.sig.fdr.int <- results.regulon.analysis %>%
  dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)

readr:::write_csv(
  x = results.regulon.analysis.sig.fdr.int,
  file = gsub("all_triplet","sig_fdr_int_triplet",file.regulon)
)


results.regulon.analysis.sig.int <- results.regulon.analysis %>%
  dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue`< 0.05)

readr:::write_csv(
  x = results.regulon.analysis.sig.int,
  file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.regulon)
)

results.regulon.analysis.sig <- results.regulon.analysis %>%
  filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))

readr:::write_csv(
  x = results.regulon.analysis.sig,
  file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.regulon)
)


#-------------------------------------------------------------------------------------------
# Final table
#------------------------------------------------------------------------------------------
writexl::write_xlsx(
  list(
    "Methreg_promoter.analysis_all" = results.promoter.analysis,
    "Methreg_distal.analysis_all" = results.distal.analysis,
    "Methreg_regulon.analysis_all" = results.regulon.analysis
  ),
  path = file.path(path.mathReg,"MethReg_regulon_distal_promoter.xlsx")
)

#-------------------------------------------------------------------------------------------
# Plot triplets
#------------------------------------------------------------------------------------------
triplets <- plyr::rbind.fill(results.promoter.analysis,results.distal.analysis,results.regulon.analysis)
triplets <- triplets[which(triplets$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.2),]

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
file.plot <- file.path(path.mathReg, paste0("plots/Regulon_distal_promoter_RLM_DNAmGroup_TF_triplet_stage_wise_adj_pvalue_less_than_0_2.pdf"))
dir.create(dirname(file.plot),showWarnings = FALSE,recursive = TRUE)
ggplot2::ggsave(
  filename = file.plot,
  plot = plots.one.page,
  width = 11,
  height = 13
)  

