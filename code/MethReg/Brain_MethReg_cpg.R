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

#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Path and libraries
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
library(MethReg)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(SummarizedExperiment)
library(coMethDMR)

path.mathReg <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/analysis_results/methReg/Brain/"
path.dropbox <- "~/TBL Dropbox/Tiago Silva/"
path.project <- file.path(path.dropbox,"/PanCancer/MethReg-useCase/")
path.data <- file.path(path.project,"data")
path.usecase <- file.path(path.project,"/UseCase2_rosmap_remap//")
path.tables.promoter <- file.path(path.mathReg,"/promoter")
path.tables.distal <- file.path(path.mathReg,"/distal")
path.ewas.analysis <-  file.path(path.dropbox,"/coMethDMR_metaAnalysis/")

dir.create(path.tables.promoter,recursive = T,showWarnings = F)
dir.create(path.tables.distal,recursive = T,showWarnings = F)

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
cpgs.prioritized <- readxl::read_xlsx(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0
)
cpgs.prioritized  <- cpgs.prioritized[[5]] %>% na.omit() %>% as.character
length(cpgs.prioritized) # 97

dmrs.prioritized <- readxl::read_xlsx(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0,sheet = 2
) 
dmrs.prioritized <- dmrs.prioritized[[4]] %>% na.omit() %>% as.character
cpgs.in.dmrs.prioritized <- coMethDMR::GetCpGsInRegion(dmrs.prioritized)

cpgs.all <- c(
  cpgs.prioritized,
  cpgs.in.dmrs.prioritized
) %>% unique

#-----------------------------------------------------------------------------
# MethReg analysis
# target gene ~ TF_activity (dorothea) + CpG + CpG * TF
#-----------------------------------------------------------------------------


#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Data from other projects
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
london.blood <- readr::read_csv("~/Dropbox (BBSR)/Tiago Silva/coMethDMR_metaAnalysis/code_validation/Meta_analysis_code/London_blood_brain_correlation_results/London_blood_brain_correlation_cpgs.csv")
colnames(london.blood) <- paste0("London_blood_",colnames(london.blood))
colnames(london.blood)[1] <- "probeID"

#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Aux functions
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-

add_annot_cpgs <- function(results){
  results$UCSC_RefGene_Group <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[results$probeID,]$UCSC_RefGene_Group
  results$UCSC_RefGene_Name <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[results$probeID,]$UCSC_RefGene_Name
  results$Relation_to_Island <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC[results$probeID,]$Relation_to_Island
  results
}


update_met_IQR <- function(results, dnam){
  iqr <- MethReg:::calculate_IQR(dnam)
  results$res.met.IQR <- results$met.IQR
  results$met.IQR <- iqr$IQR[match(results$regionID, iqr$ID)]
  results <- results %>% dplyr::relocate(dplyr::contains("IQR"), .after = last_col())
  return(results)
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
# Read Rosmap data
#-------------------------------------------------------------------------------
load(file = file.path(path.usecase,"data/matched_normalized_and_residuals_data.rda"))

iqr <- MethReg:::calculate_IQR(matched.dnam)
cpgs.iqr.higher.0.03 <- iqr$ID[iqr$IQR > 0.03]

resid.met.cpg <- residuals.matched.met[rownames(residuals.matched.met) %in% cpgs.all,]
all(colnames(residuals.matched.exp) == colnames(resid.met.cpg))

hm450.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "450k")
rownames(resid.met.cpg) <- MethReg:::make_names_from_granges(hm450.hg19[rownames(resid.met.cpg) ,])

rnaseq.tf.es.gsva.ensg <- rnaseq.tf.es.gsva
rownames(rnaseq.tf.es.gsva.ensg) <- MethReg:::map_symbol_to_ensg(rownames(rnaseq.tf.es.gsva))
rnaseq.tf.es.gsva.ensg <- rnaseq.tf.es.gsva.ensg[!is.na(rownames(rnaseq.tf.es.gsva.ensg)),]


matched.dnam.with.regions <- matched.dnam
hm450.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "450k")
rownames(matched.dnam.with.regions) <- MethReg:::make_names_from_granges(hm450.hg19[rownames(matched.dnam.with.regions) ,])

#-------------------------------------------------------------------------------
# Get triplets using remap
#-------------------------------------------------------------------------------
library(ReMapEnrich)
remapCatalog2018hg19 <- downloadRemapCatalog(path.data, assembly = "hg19")
remapCatalog <- bedToGranges(remapCatalog2018hg19)

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.promoter.ewas <- create_triplet_distance_based(
  region = hm450.hg19[cpgs.all,],
  genome = "hg19",
  target.method =  "genes.promoter.overlap",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
dim(triplet.promoter.ewas) # 69928 triplets
triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es.gsva.ensg))
dim(triplet.promoter.ewas) # 57715 triplets
triplet.promoter.ewas$probeID <- names(hm450.hg19)[match(triplet.promoter.ewas$regionID,make_names_from_granges(hm450.hg19))]

triplet.promoter.ewas$TF %>% unique %>% length # 286
triplet.promoter.ewas$regionID %>% unique %>% length # 1404
triplet.promoter.ewas$target %>% unique %>% length # 1164


file.promoter <- file.path(path.tables.promoter, "ROSMAP_and_remap_promoter_analysis_using_TF_es_gsva_all_triplet.csv")

if (!file.exists(file.promoter)) {
  cores <- 10
  results.promoter.analysis <- 
    triplet.promoter.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    ) %>% cor_dnam_target_gene(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      cores = cores,
      filter.results = FALSE, 
      min.cor.estimate = 0.2,
      min.cor.pval = 0.05
    ) %>%  interaction_model(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    )
  
  results.promoter.analysis <- results.promoter.analysis %>% add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = matched.dnam, exp = matched.exp.log2) %>%
    update_met_IQR(dnam = matched.dnam.with.regions)
  
  results.promoter.analysis$RLM_TF_fdr <- p.adjust(results.promoter.analysis$RLM_TF_pvalue,method = "fdr")
  results.promoter.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.promoter.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
  results.promoter.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.promoter.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")
  
  readr:::write_csv(
    x = results.promoter.analysis,
    file = file.promoter
  )
  results.promoter.analysis.sig.fdr.int <- results.promoter.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)
  
  # results.promoter.analysis.sig.fdr.int.with.blood <- merge(results.promoter.analysis.sig.fdr.int, london.blood)
  #readr:::write_csv(
  #  x = results.promoter.analysis.sig.fdr.int.with.blood,
  #  file = gsub("all_triplet","sig_fdr_int_triplet_with_london_blood",file.promoter)
  # )
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
} else {
  results.promoter.analysis <- readr::read_csv(
    file = file.promoter,
    col_types = c(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` = "n")
  )
  
  results.promoter.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.promoter.analysis$probeID %in% enhancer.probes
}


results.promoter.analysis.sig.int <- results.promoter.analysis.sig.int[order(results.promoter.analysis.sig.int$`quant_triplet_stage_wise_adj_pval_metGrp:es.tf`),]

#~~~~~~~~~~~~~
# Plots      
#~~~~~~~~~~~~~
plots.promoter <- plot_interaction_model(
  triplet.results = results.promoter.analysis.sig.int, 
  dnam = resid.met.cpg, 
  label.dnam = "residuals",
  label.exp = "residuals",
  exp = residuals.matched.exp,
  tf.activity.es = rnaseq.tf.es.gsva.ensg,
  genome = "hg19"
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots.promoter, nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.tables.promoter, "ROSMAP_Promoter_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf"),
  plot = plots.one.page,
  width = 12,
  height = 7
)  


#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.distal.ewas <- create_triplet_distance_based(
  region = hm450.hg19[cpgs.all,],
  genome = "hg19",
  target.method =  "nearby.genes",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500,
  target.rm.promoter.regions.from.distal.linking = TRUE
) 
dim(triplet.distal.ewas) # 1,422,184 triplets
triplet.distal.ewas <- triplet.distal.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es.gsva.ensg))
dim(triplet.distal.ewas) # 946,848  triplets
triplet.distal.ewas$probeID <- names(hm450.hg19)[match(triplet.distal.ewas$regionID,make_names_from_granges(hm450.hg19))]


triplet.distal.ewas$TF %>% unique %>% length # 288
triplet.distal.ewas$regionID %>% unique %>% length # 1756
triplet.distal.ewas$target %>% unique %>% length # 9,934


file.distal <- file.path(path.tables.distal, "ROSMAP_and_remap_distal_analysis_using_TF_es_gsva_all_triplet.csv")

if (!file.exists(file.distal)) {
  cores <- 4
  results.distal.analysis <- 
    triplet.distal.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    ) %>% cor_dnam_target_gene(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      cores = cores,
      filter.results = FALSE, 
      min.cor.estimate = 0.2,
      min.cor.pval = 0.05
    ) %>%  interaction_model(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = resid.met.cpg,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    )
  
  #results.distal.analysis <- results.distal.analysis[results.distal.analysis$probeID %in% cpgs.iqr.higher.0.03,]
  results.distal.analysis <- results.distal.analysis %>% 
    add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = matched.dnam, exp = matched.exp.log2)  %>%
    update_met_IQR(dnam = matched.dnam.with.regions)
  
  
  results.distal.analysis$RLM_TF_fdr <- p.adjust(results.distal.analysis$RLM_TF_pvalue,method = "fdr")
  results.distal.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.distal.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
  results.distal.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.distal.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")
  
  readr:::write_csv(
    x = results.distal.analysis,
    file = file.distal
  )
  
  results.distal.analysis.sig.fdr.int <- results.distal.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)
  
  # results.distal.analysis.sig.fdr.int.with.blood <- merge(results.distal.analysis.sig.fdr.int, london.blood)
  #readr:::write_csv(
  #  x = results.distal.analysis.sig.fdr.int.with.blood,
  #  file = gsub("all_triplet","sig_fdr_int_triplet_with_london_blood",file.distal)
  #)
  
  
  readr:::write_csv(
    x = results.distal.analysis.sig.fdr.int,
    file = gsub("all_triplet","sig_fdr_int_triplet",file.distal)
  )
  
  results.distal.analysis.sig.int <- results.distal.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)
  
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
}  else {
  results.distal.analysis <- readr::read_csv(
    file = file.distal,
    col_types = c(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` = "n")
  )
  
  results.distal.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.distal.analysis$probeID %in% enhancer.probes
}


results.distal.analysis.sig.int <- results.distal.analysis.sig.int[order(results.distal.analysis.sig.int$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue`),]


#~~~~~~~~~~~~~
# Plots      
#~~~~~~~~~~~~~
plots.distal <- plot_interaction_model(
  triplet.results = results.distal.analysis.sig.int, 
  dnam = resid.met.cpg, 
  exp = residuals.matched.exp,
  tf.activity.es = rnaseq.tf.es.gsva.ensg,
  label.dnam = "residuals",
  label.exp = "residuals",
  genome = "hg19",
  add.tf.vs.exp.scatter.plot = TRUE
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots.distal, nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.tables.distal, "ROSMAP_Distal_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf"),
  plot = plots.one.page,
  width = 12,
  height = 10
)  

#-------------------------------------------------------------------------------------------
# Save table
#------------------------------------------------------------------------------------------
writexl::write_xlsx(
  list(
    "ROSMAP_Methreg_promoter" = results.promoter.analysis,
    "ROSMAP_Methreg_distal" = results.distal.analysis
  ),
  path = file.path(path.mathReg,"ROSMAP_MethReg_distal_promoter.xlsx")
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(c(plots.promoter, plots.distal), nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.mathReg, "ROSMAP_Promoter_Distal_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf"),
  plot = plots.one.page,
  width = 12,
  height = 10
)  

