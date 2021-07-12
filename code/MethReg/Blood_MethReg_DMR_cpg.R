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
# MethReg analysis
# target gene ~ TF_activity (GSVA) + median CpG at DMR + median CpG at DMR * TF
#-----------------------------------------------------------------------------
path.mathReg <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/analysis_results/methReg_DMR_by_cpg"

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
)

#-----------------------------------------------------------------------------
# GET ADNI data
#-----------------------------------------------------------------------------
load("~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/datasets/Aux/ADNI_matched_rna_dnam.rda")
colnames(adni.se) <- paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE)
cpgs <- combp_AD_vs_CN$Probes %>% stringr::str_split(";") %>% unlist
se.selected <- adni.se[rownames(adni.se) %in% cpgs,]
#-----------------------------------------------------------------------------
# get residuals 
#-----------------------------------------------------------------------------
metadata.dnam <- colData(adni.se)[,c("CD8T","Mono","Neutro","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- get_residuals(
  data = log2(assay(se.selected)) / (1 - assay(se.selected)), # m-values
  metadata.samples = metadata.dnam,
  cores = 4
)

Affy_Plate <- ADNI_Gene_Expression_Metadata[7,-c(1:3)] %>% as.numeric()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

metadata.exp <- colData(adni.se)[,c("age_at_visit","PTGENDER")]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


library(xCell)
aux <- expression.matrix
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[5] <- "B_cells"
colnames(metadata.exp)[7] <- "CD4_T_cells"
colnames(metadata.exp)[8] <- "CD8_T_cells"

sample.to.rm <- rownames(metadata.exp)[is.na(metadata.exp$RIN)]
residuals.matched.exp <- get_residuals(
  data = log2(expression.matrix),
  metadata.samples = metadata.exp,
  cores = 4
)
residuals.matched.exp <- residuals.matched.exp[,which(colnames(residuals.matched.exp) != sample.to.rm)]
metadata.exp <- metadata.exp[which(rownames(metadata.exp) != sample.to.rm),]
residuals.matched.met <- residuals.matched.met[,which(colnames(residuals.matched.met) != sample.to.rm)]
metadata.dnam <- metadata.dnam[which(rownames(metadata.dnam) != sample.to.rm),]

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
library(GSVA)
file <-  "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/datasets/Aux/exp_tf_es_gsva.rda"
load(file)
all(colnames(exp.tf.es.gsva) == colnames(residuals.matched.exp))

exp.tf.es.gsva.ensg <- exp.tf.es.gsva
rownames(exp.tf.es.gsva.ensg) <- MethReg:::map_symbol_to_ensg(rownames(exp.tf.es.gsva))
exp.tf.es.gsva.ensg <- exp.tf.es.gsva.ensg[!is.na(rownames(exp.tf.es.gsva.ensg)),]

save(
  exp.tf.es.gsva,
  exp.tf.es.gsva.ensg,
  residuals.matched.exp,
  metadata.dnam,
  metadata.exp,
  residuals.matched.met,
  file = "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/datasets/Aux/methReg_data_input_DMR_by_cpg.rda"
)
load("~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/datasets/Aux/methReg_data_input_DMR_by_cpg.rda")

#-------------------------------------------------------------------------------
# Analysis
#-------------------------------------------------------------------------------
# Get triplets using remap
library(MethReg)
library(ReMapEnrich)
dir.base <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/"
dir.data.aux <- file.path(dir.base,"datasets/Aux/") 
remapCatalog2018hg19 <- downloadRemapCatalog(dir.data.aux, assembly = "hg19")
remapCatalog <- bedToGranges(remapCatalog2018hg19)

#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------

update_met_IQR <- function(results, dnam){
  iqr <- calculate_IQR(dnam)
  results$met.IQR <- iqr$IQR[match(results$probeID, iqr$ID)]
  results <- results %>% dplyr::relocate(dplyr::contains("IQR"), .after = last_col())
  return(results)
}


add_annot_cpgs <- function(results){
  
  aux <-   separate_rows(combp_AD_vs_CN,"Probes",sep = ";")
  results <- cbind(
    results,
    aux[
      match(results$probeID,aux$Probes),
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
    unique(results[,c("probeID","target")]),
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
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")

triplet.promoter.ewas <- create_triplet_distance_based(
  region =  EPIC.hg19[rownames(residuals.matched.met),],
  genome = "hg19",
  target.method =  "genes.promoter.overlap",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
triplet.promoter.ewas <- triplet.promoter.ewas[!is.na(triplet.promoter.ewas$TF),]
nrow(triplet.promoter.ewas) # 2195 triplets

triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$TF %in% rownames(exp.tf.es.gsva.ensg))
triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp))
triplet.promoter.ewas$probeID <- names(EPIC.hg19)[match(triplet.promoter.ewas$regionID,make_names_from_granges(EPIC.hg19))]

nrow(triplet.promoter.ewas) # 1233 triplets

triplet.promoter.ewas$TF %>% unique %>% length # 196
triplet.promoter.ewas$regionID %>% unique %>% length # 5
triplet.promoter.ewas$target %>% unique %>% length # 4

file.promoter <- file.path(path.mathReg, "promoter/ADNI_and_remap_promoter_analysis_using_TF_es_gsva_all_triplet.csv")
dir.create(dirname(file.promoter),showWarnings = FALSE,recursive = TRUE)
dnam <- MethReg:::map_probes_to_regions(residuals.matched.met,genome = "hg19",arrayType = "EPIC")


if (!file.exists(file.promoter)) {
  cores <- 1
  results.promoter.analysis <- 
    triplet.promoter.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = exp.tf.es.gsva.ensg,
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
      tf.activity.es = exp.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = dnam,
      exp = residuals.matched.exp,
      tf.activity.es = exp.tf.es.gsva.ensg,
      cores = cores
    )
  
  results.promoter.analysis <- results.promoter.analysis %>% add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = residuals.matched.met, exp = residuals.matched.exp) %>%
    update_met_IQR(dnam = dnam)
  
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
} else {
  results.promoter.analysis <- readr::read_csv(
    file = file.promoter,
    col_types = c(`quant_triplet_stage_wise_adj_pval_metGrp:es.tf` = "n")
  )
  
  # results.promoter.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.promoter.analysis$probeID %in% enhancer.probes
}

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.distal.ewas <- create_triplet_distance_based(
  region =  EPIC.hg19[rownames(residuals.matched.met),],
  genome = "hg19",
  target.method =  "nearby.genes",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500,
  target.rm.promoter.regions.from.distal.linking = TRUE
) 
triplet.distal.ewas <- triplet.distal.ewas[!is.na(triplet.distal.ewas$TF),]
dim(triplet.distal.ewas) # 6280 triplets
triplet.distal.ewas <- triplet.distal.ewas %>% dplyr::filter(.data$TF %in% rownames(exp.tf.es.gsva.ensg))
dim(triplet.distal.ewas) # 4300  triplets
triplet.distal.ewas$probeID <- names(EPIC.hg19)[match(triplet.distal.ewas$regionID,make_names_from_granges(EPIC.hg19))]


triplet.distal.ewas$TF %>% unique %>% length # 108
triplet.distal.ewas$regionID %>% unique %>% length # 2
triplet.distal.ewas$target %>% unique %>% length # 20


file.distal <- file.path(path.mathReg, "distal/ADNI_and_remap_distal_analysis_using_TF_es_gsva_all_triplet.csv")
dir.create(dirname(file.distal),showWarnings = FALSE,recursive = TRUE)
if (!file.exists(file.distal)) {
  cores <- 1
  results.distal.analysis <- 
    triplet.distal.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = exp.tf.es.gsva.ensg,
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
      tf.activity.es = exp.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = dnam,
      exp = residuals.matched.exp,
      tf.activity.es = exp.tf.es.gsva.ensg,
      cores = cores
    )
  
  #results.distal.analysis <- results.distal.analysis[results.distal.analysis$probeID %in% cpgs.iqr.higher.0.03,]
  results.distal.analysis <- results.distal.analysis %>% add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = residuals.matched.met, exp = residuals.matched.exp) %>%
    update_met_IQR(dnam = dnam)
  
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
}  else {
  results.distal.analysis <- readr::read_csv(
    file = file.distal,
    col_types = c(`quant_triplet_stage_wise_adj_pval_metGrp:es.tf` = "n")
  )
  
  results.distal.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.distal.analysis$probeID %in% enhancer.probes
  
}

#-------------------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    "Methreg_promoter.analysis_all" = results.promoter.analysis,
    "Methreg_distal.analysis_all" = results.distal.analysis
  ),
  path = file.path(path.mathReg,"MethReg_regulon_distal_promoter.xlsx")
)

#-------------------------------------------------------------------------------------------
# Plot triplets
#------------------------------------------------------------------------------------------
triplets <- rbind(results.promoter.analysis,results.distal.analysis)
triplets <- triplets[which(triplets$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05),]

plots <- plot_interaction_model(
  triplet.results = triplets, 
  dnam = dnam, 
  label.dnam = "residuals",
  label.exp = "residuals",
  exp = residuals.matched.exp,
  tf.activity.es = exp.tf.es.gsva.ensg,
  genome = "hg19"
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
file.plot <- file.path(path.mathReg, paste0("plots/Regulon_distal_promoter_RLM_DNAmGroup_TF_triplet_stage_wise_adj_pvalue_less_than_0_05.pdf"))
dir.create(dirname(file.plot),showWarnings = FALSE,recursive = TRUE)
ggplot2::ggsave(
  filename = file.plot,
  plot = plots.one.page,
  width = 11,
  height = 13
)  

