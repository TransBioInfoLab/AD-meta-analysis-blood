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
dir.plots <-  "plots/manhattan_plot/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
# install_github("juliedwhite/miamiplot", build_vignettes = TRUE)
library(miamiplot)
library(dplyr)
library(dplyr)
source("code/manhattan_plot/ggmiami.R")

#-------------------------------------------------------------------------------
# Get brain meta-analysis results data
#-------------------------------------------------------------------------------
cpgs.blood <- readr::read_csv("analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv")
table(cpgs.blood$pVal.final.bacon < 1e-5)
df.blood <- cpgs.blood %>% dplyr::select(c("cpg","chr","pos","pVal.final.bacon","fdr.bacon","UCSC_RefGene_Name")) %>%
  dplyr::rename(
    pVal.final = pVal.final.bacon,
    fdr = fdr.bacon
  )
df.blood$chr <- as.numeric(as.character(gsub("chr", "", df.blood$chr)))
df.blood$source <- "Blood"

cpgs.blood.cutoff <- max(
  df.blood[df.blood$fdr < 0.05, "pVal.final"], na.rm = T
)

dmrs.blood <- readr::read_csv("analysis_results/combp/AD_vs_CN_output_annotated.csv")
dmrs.blood <- dmrs.blood[dmrs.blood$n_probes > 2,]
dmrs.blood.cpgs <- coMethDMR:::GetCpGsInAllRegion(dmrs.blood$region,arrayType = "EPIC") %>% unlist %>% unique()

#-------------------------------------------------------------------------------
# Get brain meta-analysis results data
#-------------------------------------------------------------------------------
cpgs.brain <- readr::read_csv("datasets/brain_meta_analysis/meta_analysis_single_cpg_df.csv.gz")
cpgs.brain <-  cpgs.brain %>% dplyr::select(cpg, seqnames, start, pVal.final, fdr) %>%
  dplyr::rename(pos = start, chr = seqnames) %>%
  mutate(
    chr = substring(chr, 4),
    source = "Brain"
  )
cpgs.brain$chr <- as.numeric(cpgs.brain$chr)
cpgs.brain$chr <- ifelse(
  is.na(cpgs.brain$chr), 23, cpgs.brain$chr
)
cpgs.brain <- cpgs.brain %>% dplyr::filter(chr != 23)
cpgs.brain$UCSC_RefGene_Name <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other$UCSC_RefGene_Name[match(cpgs.brain$cpg,rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19::Other))]

cpgs.brain.cutoff <- max(
  cpgs.brain[cpgs.brain$fdr < 0.05, "pVal.final"], na.rm = T
)

brain_blood_df <- plyr::rbind.fill(df.blood,cpgs.brain)

plot_data <- prep_miami_data(
  data = brain_blood_df, split_by = "source", split_at = "Blood", p = "pVal.final"
)

#-------------------------------------------------------------------------------
# Label
#-------------------------------------------------------------------------------
brain_3751_cpgs <- readxl::read_xlsx("datasets/brain_meta_analysis/Supplementary Data 1-24.xlsx",skip = 3)
brain.top.20.cpgs <- brain_3751_cpgs$cpg[which(brain_3751_cpgs$pVal.final %in% head(sort(brain_3751_cpgs$pVal.final),n = 20))]

brain_119_dmrs <- readxl::read_xlsx("datasets/brain_meta_analysis/Supplementary Data 1-24.xlsx",skip = 3,sheet = 2)
brain_top20_regions <- brain_119_dmrs$DMR[which(brain_119_dmrs$pVal.final %in% head(sort(brain_119_dmrs$pVal.final),n = 20))]
dmrs.brain.cpgs <- coMethDMR:::GetCpGsInAllRegion(brain_top20_regions,arrayType = "450k") %>% unlist %>% unique()
brain_label_chr <- c(brain.top.20.cpgs,dmrs.brain.cpgs) %>% unique()

blood.top.20.cpgs <- cpgs.blood$cpg[which(cpgs.blood$pVal.final.bacon %in% head(sort(cpgs.blood$pVal.final.bacon),n = 20))]
blood_label_chr <- c(blood.top.20.cpgs,dmrs.blood.cpgs) %>% unique()
blood_label_chr <- c(blood.top.20.cpgs) %>% unique()

brain_labels <- plot_data$lower %>%
  filter(cpg %in% brain_label_chr) %>%
  select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
brain_labels <- brain_labels[brain_labels$label != "",]
brain_labels <- brain_labels[!duplicated(brain_labels$label),]

blood_labels <- plot_data$upper %>%
  filter(cpg %in% blood_label_chr) %>%
  select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()


#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------
plot <- ggmiami2(
  data = brain_blood_df, 
  split_by = "source", 
  split_at = "Blood", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "Blood",
  lower_ylab = "Brain",
  genome_line = NULL,
  top_n_hits = 10,
  suggestive_line_bottom = cpgs.brain.cutoff,
  suggestive_line_upper = cpgs.blood.cutoff,
  upper_labels_df = blood_labels,
  lower_labels_df = brain_labels
)

ggplot2::ggsave(
  plot = plot, 
  filename = file.path(dir.plots,"brain_blood_manhattan_plot.pdf"),
  width = 10,
  height = 8
)

ggplot2::ggsave(
  plot = plot, 
  filename = file.path(dir.plots,"brain_blood_manhattan_plot.png"),
  width = 7,
  height = 4
)





