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

#-------------------------------------------------------------------------------
# Libs
#-------------------------------------------------------------------------------
library(dplyr)
library(meta)
library(ggplot2)
library(cowplot)
library(grid)
#-------------------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------------------
dir.plot <- "plots/forest_plot/"
dir.draft <- "DRAFT-TABLES_FIGURES_4-17-2021"
dir.analysis.results <- "analysis_results"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-------------------------------------------------------------------------------
# cpgs
#-------------------------------------------------------------------------------
cpgs <- readxl::read_xlsx(
  file.path(dir.draft, "_Main Table 1 FDR sig CpGs in AD vs CN blood samples meta-analysis.xlsx"),
  skip = 3
)

# AIBL & ADNI
res <- readr::read_csv(
  file.path(dir.analysis.results, "meta_analysis/Logistic_regression_model/AD_vs_CN/meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated_pvalue_cut_off_0_05.csv")
)

plots <- plyr::llply(
  .data = cpgs$cpg,
  .fun = function(x){
    
    row <- res %>% dplyr::filter(cpg == x) 

    adni <- data.frame (
      cpg = x,
      study = "ADNI", 
      estimate = row$ADNI_AD_vs_CN_Estimate.bacon / 100, 
      se = row$ADNI_AD_vs_CN_StdErr.bacon / 100
    )
    aibl <- data.frame (
      cpg = x,
      study = "AIBL",
      estimate = row$AIBL_AD_vs_CN_Estimate.bacon / 100, 
      se = row$AIBL_AD_vs_CN_StdErr.bacon / 100
    )
    
    two <- rbind (adni, aibl)
    
    meta_bin <- metagen (
      TE = estimate, 
      seTE = se, 
      comb.fixed = TRUE,
      comb.random = FALSE,
      prediction = FALSE,
      sm = "OR",
      data = two, 
      studlab = study, 
      title = paste0("CpG ",x)
    )
    grid.newpage()
    forest (
      meta_bin, 
      comb.fixed = TRUE, 
      comb.random = FALSE, 
      print.tau2 = TRUE,
      sortvar = TE, # sorting variable
      digits.sd = 2,
      digits.se = 2,
      squaresize = 1, # size of the squares in the plot
      digits = 2,
      ddigits = 3, 
      digits.tau2 = 2,
      col.square = "blue",
      col.inside.fixed = "red",
      col.diamond = "red", 
      text.random = "Fixed effect model",
      leftcols = c("studlab"),
      rightcols = c("ci"),
      print.I2.ci = TRUE,
      #layout = "RevMan5", # Or JAMA
      colgap.forest.left = unit(4,"cm")
    )
    grid.grab()
  })

merged.cpg <- ggpubr::ggarrange(
  plotlist = plots,
  ncol = 1,
  nrow = 5,
  labels = paste0(
    "TOP ", 1:5," cpg ", 
    c(cpgs$cpg),
    " (", gsub("\\(-[0-9].*\\)|\\(\\+[0-9].*\\)","",cpgs$GREAT_annotation) %>% stringr::str_trim(),")"
  )
)
ggsave(
  plot = merged.cpg,
  filename = file.path(dir.plot, "Forest_plot_cpg.pdf")
)

