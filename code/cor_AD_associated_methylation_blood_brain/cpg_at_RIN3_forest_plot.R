setwd("C:/Users/lxw391/TBL Dropbox/Lily Wang/AD-meta-analysis-blood-samples/LW/ADNI_subgrp_analysis")
library(dplyr)
library (meta)
library(ggplot2)
# forest plot for cg05157625 at RIN3 gene


# AIBL & ADNI
res <- read.csv("C:/Users/lxw391/TBL Dropbox/Lily Wang/AD-meta-analysis-blood-samples/analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated_pvalue_cut_off_0_05.csv")
res <- readr::read_csv("~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated_pvalue_cut_off_0_05.csv")

one <- subset (
  res, 
  cpg == "cg05157625", 
  select = c(
    cpg, 
    ADNI_AD_vs_CN_Estimate.bacon,
    ADNI_AD_vs_CN_StdErr.bacon,
    AIBL_AD_vs_CN_Estimate.bacon,
    AIBL_AD_vs_CN_StdErr.bacon
  )
)

adni <- data.frame (
  cpg = "cg05157625", 
  study = "ADNI", 
  estimate = one$ADNI_AD_vs_CN_Estimate.bacon / 100, 
  se = one$ADNI_AD_vs_CN_StdErr.bacon / 100
)
aibl <- data.frame (
  cpg = "cg05157625", 
  study = "AIBL",
  estimate = one$AIBL_AD_vs_CN_Estimate.bacon / 100, 
  se = one$AIBL_AD_vs_CN_StdErr.bacon / 100
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
  title = "CpG cg05157625"
)
pdf ("CpG cg05157625_blood.pdf",height = 2.5,width = 7)
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
dev.off()

#-------------------------------------------------------------------------------
# ROSMAP, Mt Sinai, London, Gasparoni 
#-------------------------------------------------------------------------------
res.blood <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/code_validation/Meta_analysis_code/NatComm_revision/DATASETS/meta_analysis_single_cpg_bacon_df.csv"
)

cpg <- res.blood %>% 
  dplyr::filter(cpg == "cg05157625") %>% 
  dplyr::select(c("cpg",contains("bacon")))

gasparoni <- data.frame (
  cpg = "cg05157625", 
  study = "GASPARONI", 
  estimate = cpg$GASPARONI_Estimate.bacon / 100.00, 
  se = cpg$GASPARONI_StdErr.bacon / 100.00
)
mtsinai <- data.frame (
  cpg = "cg05157625", 
  study = "MTSINAI",
  estimate = cpg$MTSINAI_Estimate.bacon / 100.00, 
  se = cpg$MTSINAI_StdErr.bacon / 100.00
)

london <- data.frame (
  cpg = "cg05157625", 
  study = "LONDON",
  estimate = cpg$LONDON_Estimate.bacon  / 100.00,
  se = cpg$LONDON_StdErr.bacon  / 100.00
)
rosmap <- data.frame (
  cpg = "cg05157625", 
  study = "ROSMAP",
  estimate = cpg$ROSMAP_Estimate.bacon  / 100.00, 
  se = cpg$ROSMAP_StdErr.bacon  / 100.00
)

brain <- rbind (rosmap, london, mtsinai, gasparoni)

meta_bin <- metagen (
  TE = estimate, 
  seTE = se, 
  comb.fixed = TRUE,
  comb.random = FALSE,
  prediction = FALSE,
  sm="SMD",
  data = brain, 
  studlab = study, 
  title = "CpG cg05157625"
)
pdf ("CpG cg05157625_brain.pdf",height = 3,width = 7)
forest (
  meta_bin, 
  comb.fixed = TRUE, 
  comb.random = FALSE, 
  print.tau2 = TRUE,
  sortvar = TE, # sorting variable
  digits.sd = 2,
  digits.se = 2,
  squaresize = 0.8, # size of the squares in the plot
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
  colgap.forest.left = unit(5,"cm")
)
dev.off()
