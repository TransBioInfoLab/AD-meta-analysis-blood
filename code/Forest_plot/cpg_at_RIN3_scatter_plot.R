library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)
# forest plot for cg05157625 at RIN3 gene
cpg <- "cg05157625"
gene <- "ENSG00000100599"

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ADNI
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
load("datasets/Aux/ADNI_matched_rna_dnam_residuals.rda")

df.adni <- data.frame(
  "DNAm"  = residuals.matched.met[cpg,],
  "RNA" =  residuals.matched.exp[gene,],
  "AD_status" = metadata.dnam$DX
)

rlm <- MASS::rlm(
  RNA ~ DNAm + AD_status, 
  data = df.adni,
  psi = MASS::psi.bisquare,
  maxit = 100
) %>% summary %>% coef %>% data.frame

degrees.freedom.value <- nrow(df.adni) - 3
rlm$pval <- 2 * (1 - pt( abs(rlm$t.value), df = degrees.freedom.value) )

quant.pval <- rlm[-1,4,drop = FALSE] %>%
  t %>%
  as.data.frame()
colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")

quant.estimate <- rlm[-1,1,drop = FALSE] %>%
  t %>%
  as.data.frame()
colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")

plot.adni <- ggpubr::ggscatter(
  data = df.adni,
  x = "DNAm",
  title = "ADNI blood samples",
  y = "RNA",
  xlab = "DNAm residuals\ncg05157625",
  ylab = "RNA residuals\nRIN3 (ENSG00000100599)"
) + geom_smooth(
  method = MASS::rlm,
  method.args = list(psi = "psi.bisquare"),
  se = FALSE
) + annotate(
  "text", 
  x = -1.2, y = 3, 
  label = paste0(
    "Effect estimate = ", 
    quant.estimate[1] %>%  as.numeric() %>% formatC(digits = 2),
    "\nP-value = ",
    quant.pval[1] %>% as.numeric() %>% formatC(digits = 2)
  )
)


#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
load("~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/DNAm_RNA/data/matched_data.rda")
dim(matched.dnam)
dim(matched.exp)

# 1) remove confounding effects in DNAm data: 
resid_met <- coMethDMR:::GetResiduals(
  dnam = matched.dnam[rownames(matched.dnam) %in% "cg05157625",,drop = FALSE],
  betaToM = TRUE, #converts to M-values for fitting linear model 
  pheno_df = matched.phenotype,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","msex","age_death"), 
  nCores_int = 1,
  progressbar = TRUE  
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
matched.phenotype$rnaseq_id  <- map$rnaseq_id[match(matched.phenotype$Sample,map$mwas_id)]
residuals.matched.exp <- plyr::adply(
  .data = matched.exp.log2[gene,],
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
rownames(residuals.matched.exp) <- gene



df <- data.frame(
  "DNAm"  = resid_met %>% as.numeric(),
  "RNA" =  residuals.matched.exp %>% as.numeric(),
  "Braak_stage" = matched.phenotype$braaksc
)


rlm <- MASS::rlm(
  RNA ~ DNAm + Braak_stage, 
  data = df,
  psi = MASS::psi.bisquare,
  maxit = 100
) %>% summary %>% coef %>% data.frame

degrees.freedom.value <- nrow(df) - 3
rlm$pval <- 2 * (1 - pt( abs(rlm$t.value), df = degrees.freedom.value) )

quant.pval <- rlm[-1,4,drop = FALSE] %>%
  t %>%
  as.data.frame()
colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")

quant.estimate <- rlm[-1,1,drop = FALSE] %>%
  t %>%
  as.data.frame()
colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")


plot.rosmap <- ggpubr::ggscatter(
  data = df,
  title = "ROSMAP brain samples",
  x = "DNAm",
  y = "RNA",
  xlab = "DNAm residuals\ncg05157625",
  ylab = "RNA residuals\nRIN3 (ENSG00000100599)"
) + geom_smooth(
  method = MASS::rlm,
  method.args = list(psi = "psi.bisquare"),
  se = FALSE
) + annotate(
  "text", 
  x = -1.0, y = 1.1, 
  label = paste0(
    "Effect estimate = ", 
    quant.pval[1] %>%  as.numeric() %>% formatC(digits = 3),
    "\nP-value = ",
    quant.estimate[1] %>% as.numeric() %>% formatC(digits = 3)
  )
)

p <- grid.arrange(plot.rosmap,plot.adni,nrow = 1)
ggsave(p,filename = "plots/RIN3_cg05157625_scatter.pdf",width = 10,height = 5)
