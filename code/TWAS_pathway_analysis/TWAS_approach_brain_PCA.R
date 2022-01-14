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
# Date: 21 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article section: 
# Integrative analysis revealed gene expressions associated with DNA methylation 
# differences in the blood and the brain converges in biological pathways	
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Libs
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(readr)
library(dplyr)
library(fgsea)
library(msigdbr)
library(SummarizedExperiment)

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
cpgs.brain <- readxl::read_xlsx("datasets/brain_meta_analysis/Supplementary Data 1-24.xlsx",skip = 3)
dim(cpgs.brain)
cpgs.all <- cpgs.brain$cpg

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP DNAm data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
load("datasets/brain_meta_analysis/ROSMAP_matched_data.rda")
dim(matched.dnam)
dim(matched.exp)
matched.exp <- matched.exp[rowSums(matched.exp) > 0,]


# 1) remove confounding effects in DNAm data: 
resid_met <- coMethDMR:::GetResiduals(
  dnam = matched.dnam[rownames(matched.dnam) %in% cpgs.all,],
  betaToM = TRUE, #converts to Mvalues for fitting linear model 
  pheno_df = matched.phenotype,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","msex","age_death"), 
  nCores_int = 6,
  progressbar = TRUE  
)

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

doParallel::registerDoParallel(cores = 6)
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
  .parallel = TRUE)
rownames(residuals.matched.exp) <- rownames(matched.exp.log2)

save(
  residuals.matched.exp,
  matched.exp.log2,
  resid_met,cpgs.all,
  file = "code/TWAS_pathway_analysis/brain_residuals.rda"
)

#-----------------------------------------------------------------------------
# - compute PC1 for the 3751 cpgs
#-----------------------------------------------------------------------------
load("code/TWAS_pathway_analysis/brain_residuals.rda")
pca <- prcomp(resid_met[rownames(resid_met) %in% cpgs.all,] %>% t,center = TRUE,scale = TRUE)
PC1 <- pca$x[,"PC1"]

# - test association between each gene with PC1 as predictor, adjusting for age, sex, batch, cell types
# - use Limma or linear model
res <- plyr::adply(residuals.matched.exp %>% as.matrix ,.margins = 1,.fun = function(gene){
  res <- lm(gene ~ PC1) %>% summary() %>% coef
  res["PC1",]
},.progress = "time")

# - pathway analysis using FGSEA, test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
# - use |t-stat| to rank
res$Gene_symbol <- MethReg:::map_ensg_to_symbol(res$X1)
res$abs_tValue <- abs(res$`t value`) ## LW added
res <- subset (res, !is.na(res$Gene_symbol))  ## LW added
res <- res[!duplicated(res$Gene_symbol),]
rank <- res$abs_tValue
names(rank) <- res$Gene_symbol

library(fgsea)
library(msigdbr)

set.seed(42)

# - test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
c2.cp.biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
c2.cp.kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
c2.cp.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
c2.cp.wiki     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "WIKIPATHWAYS")
c2.cp.cgp    <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
c2.cp <- rbind (c2.cp.biocarta, c2.cp.kegg, c2.cp.reactome, c2.cp.wiki, c2.cp.cgp)

c2.cp.list <- split(x = c2.cp$human_gene_symbol, f = c2.cp$gs_name)

c2.cp.fgsea <- fgsea(
  pathways = c2.cp.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 700,
  eps = 0,
  scoreType = "pos"
)

topPathways <- c2.cp.fgsea[padj < 0.05][order(pval), pathway]
c2.cp.tab <- plotGseaTable(
  c2.cp.list[topPathways], 
  rank, 
  c2.cp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)


c2.cp.fgsea.collapsedPathways <- collapsePathways(
  fgseaRes = c2.cp.fgsea[pathway %in% topPathways],
  pathways = c2.cp.list, 
  stats    = rank
)
c2.cp.tab.collapsedPathways <- plotGseaTable(
  c2.cp.list[c2.cp.fgsea.collapsedPathways$mainPathways], 
  rank, 
  c2.cp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)


c5.bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
c5.bp.list <- split(x = c5.bp$human_gene_symbol, f = c5.bp$gs_name)

c5.bp.fgsea <- fgsea(
  pathways = c5.bp.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 500,
  eps = 0,
  scoreType = "pos"
)


topPathways <- c5.bp.fgsea[padj < 0.05][order(pval), pathway]
c5.bp.tab <- plotGseaTable(
  c5.bp.list[topPathways], 
  rank, 
  c5.bp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)

c5.bp.fgsea.collapsedPathways <- collapsePathways(
  fgseaRes = c5.bp.fgsea[pathway %in% topPathways],
  pathways = c5.bp.list, 
  stats    = rank
)

c5.bp.tab.collapsedPathways <- plotGseaTable(
  c5.bp.list[c5.bp.fgsea.collapsedPathways$mainPathways], 
  rank, 
  c5.bp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)

plots <- gridExtra::marrangeGrob(
  list(
    "C5 BP" = c5.bp.tab,
    "C2 CP" = c2.cp.tab,
    "C2 CP collapsedPathways" = c2.cp.tab.collapsedPathways,
    "C5 BP collapsedPathways" = c5.bp.tab.collapsedPathways
  ),
  top = c(""),
  ncol = 1,nrow = 1
)

ggplot2::ggsave(
  plot = plots,
  filename = "plots/pathway_analysis_PC_analysis_rank_abs_t_stast.pdf",
  width = 20,
  height = 8
)

pca.df <-  pca$x %>% as.data.frame
pca.df <- cbind(rownames(pca.df),pca.df)
colnames(pca.df)[1] <- "Samples"

c5.bp.fgsea$leadingEdge <- lapply(c5.bp.fgsea$leadingEdge,FUN = function(x)paste(x,collapse = ",")) %>% unlist
c2.cp.fgsea$leadingEdge <- lapply(c2.cp.fgsea$leadingEdge,FUN = function(x)paste(x,collapse = ","))  %>% unlist


library(mygene)

fields <- c(
  "alias","ensembl.type_of_gene","name","genomic_pos.chr",
  "genomic_pos.start","genomic_pos.end","genomic_pos.strand","summary"
)
gene.info <- getGenes(res$X1, fields = fields)
gene.info$X1 <- gene.info$query
gene.info[,1:3] <- NULL
res.with.info <- dplyr::left_join(res,gene.info %>% as.data.frame())

load("datasets/Aux/great_EPIC_array_annotation.rda")
pc.with.info <- pca$rotation[,1:3] %>% tibble::as_tibble(rownames = "cpg")
pc.with.info <- pc.with.info[order(-abs(pc.with.info$PC1)),]
pc.with.info <- dplyr::left_join(pc.with.info,great)

pc.samples <- pca$x[,1:3] %>% tibble::as_tibble(rownames = "Samples")
pc.samples$braaksc <-  matched.phenotype$braaksc
pc.samples$braaksc_merged <-  matched.phenotype$stage3  %>% as.factor
writexl::write_xlsx(
  list(
    "Gene ~ PC1" = res.with.info,
    "PCA variable loadings" = pc.with.info,
    "c5.bp.fgsea" = c5.bp.fgsea,
    "c5.bp.fgsea main pathways" = c5.bp.fgsea[c5.bp.fgsea$pathway %in% c5.bp.fgsea.collapsedPathways$mainPathways],
    "c2.cp.fgsea" = c2.cp.fgsea,
    "c2.cp.fgsea main pathways" = c2.cp.fgsea[c2.cp.fgsea$pathway %in% c2.cp.fgsea.collapsedPathways$mainPathways],    
    "PCA variable loadings - Samples" = pc.samples
  ),
  path = "analysis_results/pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx"
)


plot <- ggpubr::ggboxplot(pc.samples,x = "braaksc", y ="PC1", fill = "braaksc", add = "jitter",title = "Pearson correlation PC1 vs braaksc: r = -0.15, p-value = 0.0003387")

ggplot2::ggsave(
  plot = plot, 
  filename = "analysis_results/pathway_analysis/PC1_boxplot_brain.pdf",
  width = 4,
  height = 4
)

plot.merged <- ggpubr::ggboxplot(pc.samples,x = "braaksc_merged", y ="PC1", fill = "braaksc_merged", add = "jitter",title = "Pearson correlation PC1 vs braaksc: r = -0.15, p-value = 0.0003387") 

ggplot2::ggsave(
  plot = plot.merged, 
  filename = "analysis_results/pathway_analysis/PC1_boxplot_brain_merged.pdf",
  width = 4,
  height = 4
)
stats::cor.test(pc.samples$PC1,pc.samples$braaksc)
