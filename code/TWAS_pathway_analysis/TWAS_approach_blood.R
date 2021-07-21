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

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted-V2.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

# - use ADNI dataset with matched DNAm and RNA data 
# Only CN and Dementia samples
# CN 181, Dementia 84
load("datasets/Aux/ADNI_matched_rna_dnam_residuals.rda")
#-----------------------------------------------------------------------------
# - compute PC1 for the 50 AD vs. CN CpGs + cpgs in 9 comp-b DMRs (prcomp function can do it) 
#-----------------------------------------------------------------------------
pca <- prcomp(residuals.matched.met[rownames(residuals.matched.met) %in% cpgs.all,] %>% t,center = TRUE,scale = TRUE)

PC1 <- pca$x[,"PC1"]

# - test association between each gene with PC1 as predictor, adjusting for age, sex, batch, cell types
# - use Limma or linear model
res <- plyr::adply(residuals.matched.exp,.margins = 1,.fun = function(gene){
  res <- lm(gene ~ PC1) %>% summary() %>% coef
  res["PC1",]
},.progress = "time")

# - pathway analysis using FGSEA, test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
# - use |t-stat| to rank
res$Gene_symbol <- MethReg:::map_ensg_to_symbol(res$X1)

res$abs_tValue <- abs(res$`t value`) ## LW added
res <- subset (res, !is.na(res$Gene_symbol))  ## LW added

rank <- res$abs_tValue
names(rank) <- res$Gene_symbol

library(fgsea)
library(msigdbr)

set.seed(42)

# - test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
c2.cp.biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="BIOCARTA")
c2.cp.kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="KEGG")
c2.cp.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="REACTOME")
c2.cp.wiki     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="WIKIPATHWAYS")

c2.cp <- rbind (c2.cp.biocarta, c2.cp.kegg, c2.cp.reactome, c2.cp.wiki )

c2.cp.list <- split(x = c2.cp$human_gene_symbol, f = c2.cp$gs_name)

c2.cp.fgsea <- fgsea(
  pathways = c2.cp.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 500,
  scoreType = "pos"
)

topPathways <- c2.cp.fgsea[padj < 0.05][order(pval), pathway]


c2.cp.fgsea.collapsedPathways <- collapsePathways(
  fgseaRes = c2.cp.fgsea[pathway %in% topPathways],
  pathways = c2.cp.list, 
  stats    = rank
)

c2.cp.tab <- plotGseaTable(
  pathways = c2.cp.list[topPathways], 
  stats = rank, 
  fgseaRes = c2.cp.fgsea, 
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
  scoreType = "pos"
)

topPathways <- c5.bp.fgsea[ES > 0 & padj <0.05][order(pval), pathway]

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
gene.info <- getGenes(res$X1, fields=fields)
res.with.info <- cbind(res,gene.info[,fields])

load("datasets/Aux/great_EPIC_array_annotation.rda")
pc.with.info <- pca$rotation[,1:3] %>% as_tibble(rownames = "CpG")
pc.with.info <- pc.with.info[order(-abs(pc.with.info$PC1)),]
great$CpG <- great$cpg
great$cpg <- NULL
pc.with.info <- left_join(pc.with.info,great)

pca.samples.df <- pca$x[,1:3] %>% tibble::as_tibble(rownames = "sample")
pca.samples.df <- dplyr::left_join(pca.samples.df, metadata.dnam[,c("DX"),drop = F] %>% tibble::as_tibble(rownames = "sample"))

writexl::write_xlsx(
  list(
    "Gene ~ PC1" = res.with.info,
    "PCA rotated data" = pca.samples.df,
    "PCA variable loadings" = pc.with.info,
    "c5.bp.fgsea" = c5.bp.fgsea[order(c5.bp.fgsea$pval),],
    "c5.bp.fgsea main pathways" = c5.bp.fgsea[c5.bp.fgsea$pathway %in% c5.bp.fgsea.collapsedPathways$mainPathways],
    "c2.cp.fgsea" = c2.cp.fgsea[order(c2.cp.fgsea$pval),]
  ),
  path = "code/pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx"
)

plot <- ggpubr::ggboxplot(
  pca.samples.df,
  x = "DX", 
  y = "PC1", 
  fill = "DX", 
  add = "jitter"
) + ggpubr::stat_compare_means(method = "wilcox.test")

ggplot2::ggsave(
  plot = plot, 
  filename = "code/pathway_analysis/PC1_boxplot_blood.pdf",
  width = 4,
  height = 4
)
