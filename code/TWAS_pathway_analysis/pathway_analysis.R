library(fgsea)
library(msigdbr)

set.seed(42)

ranked.genes <- read.csv("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/LW_GSEA-pathway-analysis/blood/blood-cpgs-dmrs_neglog10Pvals.csv")
rank <- ranked.genes$resid_neglog10pval
names(rank) <- ranked.genes$Gene

# - test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
c2.cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
c2.bp.list <- split(x = c2.cp$gene_symbol, f = c2.cp$gs_name)

c2.cp.fgsea <- fgsea(
  pathways = c2.bp.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 500
)



topPathwaysUp <- c2.cp.fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- c2.cp.fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
c2.cp.tab <- plotGseaTable(
  c2.bp.list[topPathways], 
  rank, 
  c2.cp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)



c5.bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
c5.bp.list <- split(x = c5.bp$gene_symbol, f = c5.bp$gs_name)

c5.bp.fgsea <- fgsea(
  pathways = c5.bp.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 500
)


topPathwaysUp <- c5.bp.fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- c5.bp.fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
c5.bp.tab <- plotGseaTable(
  c5.bp.list[topPathways], 
  rank, 
  c5.bp.fgsea, 
  gseaParam = 0.5,
  render = FALSE
)


c7.IMMUNESIGDB <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB")
c7.IMMUNESIGDB.list <- split(x = c7.IMMUNESIGDB$gene_symbol, f = c7.IMMUNESIGDB$gs_name)

c7.IMMUNESIGDBfgsea <- fgsea(
  pathways = c7.IMMUNESIGDB.list, 
  stats    = rank,
  minSize  = 15,
  maxSize  = 500
)

topPathwaysUp <- c7.IMMUNESIGDBfgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- c7.IMMUNESIGDBfgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
c7.IMMUNESIGDB.tab <- plotGseaTable(
  c7.IMMUNESIGDB.list[topPathways], 
  rank, 
  c7.IMMUNESIGDBfgsea, 
  gseaParam = 0.5,
  render = FALSE
)

plots <- gridExtra::marrangeGrob(
  list(
    "c7.IMMUNESIGDB" = c7.IMMUNESIGDB.tab,
    "C5 BP" = c5.bp.tab,
    "C2 CP" = c2.cp.tab
  ),
  top = c(""),
  ncol = 1,nrow = 1
)
ggplot2::ggsave(plot = plots,filename = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/plots/pathway_analysis.pdf",width = 20,height = 8)

c5.bp.fgsea$leadingEdge <- lapply(c5.bp.fgsea$leadingEdge,FUN = function(x)paste(x,collapse = ",")) %>% unlist
c2.cp.fgsea$leadingEdge <- lapply(c2.cp.fgsea$leadingEdge,FUN = function(x)paste(x,collapse = ","))  %>% unlist
c7.IMMUNESIGDBfgsea$leadingEdge <- lapply(c7.IMMUNESIGDBfgsea$leadingEdge,FUN = function(x)paste(x,collapse = ","))  %>% unlist


writexl::write_xlsx(
  list(
    "c5.bp.fgsea" = c5.bp.fgsea,
    "c2.cp.fgsea" = c2.cp.fgsea,
    "c7.IMMUNESIGDBfgsea" = c7.IMMUNESIGDBfgsea
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/code/pathway_analysis/pathway_analysis.xlsx"
)