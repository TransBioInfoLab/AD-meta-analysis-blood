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

# The leading edge analysis allows for the GSEA to determine which subsets 
# (referred to as the leading edge subset) of genes contributed the most to 
# the enrichment signal of a given gene set's leading edge or core enrichment 

#------------------------------------------------------------------------------
# C5
# ----------------------------------------------------------------------------
library(dplyr)
blood <- readxl::read_xlsx(
  path = "analysis_results//pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx",sheet = 4
)
blood <- blood %>% dplyr::filter(padj < 0.05)


brain <- readxl::read_xlsx(
  path = "analysis_results/pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx",sheet = 3
)
brain <- brain %>% dplyr::filter(padj < 0.05)

brain <- brain[brain$pathway %in% intersect(brain$pathway,blood$pathway),]
blood <- blood[blood$pathway %in% intersect(brain$pathway,blood$pathway),]


merged.c5 <- dplyr::full_join(
  brain %>% dplyr::rename_with(.fn = function(x) paste0("Brain_",x)) %>% dplyr::rename(pathway = Brain_pathway),
  blood %>% dplyr::rename_with(.fn = function(x) paste0("Blood_",x)) %>% dplyr::rename(pathway = Blood_pathway)
)

brain.genes.contributed.the.most.to.enrichment <- stringr::str_split(brain$leadingEdge,pattern = ",")
names(brain.genes.contributed.the.most.to.enrichment) <- paste0("Brain_C5_",brain$pathway)

blood.genes.contributed.the.most.to.enrichment <- stringr::str_split(blood$leadingEdge,pattern = ",")
names(blood.genes.contributed.the.most.to.enrichment) <- paste0("Blood_C5_",blood$pathway)


doParallel::registerDoParallel(cores = 4)
Jaccard_df <-  plyr::ldply(
  .data = combn(c(brain.genes.contributed.the.most.to.enrichment,blood.genes.contributed.the.most.to.enrichment), 2, simplify = FALSE), 
  .fun = function(x) {
    data.frame(
      "Pathway_1" = names(x),
      "Pathway_2" = rev(names(x)),
      "Jaccard_index" = rep(length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])),2)
    )
  },.progress = "time",.parallel = TRUE
) %>% rbind(
  data.frame( 
    "Pathway_1" = c(names(brain.genes.contributed.the.most.to.enrichment),names(blood.genes.contributed.the.most.to.enrichment)), 
    "Pathway_2" = c(names(brain.genes.contributed.the.most.to.enrichment),names(blood.genes.contributed.the.most.to.enrichment)),
    "Jaccard_index" = 1
  )
)

library(tidyr)
Jaccard_matrix <- pivot_wider(Jaccard_df, names_from = c(Pathway_2), values_from = Jaccard_index) %>% as.data.frame()
rownames(Jaccard_matrix) <- Jaccard_matrix$Pathway_1
Jaccard_matrix <- Jaccard_matrix[,Jaccard_matrix$Pathway_1]  %>% as.matrix()
Jaccard_matrix.C5 <- Jaccard_matrix

merged.c5$Jaccard_index_brain_blood <- Jaccard_df[match(paste0("Blood_C5_",merged.c5$pathway,"Brain_C5_",merged.c5$pathway),paste0(Jaccard_df$Pathway_1,Jaccard_df$Pathway_2)),"Jaccard_index"]

# table(is.na(Jaccard_matrix))
library(ComplexHeatmap)
library(dendextend)

rows.order <- vegan::vegdist(Jaccard_matrix,method = "euclidean") %>% hclust(method = "average")
col_dend <- as.dendrogram(rows.order)
col_dend <- color_branches(col_dend, h = 1.5)
split <- data.frame(cutree(rows.order, h = 1.5))
splitNum <- max(unique(split[,1]))
order <- rows.order$order

ht <- ComplexHeatmap::Heatmap(
  matrix = Jaccard_matrix,
  row_order = order,
  cluster_columns = col_dend,
  column_split = splitNum,
  cluster_rows  = col_dend,
  row_split = splitNum,
  column_order = order,
  width = unit(8,"cm"),
  height = unit(8,"cm"),
  name = "Jaccard index",
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4)
) 

column.order <- column_order(ht)
cluster.df.c5 <- plyr::ldply(column.order,.fun = function(cluster){
  data.frame("Sample" = colnames(Jaccard_matrix)[cluster])
},.id = "cluster")

pdf("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/pathway_analysis/plots/Jaccard_matrix_C5.pdf",width = 10,height = 10)
ht %>% draw( heatmap_legend_side = "left")
dev.off()


#------------------------------------------------------------------------------
# C2
# ----------------------------------------------------------------------------
library(dplyr)
blood <- readxl::read_xlsx(
  path = "analysis_results/pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx",sheet = 6
)
blood <- blood %>% dplyr::filter(padj < 0.05)


brain <- readxl::read_xlsx(
  path = "analysis_results/pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx",sheet = 5
)
brain <- brain %>% dplyr::filter(padj < 0.05)

brain <- brain[brain$pathway %in% intersect(brain$pathway,blood$pathway),]
blood <- blood[blood$pathway %in% intersect(brain$pathway,blood$pathway),]

merged.c2 <- dplyr::full_join(
  brain %>% dplyr::rename_with(.fn = function(x) paste0("Brain_",x)) %>% dplyr::rename(pathway = Brain_pathway),
  blood %>% dplyr::rename_with(.fn = function(x) paste0("Blood_",x)) %>% dplyr::rename(pathway = Blood_pathway)
)
brain.genes.contributed.the.most.to.enrichment <- stringr::str_split(brain$leadingEdge,pattern = ",")
names(brain.genes.contributed.the.most.to.enrichment) <- paste0("Brain_C2_",brain$pathway)

blood.genes.contributed.the.most.to.enrichment <- stringr::str_split(blood$leadingEdge,pattern = ",")
names(blood.genes.contributed.the.most.to.enrichment) <- paste0("Blood_C2_",blood$pathway)


doParallel::registerDoParallel(cores = 4)
Jaccard_df <-  plyr::ldply(
  .data = combn(c(brain.genes.contributed.the.most.to.enrichment,blood.genes.contributed.the.most.to.enrichment), 2, simplify = FALSE), 
  .fun = function(x) {
    data.frame(
      "Pathway_1" = names(x),
      "Pathway_2" = rev(names(x)),
      "Jaccard_index" = rep(length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])),2)
    )
  },.progress = "time",.parallel = TRUE
) %>% rbind(
  data.frame( 
    "Pathway_1" = c(names(brain.genes.contributed.the.most.to.enrichment),names(blood.genes.contributed.the.most.to.enrichment)), 
    "Pathway_2" = c(names(brain.genes.contributed.the.most.to.enrichment),names(blood.genes.contributed.the.most.to.enrichment)),
    "Jaccard_index" = 1
  )
)

library(tidyr)
Jaccard_matrix <- pivot_wider(Jaccard_df, names_from = c(Pathway_2), values_from = Jaccard_index) %>% as.data.frame()
rownames(Jaccard_matrix) <- Jaccard_matrix$Pathway_1
Jaccard_matrix <- Jaccard_matrix[,Jaccard_matrix$Pathway_1]  %>% as.matrix()
Jaccard_matrix.C2 <- Jaccard_matrix

# table(is.na(Jaccard_matrix))
library(ComplexHeatmap)
library(dendextend)

rows.order <- vegan::vegdist(Jaccard_matrix,method = "euclidean") %>% hclust(method = "average")
col_dend <- as.dendrogram(rows.order)
col_dend <- color_branches(col_dend, h = 1.5)
split <- data.frame(cutree(rows.order, h = 1.5))
splitNum <- max(unique(split[,1]))
order <- rows.order$order

ht <- ComplexHeatmap::Heatmap(
  matrix = Jaccard_matrix,
  row_order = order,
  cluster_columns = col_dend,
  column_split = splitNum,
  cluster_rows  = col_dend,
  row_split = splitNum,
  column_order = order,
  width = unit(8,"cm"),
  height = unit(8,"cm"),
  name = "Jaccard index",
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4)
) 

column.order <- column_order(ht)
names(column.order) <- paste0(1:length(column.order))
cluster.df.c2 <- plyr::ldply(column.order,.fun = function(cluster){
  data.frame("Sample" = colnames(Jaccard_matrix)[cluster])
},.id = "cluster")


pdf("analysis_results/pathway_analysis/plots/Jaccard_matrix_C2.pdf",width = 10,height = 10)
ht %>% draw( heatmap_legend_side = "left")
dev.off()

merged.c2$Jaccard_index_brain_blood <- Jaccard_df[match(paste0("Blood_C2_",merged.c2$pathway,"Brain_C2_",merged.c2$pathway),paste0(Jaccard_df$Pathway_1,Jaccard_df$Pathway_2)),"Jaccard_index"]

writexl::write_xlsx(
  list(
    "C2 - Brain and Blood" = merged.c2,
    "C5 - Brain and Blood" = merged.c5, 
    "C2 - Jaccard matrix" = Jaccard_matrix.C2 %>% as_tibble(rownames = "Pathways"),
    "C5 - Jaccard matrix"= Jaccard_matrix.C5 %>% as_tibble(rownames = "Pathways"),
    "C2 - Cluster Jaccard matrix" = cluster.df.c2,
    "C5 - Cluster Jaccard matrix" = cluster.df.c5
  ),
  path = "analysis_results/pathway_analysis/Jaccard_matrix.xlsx"
)
