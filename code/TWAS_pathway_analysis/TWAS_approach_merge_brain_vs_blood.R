library(dplyr)

#-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=
#  GO
#-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=

# Sheets
# 1 - "Gene ~ PC1" 
# 2 - "PCA rotated data" 
# 3 - "PCA variable loadings" 
# 4 - "c5.bp.fgsea" 
# 5 - "c5.bp.fgsea main pathways" 
# 6 - "c2.cp.fgsea" 
blood <- readxl::read_xlsx(
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx",sheet = 4
)
blood <- blood %>% dplyr::filter(padj < 0.05)
blood <- blood[,c("pathway","pval","padj","NES")]
colnames(blood)[-1] <- paste0("blood_",colnames(blood)[-1])

# Sheets
# 1 - "Gene ~ PC1" 
# 2 - "PCA variable loadings" 
# 3 - "c5.bp.fgsea" = c5.bp.fgsea,
# 4 - "c5.bp.fgsea main pathways" 
# 5 - "c2.cp.fgsea" = c2.cp.fgsea,
# 6 - "c2.cp.fgsea main pathways" 
# 7 - "PCA variable loadings - Samples" 
brain <- readxl::read_xlsx(
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx",sheet = 3
)
brain <- brain %>% dplyr::filter(padj < 0.05)
brain <- brain[,c("pathway","pval","padj","NES")]
colnames(brain)[-1] <- paste0("brain_",colnames(brain)[-1])

library(ggVennDiagram)
library(ggplot2)

pdf("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/plots/GO_Brain_vs_Blood.pdf")
ggVennDiagram(
  list("Blood pathways(P.Adj < 0.05)" = blood$pathway, 
       "Brain pathways (P.Adj < 0.05)" = brain$pathway), 
  label_alpha = 0,
  label_color = "white"
) + theme(legend.position = "none")
dev.off()

merged <- dplyr::full_join(brain,blood)

writexl::write_xlsx(
  list(
    "Brain and Blood" = merged,
    "overlap" = data.frame(
      "Pathways" = unique(c(blood$pathway,brain$pathway)),
      "Brain FDR sig" = unique(c(blood$pathway,brain$pathway)) %in% brain$pathway, 
      "Blood FDR sig" = unique(c(blood$pathway,brain$pathway)) %in% blood$pathway
    )
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/Brain_blood_c5.bp.fgsea.xlsx"
)


#-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=
# 
#-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=--=

library(dplyr)
blood <- readxl::read_xlsx(
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx",sheet = 6
)
blood <- blood %>% dplyr::filter(padj < 0.05)
colnames(blood)[-1] <- paste0("blood_",colnames(blood)[-1])


brain <- readxl::read_xlsx(
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx",sheet = 5
)
brain <- brain %>% dplyr::filter(padj < 0.05)
colnames(brain)[-1] <- paste0("brain_",colnames(brain)[-1])

library(ggVennDiagram)
library(ggplot2)

pdf("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/plots/GO_Brain_vs_Blood.pdf")
ggVennDiagram(
  list("Blood pathways(P.Adj < 0.05)" = blood$pathway, 
       "Brain pathways (P.Adj < 0.05)" = brain$pathway), 
  label_alpha = 0,
  label_color = "white"
) + theme(legend.position = "none")
dev.off()

merged <- dplyr::full_join(brain,blood)

writexl::write_xlsx(
  list(
    "Brain and Blood" = merged,
    "overlap" = data.frame(
      "Pathways" = unique(c(blood$pathway,brain$pathway)),
      "Brain FDR sig" = unique(c(blood$pathway,brain$pathway)) %in% brain$pathway, 
      "Blood FDR sig" = unique(c(blood$pathway,brain$pathway)) %in% blood$pathway
    )
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/Brain_blood_c2.cp.fgsea.xlsx"
)


merged <- dplyr::full_join(brain,blood)

writexl::write_xlsx(
  list(
    "Brain leading Edge genes" = strsplit(brain[brain$pathway == "REACTOME_NEUTROPHIL_DEGRANULATION",]$brain_leadingEdge,",") %>% 
      unlist %>% as.data.frame %>% dplyr::rename( "Brain leading Edge genes REACTOME_NEUTROPHIL_DEGRANULATION"= "."),
    
    "Blood leading Edge genes" = strsplit(blood[blood$pathway == "REACTOME_NEUTROPHIL_DEGRANULATION",]$blood_leadingEdge,",") %>% 
      unlist %>% as.data.frame %>% dplyr::rename( "Blood leading Edge genes REACTOME_NEUTROPHIL_DEGRANULATION"= ".")
  ),
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/REACTOME_NEUTROPHIL_DEGRANULATION_Brain_vs_Blood.xlsx"
)




pdf("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/plots/REACTOME_NEUTROPHIL_DEGRANULATION_Brain_vs_Blood.pdf")
ggVennDiagram(
  list("Brain leading Edge genes" = strsplit(brain[brain$pathway == "REACTOME_NEUTROPHIL_DEGRANULATION",]$brain_leadingEdge,",") %>% unlist, 
       "Blood leading Edge genes" = strsplit(blood[blood$pathway == "REACTOME_NEUTROPHIL_DEGRANULATION",]$blood_leadingEdge,",") %>% unlist
), 
  label_alpha = 0,
  label_color = "white",
) + ggtitle("REACTOME_NEUTROPHIL_DEGRANULATION")  + 
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
dev.off()