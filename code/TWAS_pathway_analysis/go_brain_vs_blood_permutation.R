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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Brain
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
load("code/pathway_analysis/brain_residuals.rda")
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
rank.brain <- rank
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Blood
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted-V2.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

cpgs.all <- c(
  cpgs.ad.cn
) %>% unique


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
rank.blood <- rank


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Gene sets
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(msigdbr)
# - test genesets in MSigDB C2:CP, C5:BP, C7:IMMUNESIGDB
c2.cp.biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="BIOCARTA")
c2.cp.kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="KEGG")
c2.cp.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="REACTOME")
c2.cp.wiki     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="WIKIPATHWAYS")

c2.cp <- rbind (c2.cp.biocarta, c2.cp.kegg, c2.cp.reactome, c2.cp.wiki )

c2.cp.list <- split(x = c2.cp$human_gene_symbol, f = c2.cp$gs_name)

c5.bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
c5.bp.list <- split(x = c5.bp$human_gene_symbol, f = c5.bp$gs_name)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Permutation
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(fgsea)

set.seed(42)

sample(names(rank.blood))  %>% head
sample(names(rank.blood))  %>% head

# permute ranked gene lists in brain, and in blood 
results.c5 <- plyr::adply(1:1000,.margins = 1,.fun = function(x){
  rank.blood.perm <- rank.blood
  names(rank.blood.perm) <- sample(names(rank.blood)) # shuffle genes names
  rank.brain.perm <- rank.brain
  names(rank.brain.perm) <- sample(names(rank.brain)) # shuffle genes names
  
  
  fgsea.blood.permu <- fgsea(
    pathways = c5.bp.list, 
    stats    = rank.blood.perm,
    minSize  = 15,
    maxSize  = 500,
    eps = 0,
    scoreType = "pos"
  )
  fgsea.blood.permu.sig <- fgsea.blood.permu[fgsea.blood.permu$padj < 0.05,]
  
  
  fgsea.brain.permu <- fgsea(
    pathways = c5.bp.list, 
    stats    = rank.brain.perm,
    minSize  = 15,
    maxSize  = 500,
    eps = 0,
    scoreType = "pos"
  )
  fgsea.brain.permu.sig <- fgsea.brain.permu[fgsea.brain.permu$padj < 0.05,]
  
  data.frame(
    nSigBlood = nrow(fgsea.blood.permu.sig),
    nSigBrain = nrow(fgsea.brain.permu.sig),
    nOverlap = length(intersect(fgsea.blood.permu.sig$pathway,fgsea.brain.permu.sig$pathway))
  )
  
},.progress = "time",.id = "rep")
writexl::write_xlsx(results.c5,path = "code/pathway_analysis/C5_1000_permutations.xlsx")



# permute ranked gene lists in brain, and in blood 
results.c2 <- plyr::adply(1:1000,.margins = 1,.fun = function(x){
  rank.blood.perm <- rank.blood
  names(rank.blood.perm) <- sample(names(rank.blood)) # shuffle genes names
  rank.brain.perm <- rank.brain
  names(rank.brain.perm) <- sample(names(rank.brain)) # shuffle genes names
  
  fgsea.blood.permu <- fgsea(
    pathways = c2.cp.list, 
    stats    = rank.blood.perm,
    minSize  = 15,
    maxSize  = 500,
    eps = 0,
    scoreType = "pos"
  )
  fgsea.blood.permu.sig <- fgsea.blood.permu[fgsea.blood.permu$padj < 0.05,]
  
  
  fgsea.brain.permu <- fgsea(
    pathways = c2.cp.list, 
    stats    = rank.brain.perm,
    minSize  = 15,
    maxSize  = 500,
    eps = 0,
    scoreType = "pos"
  )
  fgsea.brain.permu.sig <- fgsea.brain.permu[fgsea.brain.permu$padj < 0.05,]
  
  data.frame(
    nSigBlood = nrow(fgsea.blood.permu.sig),
    nSigBrain = nrow(fgsea.brain.permu.sig),
    nOverlap = length(intersect(fgsea.blood.permu.sig$pathway,fgsea.brain.permu.sig$pathway))
  )
  
},.progress = "time",.id = "rep")
writexl::write_xlsx(results.c2,path = "code/pathway_analysis/C2_1000_permutations.xlsx")


