#-----------------------------------------------------------------------------
# Libs
#-----------------------------------------------------------------------------
library(readr)
library(dplyr)
library(matrixStats)
library(SummarizedExperiment)
library(MethReg)
library(DelayedMatrixStats)
library(GenomicRanges)
library(doParallel)


########################    
# TWAS sig genes
######################## 
TWAS_brain <- readxl::read_xlsx( "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx")
TWAS_brain.sig <- TWAS_brain[TWAS_brain$`Pr(>|t|)` < 0.05,]

TWAS_blood <- readxl::read_xlsx( "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx")
TWAS_blood.sig <- TWAS_blood[TWAS_blood$`Pr(>|t|)` < 0.05,]

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
rosmap.data.rnaseq.path <- "~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/DNAm_RNA/data/"
rosmap.data.dnam.path <- "~/TBL Dropbox/Tiago Silva/coMethDMR_metaAnalysis/code_validation/Meta_analysis_code/DATASETS/ROSMAP/"
#------------------------------------------------------------------------------
# 1) RNA-SEQ
# Read and add gene name
#------------------------------------------------------------------------------

# retrieve gene information to map from ID to Symbol
gene.info <- TCGAbiolinks::get.GRCh.bioMart("hg19")

rna.seq <-  file.path(rosmap.data.rnaseq.path,"ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv") %>% 
  readr::read_tsv() %>%
  dplyr::full_join(
    readr::read_tsv(file.path(rosmap.data.rnaseq.path,"ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv")) 
  )
# rna.seq$geneName <-  gene.info$external_gene_name[match(gsub("\\.[0-9]*", "", rna.seq$gene_id),  gene.info$ensembl_gene_id)]



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Map RNA-seq and DNA methylation
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map <- readr::read_csv(file.path(rosmap.data.rnaseq.path,"ROSMAP_IDkey.csv"))
# Only keep samples with DNA methytlation and gene expression
map <- na.omit(unique(map[, c("projid","mwas_id", "rnaseq_id")]))
map <- map[map$rnaseq_id %in% gsub("_[0-9]$", "", c(colnames(rna.seq))),]

rna.seq <- rna.seq %>% as.data.frame()
rownames(rna.seq) <- gsub("\\.[0-9]*$", "", rna.seq$gene_id)

rna.seq <- rna.seq[,gsub("_[0-9]$", "", colnames(rna.seq)) %in% map$rnaseq_id]
map <- map[match(gsub("_[0-9]$", "", colnames(rna.seq)),map$rnaseq_id ),]


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Clinical
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
clinical <- readr::read_csv(file.path(rosmap.data.rnaseq.path,"ROSMAP_clinical.csv"), col_types = readr::cols())
clinical <- dplyr::right_join(clinical,map)
clinical <- clinical[match(gsub("_[0-9]$", "", colnames(rna.seq)),clinical$rnaseq_id ),] %>% as.data.frame()

#-----------------  CHECKS   ----------------------
nrow(clinical)  == ncol(rna.seq)
#-------------------------------------------

#------------------------------------------------------------------------------
# Cell type markers
#------------------------------------------------------------------------------
genes <- c("ENO2", "OLIG2", "CD34", "CD68", "GFAP")
gene.ensg <- gene.info$ensembl_gene_id[match(genes, gene.info$external_gene_name)]

cell.type.markers  <- log2(rna.seq + 1)[match(gene.ensg, rownames(rna.seq)),] %>% t
colnames(cell.type.markers) <- genes

cell.type.markers <- cell.type.markers %>%  tibble::as_tibble(rownames = "rnaseq_id")

cell.type.markers$Sample <- map$mwas_id[match(gsub("_[0-9]$", "", cell.type.markers$rnaseq_id), map$rnaseq_id)]

clinical$ENO2 <- cell.type.markers$ENO2[match(clinical$rnaseq_id, gsub("_[0-9]$", "", cell.type.markers$rnaseq_id))]
clinical$OLIG2 <- cell.type.markers$OLIG2[match(clinical$rnaseq_id, gsub("_[0-9]$", "", cell.type.markers$rnaseq_id))]
clinical$CD34 <- cell.type.markers$CD34[match(clinical$rnaseq_id, gsub("_[0-9]$", "", cell.type.markers$rnaseq_id))]
clinical$CD68 <- cell.type.markers$CD68[match(clinical$rnaseq_id, gsub("_[0-9]$", "", cell.type.markers$rnaseq_id))]
clinical$GFAP <- cell.type.markers$GFAP[match(clinical$rnaseq_id, gsub("_[0-9]$", "", cell.type.markers$rnaseq_id))]

##### Modeling rnaseq ~ AD stage to find DE genes #############################
library(dplyr)
library(TCGAbiolinks)

# Log2(counts + 1) ~ braaksc + age_death + msex + Log2(ENO2 + 1) + Log2(OLIG2 + 1) + Log2(CD34 + 1)  + Log2(CD68 + 1)  +  Log2(GFAP + 1) 
doParallel::registerDoParallel(cores = 4)

matrix <-  log2(rna.seq + 1)
matrix.twas.genes <- matrix[rownames(matrix) %in% base::union(TWAS_brain.sig$X1,TWAS_blood.sig$X1),]

tab <- plyr::adply(
  .data = matrix.twas.genes,
  .margins = 1,
  .fun = function(row) {
    
    tryCatch({
      val <- t(row)
      colnames(val) <- "val"
      dat <- cbind(val, clinical)
     
      results.all <- lm(
        val ~ braaksc + age_death + msex + ENO2 + OLIG2 + CD34 + CD68 + GFAP,
        data = dat
      )
      results.all.estimate <- summary(results.all)$coefficients["braaksc", "Estimate"]
      results.all.pval <- summary(results.all)$coefficients["braaksc", "Pr(>|t|)"]
      
      data.frame(
        ensembl_gene_id = row.names(row),
        stage_estimate = results.all.estimate,
        stage_pVal = results.all.pval,
        row.names = FALSE,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      print(row)
      return()
    })
    
  },
  .parallel = FALSE,
  .progress = "time"
)

tab$stage_fdr <- p.adjust(tab$stage_pVal, method = "fdr")

tab$TWAS_sig[tab$ensembl_gene_id %in% TWAS_blood.sig$X1] <- "Sig in Blood"
tab$TWAS_sig[tab$ensembl_gene_id %in% TWAS_brain.sig$X1] <- "Sig in Brain"
tab$TWAS_sig[tab$ensembl_gene_id %in% intersect(TWAS_brain.sig$X1,TWAS_blood.sig$X1)] <- "Sig in Brain and Blood"

aux <- plyr::rbind.fill(TWAS_blood.sig,TWAS_brain.sig)
table(is.na(match(tab$ensembl_gene_id,aux$X1)))
tab.annotated <- cbind(
  tab,
  aux[match(tab$ensembl_gene_id,aux$X1),
      c("Gene_symbol","alias","ensembl.type_of_gene","name","genomic_pos.chr","genomic_pos.start",
        "genomic_pos.end", "genomic_pos.strand","summary")
  ])


tab_final <- tab.annotated[order(tab.annotated$stage_pVal),]
tab_final.sig <- tab_final[tab_final$stage_fdr < 0.05,]


writexl::write_xlsx(
  list("Genes FDR < 0.05" = tab_final.sig,"All genes" = tab_final), 
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/pathway_analysis/TWAS_geens_diff_expreessed_analysis_brain.xlsx"
)

