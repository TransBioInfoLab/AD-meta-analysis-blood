#-----------------------------------------------------------------------------
# Libs
#-----------------------------------------------------------------------------
library(readr)
library(dplyr)
library(matrixStats)
library(SummarizedExperiment)
library(MethReg)

########################    
# TWAS sig genes
######################## 
TWAS_brain <- readxl::read_xlsx( "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_brain.xlsx")
TWAS_brain.sig <- TWAS_brain[TWAS_brain$`Pr(>|t|)` < 0.05,]

TWAS_blood <- readxl::read_xlsx( "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results//pathway_analysis/TWAS_analysis_only_cpg_blood.xlsx")
TWAS_blood.sig <- TWAS_blood[TWAS_blood$`Pr(>|t|)` < 0.05,]

#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
cohort <- "ADNI"
dir.base <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/"
dir.data <- file.path(dir.base,"datasets/",cohort,"/data/DNA_methylation") 
dir.data.clinical <- file.path(dir.base,"datasets/",cohort,"/data/Clinical") 

#-----------------------------------------------------------------------------
# GET ADNI clinical data
#-----------------------------------------------------------------------------
clinical <- readr::read_csv(file.path(dir.data.clinical,"ADNIMERGE_downloaded_3-1-2021.csv"))
clinical$ID <- paste0(clinical$RID,clinical$COLPROT,clinical$EXAMDATE)
clinical$ID2 <- paste0(clinical$RID,clinical$COLPROT,clinical$VISCODE)

demo <- readr::read_csv(file.path(dir.data.clinical,"PTDEMOG_3-1-2021.csv"))
clinical$birth_month <- demo$PTDOBMM[match(clinical$RID, demo$RID)]
clinical$birth_year <- demo$PTDOBYY[match(clinical$RID, demo$RID)]
# we don't have the day so we set to one

library(lubridate)
clinical$age_at_visit <- interval(
  as.Date(paste0(clinical$birth_month,"/1/",clinical$birth_year), "%m/%d/%Y"),
  clinical$EXAMDATE
) %>% time_length(unit = "years")

clinical$RID_COLPROT_VISCODE <- paste0(clinical$RID,"_",clinical$COLPROT,"_",clinical$VISCODE)
clinical <- clinical[clinical$DX %in% c("CN","Dementia"),]

#-----------------------------------------------------------------------------
# GET ADNI expression data
#-----------------------------------------------------------------------------
# Gene expression: Affymetrix Human Genome U 219 array 
ADNI_Gene_Expression_Profile <- read_csv(
  file.path(dir.base,"datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),
  skip = 8
)
ADNI_Gene_Expression_Profile$X748 <- NULL

# We have our gene expression as below:
#   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
# 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
# 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
# 3 11715106_x_at NA        CTAGE6 || CTAGE15     2.43  2.27  2.33  2.26  2.33  2.45  2.17
# The row (3) must be break in to two
#   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
# 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
# 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
# 3 11715106_x_at NA        CTAGE6                2.43  2.27  2.33  2.26  2.33  2.45  2.17
# 4 11715106_x_at NA        CTAGE15               2.43  2.27  2.33  2.26  2.33  2.45  2.17
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile %>% tidyr::separate_rows("Symbol")

library(hgu219.db)
x <- hgu219ENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID, this is required by MethReg
probe.to.ensg <- as.data.frame(as.list(x) %>% unlist) %>% na.omit
ADNI_Gene_Expression_Profile$ENGS <- probe.to.ensg[ADNI_Gene_Expression_Profile$ProbeSet,]
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[!is.na(ADNI_Gene_Expression_Profile$ENGS),]

nrow(ADNI_Gene_Expression_Profile) # 43893

# dropping genes in bottom 10 percentile for over 80% of the samples
genes.low.expressed.samples.count <- plyr::aaply(
  ADNI_Gene_Expression_Profile[,grep("^X",colnames(ADNI_Gene_Expression_Profile))] ,
  .margins = 2, # for each sample set get all genes
  .fun = function(sample.genes){
    # for each sample, mark the bottom 10% less expressed genes as 1
    sample.genes[[1]]  <= quantile(sample.genes[[1]] , probs = c(.10),type = 3)
  }) %>% colSums()
genes.idx <- which(unname(genes.low.expressed.samples.count) > (length(grep("^X",colnames(ADNI_Gene_Expression_Profile))) * 0.8))
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[-c(genes.idx),]

nrow(ADNI_Gene_Expression_Profile) # 41681

# Since we have multiple probes mapping to the same gene in the array 
# we will take the median of the same genes
# as suggest in https://www.nature.com/articles/s41598-020-60595-1
expression.matrix <- plyr::aaply(unique(ADNI_Gene_Expression_Profile$ENGS),.margins = 1,.fun = function(gene){
  dat <- ADNI_Gene_Expression_Profile[ADNI_Gene_Expression_Profile$ENGS == gene,grep("^X",colnames(ADNI_Gene_Expression_Profile))]
  colMedians(dat %>% as.matrix)
},.progress = "time") 
rownames(expression.matrix) <- unique(ADNI_Gene_Expression_Profile$ENGS)

ADNI_Gene_Expression_Metadata <- read_csv(
  file.path(dir.base, "datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),
  skip = 0,
  col_names = FALSE,
  n_max = 7
)
ADNI_Gene_Expression_Metadata$X748 <- ADNI_Gene_Expression_Metadata$X3 <- ADNI_Gene_Expression_Metadata$X2 <- NULL

gene.exp.IDs <- apply(
  ADNI_Gene_Expression_Metadata[,2:ncol(ADNI_Gene_Expression_Metadata)],MARGIN = 2, 
  FUN = function(col) {
    paste0(
      stringr::str_extract(pattern = "[0-9]*$",string = col[3]),
      "_", col[1],"_",
      col[2]
    )
  }
)

# The expression data has the visit code in the old form
DXSUM <- readr::read_csv("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/datasets/ADNI/data/Clinical/DXSUM_PDXCONV_ADNIALL_downloaded_2-19-2021.csv")
DXSUM <- DXSUM[match(gene.exp.IDs,paste0(stringr::str_extract(pattern = "[0-9]*$",DXSUM$PTID),"_",DXSUM$Phase,"_",DXSUM$VISCODE)),]
gene.exp.IDs <- paste0(DXSUM$RID,"_",DXSUM$Phase,"_",DXSUM$VISCODE2)
colnames(expression.matrix) <- gene.exp.IDs
expression.matrix <- expression.matrix[,colnames(expression.matrix) != "NA_NA_NA"]

#-----------------------------------------------------------------------------
# Covariates information
#-----------------------------------------------------------------------------
ad.cn.idx <- which(colnames(expression.matrix) %in% clinical$RID_COLPROT_VISCODE)
expression.matrix.AD.CN <- expression.matrix[,ad.cn.idx]
ADNI_Gene_Expression_Metadata <- ADNI_Gene_Expression_Metadata %>% as.data.frame()
rownames(ADNI_Gene_Expression_Metadata) <- ADNI_Gene_Expression_Metadata$X1
ADNI_Gene_Expression_Metadata$X1 <- NULL
ADNI_Gene_Expression_Metadata.AD.CN <- ADNI_Gene_Expression_Metadata[,ad.cn.idx]
clinical.AD.CN <- clinical[match(colnames(expression.matrix.AD.CN),clinical$RID_COLPROT_VISCODE),]

Affy_Plate <- ADNI_Gene_Expression_Metadata.AD.CN["Affy Plate",] %>% as.numeric()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.AD.CN)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.AD.CN["Affy Plate",] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata.AD.CN["RIN",] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.AD.CN)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.AD.CN["RIN",] %>% as.character()))
]


metadata <- clinical.AD.CN[,c("age_at_visit","PTGENDER","DX")]
metadata$Affy_Plate <- Affy_Plate
metadata$RIN <- RIN

#library(xCell)
devtools::load_all("~/Downloads/xCell-master/")
aux <- expression.matrix.AD.CN
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.with.xcel <- cbind(metadata,t(xcell))
colnames(metadata.with.xcel)[6] <- "B_cells"
colnames(metadata.with.xcel)[8] <- "CD4_T_cells"
colnames(metadata.with.xcel)[9] <- "CD8_T_cells"



library(dplyr)
library(TCGAbiolinks)

matrix <- log2(expression.matrix.AD.CN)
matrix.twas.genes <- matrix[rownames(matrix) %in% base::union(TWAS_brain.sig$X1,TWAS_blood.sig$X1),]
doParallel::registerDoParallel(cores = 4)
tab <- plyr::adply(
  .data = matrix.twas.genes,
  .margins = 1,
  .fun = function(row) {
    
    tryCatch({
      val <- as.data.frame(row)
      colnames(val) <- "val"
      dat <- cbind(val, metadata.with.xcel)
      
      results.all <- lm(
        val ~ DX + age_at_visit + PTGENDER + B_cells + CD4_T_cells + CD8_T_cells + NKT + Monocytes + Neutrophils,
        data = dat
      )
      results.all.estimate <- summary(results.all)$coefficients["DXDementia", "Estimate"]
      results.all.pval <- summary(results.all)$coefficients["DXDementia", "Pr(>|t|)"]
      
      data.frame(
        DX_estimate = results.all.estimate,
        DX_pVal = results.all.pval
      )
      
    }, error = function(e) {
      print(e)
      return()
    })
    
  },
  .parallel = FALSE,
  .progress = "time",
  .id = "ensembl_gene_id"
)
tab <- tab %>% dplyr::rename(ensembl_gene_id = X1)
tab$DX_fdr <- p.adjust(tab$DX_pVal, method = "fdr")

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



tab_final <- tab.annotated[order(tab.annotated$DX_pVal),]
tab_final.sig <- tab_final[tab_final$DX_fdr < 0.05,]

writexl::write_xlsx(
  list("Genes FDR < 0.05" = tab_final.sig,"All genes" = tab_final), 
  path = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/analysis_results/pathway_analysis/TWAS_genes_diff_expreessed_analysis_blood.xlsx"
)


