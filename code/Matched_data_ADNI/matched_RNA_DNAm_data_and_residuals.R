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
# Date: 12 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
library(matrixStats)
library(SummarizedExperiment)
library(MethReg)
library(xCell)
library(dorothea)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# GET ADNI data
#-----------------------------------------------------------------------------
# DNA methylation
cohort <- "ADNI"
dir.base <- "./"
dir.data <- file.path(dir.base,"datasets/",cohort,"/data/DNA_methylation") 
dir.data.pca <- file.path(dir.data,"/pca_filtering/") 
adni.se <- readRDS(file.path(dir.data.pca, "ADNI_QNBMIQ_PCfiltered_min_age_at_visit_65.RDS"))
# adni.se <- adni.se[,adni.se$DX %in% c("CN","Dementia")]

# Gene expression: Affymetrix Human Genome U 219 array 
ADNI_Gene_Expression_Profile <- read_csv("datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv",skip = 8)
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
genes.idx <- which(genes.low.expressed.samples.count > (length(grep("^X",colnames(ADNI_Gene_Expression_Profile))) * 0.8))
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

ADNI_Gene_Expression_Metadata <- read_csv("datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv",skip = 0,col_names = FALSE,n_max = 7)
ADNI_Gene_Expression_Metadata$X748 <- NULL

gene.exp.IDs <- apply(
  ADNI_Gene_Expression_Metadata[,4:ncol(ADNI_Gene_Expression_Metadata)],MARGIN = 2, 
  FUN = function(col) {
    paste0(
      stringr::str_extract(pattern = "[0-9]*$",string = col[3]),
      "_", col[1],"_",
      col[2]
    )
  }
)

DXSUM <- readr::read_csv("datasets/ADNI/data/Clinical/DXSUM_PDXCONV_ADNIALL_downloaded_2-19-2021.csv")
DXSUM <- DXSUM[match(gene.exp.IDs,paste0(stringr::str_extract(pattern = "[0-9]*$",DXSUM$PTID),"_",DXSUM$Phase,"_",DXSUM$VISCODE)),]
gene.exp.IDs <- paste0(DXSUM$RID,"_",DXSUM$Phase,"_",DXSUM$VISCODE2)
colnames(expression.matrix) <- gene.exp.IDs
expression.matrix <- expression.matrix[,colnames(expression.matrix) != "NA_NA_NA"]

#-----------------------------------------------------------------------------
# get matched DNAm and Gene expression
#-----------------------------------------------------------------------------
dnam.IDs <- paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE)
table(gene.exp.IDs %in% dnam.IDs) 
# FALSE  TRUE 
# 258    265 
common.ids <- base::intersect(dnam.IDs %>% as.character,gene.exp.IDs %>% as.character)
adni.se <- adni.se[,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE) %in% common.ids]
adni.se <- adni.se[,match(common.ids,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))]
expression.matrix <- expression.matrix[,colnames(expression.matrix) %in% common.ids]
expression.matrix <- expression.matrix[,match(common.ids,colnames(expression.matrix))]

table(colnames(expression.matrix) == paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))

save(
  adni.se,
  expression.matrix,
  ADNI_Gene_Expression_Metadata,
  file = "datasets/Aux/ADNI_with_MCI_matched_rna_dnam.rda"
)

#-----------------------------------------------------------------------------
# get residuals 
#-----------------------------------------------------------------------------
library(readr)
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

# - fdr significant CpGs prioritized in cross-tissue meta-analysis-CpGs
cross_tissue_meta_analysis_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0
)
cpgs_prioritized  <- cross_tissue_meta_analysis_AD_vs_CN[[5]] %>% na.omit() %>% as.character
length(cpgs_prioritized)

# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",sheet = 2)
prioritized.dmrs <- prioritized.dmrs[[4]] %>% na.omit %>% as.character
length(prioritized.dmrs) # 10

prioritized.dmrs.probes <- coMethDMR:::GetCpGsInRegion(
  prioritized.dmrs,
  arrayType = "EPIC"
) # 19 DMRs
prioritized.dmrs <- as.data.frame(prioritized.dmrs)


se.selected.cpgs <- adni.se[rownames(adni.se) %in% c(cpgs_prioritized,cpgs.ad.cn, prioritized.dmrs.probes, strsplit(combp_AD_vs_CN$Probes,";") %>% unlist()),]
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","CD8T","Mono","Neutro","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- MethReg::get_residuals(
  data = log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs))), # m-values
  metadata.samples = metadata.dnam[,c("CD8T","Mono","Neutro","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")],
  cores = 4
)

Affy_Plate <- ADNI_Gene_Expression_Metadata[7,-c(1:3)] %>% as.numeric()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit","PTGENDER")]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[5] <- "B_cells"
colnames(metadata.exp)[7] <- "CD4_T_cells"
colnames(metadata.exp)[8] <- "CD8_T_cells"

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix),
  metadata.samples = metadata.exp,
  cores = 4
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

#-----------------------------------------------------------------------------
# Save data
#-----------------------------------------------------------------------------

save(
  metadata.dnam,
  se.selected.cpgs,
  residuals.matched.met,
  metadata.exp,
  rnaseq.tf.es,
  residuals.matched.exp,
  xcell,
  file = "datasets/Aux/ADNI_with_MCI_matched_rna_dnam_residuals.rda"
)


#-------------------------------------------------------------------------------
# MethReg DMR data
#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

#  in fdr significant DMRs prioritized in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(path = "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",sheet = 2)
prioritized.dmrs <- prioritized.dmrs[[4]] %>% na.omit %>% as.character
length(prioritized.dmrs)

library(coMethDMR)
prioritized.dmrs.probes <- GetCpGsInAllRegion(
  prioritized.dmrs,
  arrayType = "EPIC"
) # 19 DMRs
prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(prioritized.dmrs.probes,FUN = function(x) paste(x,collapse = ";"))

regions <- rbind(combp_AD_vs_CN[,c("DMR","Probes")],prioritized.dmrs[,c("DMR","Probes")])

load("datasets/Aux/ADNI_with_MCI_matched_rna_dnam.rda")

#-----------------------------------------------------------------------------
# get residuals 
#-----------------------------------------------------------------------------
colnames(adni.se) <- paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE)
# adni.se.cn.ad <- adni.se[,adni.se$DX %in% c("CN","Dementia")]
se.selected.median <- plyr::ldply(
  regions$Probes,
  .fun = function(x){
    cpgs <- stringr::str_split(x,";") %>% unlist; 
    colMedians(assay(adni.se)[rownames(adni.se) %in% cpgs,])
  },.id = NULL
)
colnames(se.selected.median) <- colnames(adni.se)
rownames(se.selected.median) <- regions$DMR

metadata.dnam <- colData(adni.se) %>% as.data.frame() %>% 
  dplyr::select(
    c("DX","CD8T","Mono",
      "Neutro","CD4T","NK","B","age_at_visit",
      "PTGENDER","PlateNumber"
    )
  )

residuals.matched.met <- MethReg::get_residuals(
  data = log2(se.selected.median) / (1 - se.selected.median), # m-values
  metadata.samples = metadata.dnam,
  cores = 4
)

residuals.matched.met <- residuals.matched.met[,colnames(residuals.matched.met) %in% colnames(residuals.matched.exp)]
metadata.dnam <- metadata.dnam[rownames(metadata.dnam) %in% colnames(residuals.matched.exp),]

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))
all(rownames(metadata.dnam) == colnames(residuals.matched.met))
all(colnames(residuals.matched.exp) == colnames(rnaseq.tf.es))

save(
  metadata.dnam,
  residuals.matched.met,
  metadata.exp,
  rnaseq.tf.es,
  residuals.matched.exp,
  xcell,
  regions,
  file = "datasets/Aux/ADNI_with_MCI_matched_rna_dnam_residuals_DMR.rda"
)


