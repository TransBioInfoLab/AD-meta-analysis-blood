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
library(SummarizedExperiment)
library(coMethDMR)
library(GenomicRanges)
library(readr)
library(readxl)

#-----------------------------------------------------------------------------
path.mathReg <- "analysis_results/RNA_vs_DNAm"
dir.create("tables")
#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# CpGs with P<1E- 5 in AD vs. CN comparison
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

cpgs.prioritized <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",skip = 0
)
cpgs.prioritized  <- cpgs.prioritized[[5]] %>% na.omit() %>% as.character
length(cpgs.prioritized)

cpgs.all <- c(
  cpgs.prioritized,
  cpgs.ad.cn
) %>% unique

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/DMRs-Combp-AD_vs_CN_output_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)

#  in fdr significant DMRs in cross-tissue meta-analysis
prioritized.dmrs <- readxl::read_xlsx(
  path = "DRAFT-TABLES_FIGURES_4-17-2021/prioritization-cpgs-dmrs_5-3-2021.xlsx",
  sheet = 2
)
prioritized.dmrs <- prioritized.dmrs[[4]] %>% na.omit %>% as.character
length(prioritized.dmrs) # 10

prioritized.dmrs.probes <- plyr::alply(
  prioritized.dmrs %>% unlist,
  .margins = 1,
  .fun = function(x) GetCpGsInRegion(x, genome = "hg19", arrayType = "EPIC")
)
prioritized.dmrs <- as.data.frame(prioritized.dmrs)
colnames(prioritized.dmrs)[1] <- "DMR"
prioritized.dmrs$Probes <- sapply(
  prioritized.dmrs.probes,
  FUN = function(x) paste(x,collapse = ";")
)

regions <- rbind(combp_AD_vs_CN[,c("DMR","Probes")],prioritized.dmrs[,c("DMR","Probes")])
regions.cpgs <- strsplit(regions$Probes,split = ";") %>% unlist %>% unique
regions.gr <- regions$DMR %>% make_granges_from_names() 
df <- rbind(
  data.frame(
    "cpg" = cpgs.prioritized,
    "cpg is from" = "prioritized cpg"
  ),
  data.frame(
    "cpg" = cpgs.ad.cn,
    "cpg is from" = "ad vs. cn"
  )
)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
epic.hg19 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
epic.hg19.gr <- epic.hg19 %>% makeGRangesFromDataFrame(
  start.field = "pos", end.field = "pos",keep.extra.columns = T
)

library(dplyr)
library(GenomicRanges)

#--------------------------------------------------------------------------
# (1) Overlap with GWAS loci
#--------------------------------------------------------------------------
SNPregions <- read.csv("datasets/Aux//AD_SNPregions.csv")
SNPregions <- subset (SNPregions, select = -c(LDchr500,	LDstart500kb,	LDend500kb) )
SNPregions_ranges <- GRanges(
  seqnames = paste0("chr",SNPregions$LDchr),
  ranges = IRanges(SNPregions$LDstart, SNPregions$LDend)
) + 250000


# CpG overlap
overlap <- findOverlaps(SNPregions_ranges, epic.hg19.gr[df$cpg])
length(overlap) 

overlap_df <- as.data.frame(overlap)
SNPsOverlap <- SNPregions[overlap_df$queryHits,]
sigCpGsOverlap <- epic.hg19.gr[df$cpg][overlap_df$subjectHits,]
SNPsigCpGsOverlap_500000 <- cbind(SNPsOverlap, sigCpGsOverlap)

# DMR overlap
overlap <- findOverlaps(SNPregions_ranges,regions.gr)
length(overlap)  # 0

SNPregions <- read.csv("datasets/Aux//AD_SNPregions.csv")
SNPregions <- subset (SNPregions, select = -c(LDchr500,	LDstart500kb,	LDend500kb) )
SNPregions_ranges <- GRanges(
  seqnames = paste0("chr",SNPregions$LDchr),
  ranges = IRanges(SNPregions$LDstart, SNPregions$LDend)
) 


overlap <- findOverlaps(SNPregions_ranges, epic.hg19.gr[df$cpg])
length(overlap) 

overlap_df <- as.data.frame(overlap)
SNPsOverlap <- SNPregions[overlap_df$queryHits,]
sigCpGsOverlap <- epic.hg19.gr[df$cpg][overlap_df$subjectHits,]
SNPsigCpGsOverlap <- cbind(SNPsOverlap, sigCpGsOverlap)

# DMR overlap
overlap <- findOverlaps(SNPregions_ranges,regions.gr)
length(overlap) 

# Save
writexl::write_xlsx(
  list(
    "SNP_146_sigCpGs_Overlap_500000" = SNPsigCpGsOverlap_500000,
    "SNP_146_sigCpGs_Overlap" = SNPsigCpGsOverlap
  ),
  path = "tables/Overlap_with_AD_associate_genetics_loci.xlsx"
)

#--------------------------------------------------------------------------
# (2) Overlap with CpGs / DMRs that are associated with genetic risk scores
#--------------------------------------------------------------------------
files <- dir("datasets/Aux/Walker_DNAm_GRS_2020/",full.names = T)
data <- plyr::alply(files,.margins = 1,.fun = function(f){
  data <- readxl::read_xlsx(f)  
  if("MAPINFO" %in% colnames(data)) {
    data$CHR <- paste0("chr",data$CHR)
    data <-data %>% makeGRangesFromDataFrame(seqnames.field = "CHR",start.field = "MAPINFO",end.field = "MAPINFO")
  }
  if("BP" %in% colnames(data)) {
    data$CHR <- paste0("chr",data$CHR)
    data <-data %>% makeGRangesFromDataFrame(seqnames.field = "CHR",start.field = "BP",end.field = "BP")
  }
  
  if("End" %in% colnames(data)) {
    data <- data[!is.na(data$Chr.),]
    data$CHR <- paste0("chr",data$Chr.)
    data <- data %>% makeGRangesFromDataFrame(seqnames.field = "CHR",start.field = "Start",end.field = "End")
    
  }
  data
} )

all.cpg <- lapply(data, function(gr){
  overlap <- findOverlaps(gr, epic.hg19.gr[df$cpg])
  overlap_df <- as.data.frame(overlap)
  gr <- gr[overlap_df$queryHits]
  sigCpGsOverlap <- epic.hg19.gr[df$cpg][overlap_df$subjectHits,]
  cbind(gr %>% as.data.frame(), sigCpGsOverlap %>% as.data.frame())
})
names(all.cpg) <- paste0(gsub("\\.xlsx","",basename(files)),"_vs_146_sig_cpgs")

all.dmr <- lapply(data, function(gr){
  overlap <- findOverlaps(gr, regions.gr)
  overlap_df <- as.data.frame(overlap)
  gr <- gr[overlap_df$queryHits]
  sigDMROverlap <- regions.gr[overlap_df$subjectHits]
  cbind(gr %>% as.data.frame(), sigDMROverlap %>% as.data.frame())
})
names(all.dmr) <- paste0(gsub("\\.xlsx","",basename(files)),"_vs_19_sig_DMR")


writexl::write_xlsx(
  all.cpg,
  path = "tables/Overlap_with_AD_associate_genetics_risk_scores.xlsx"
)
