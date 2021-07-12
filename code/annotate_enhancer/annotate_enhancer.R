data <- readr::read_tsv("datasets/nasser_2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
CellType.selected <- readxl::read_xlsx("code/annotate_enhancer/41586_2021_3446_MOESM4_ESM_131-biosamples.xlsx",col_names = FALSE) %>% dplyr::pull(1)
data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
  dplyr::filter(!isSelfPromoter)  %>% 
  dplyr::filter(class != "promoter")
data.gr <- data.filtered %>% makeGRangesFromDataFrame(start.field = "start",end.field = "end",seqnames.field = "chr",keep.extra.columns = TRUE)

#-----------------------------------------------------------------------------
# overlap with the 50 cpgs + 9 dmrs sig in AD vs. CN analysis
#-----------------------------------------------------------------------------
AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 2 final_AD_vs_CN-selcted-columns-formatted.xlsx",skip = 3
)
cpgs.ad.cn <- AD_vs_CN$cpg
length(cpgs.ad.cn) # 50

cpg.gr <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations[cpgs.ad.cn,] %>%
  GenomicRanges::makeGRangesFromDataFrame(start.field = "pos",end.field = "pos")
cpg.gr <- cpg.gr + 250
cpg.gr$cpg <- cpgs.ad.cn
hits <- findOverlaps(cpg.gr,data.gr) %>% as.data.frame()

cpgs.ad.is.enahncer.nasser <- data.frame(
  "Cpg" = cpg.gr[hits$queryHits,]$cpg,
  "Cell_type" = data.gr$CellType[hits$subjectHits]
) %>% unique %>% dplyr::group_by(Cpg) %>% summarise("Cell_type" = paste(Cell_type,collapse = ";"))

# - CpGs within significant DMRs (with length > 3cpgs) identified by combp
combp_AD_vs_CN <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 2 DMRs-Combp-AD_vs_CN_annotated.xlsx",skip = 1
) # 9 DMRs
nrow(combp_AD_vs_CN)
dmr.gr <- combp_AD_vs_CN$DMR %>% MethReg::make_granges_from_names()
dmr.gr <- dmr.gr + 250
hits <- findOverlaps(dmr.gr,data.gr) %>% as.data.frame()
dmrs.ad.is.enahncer.nasser <- unique(combp_AD_vs_CN$DMR[hits$queryHits])


dmrs.ad.is.enahncer.nasser <- data.frame(
  "DMR" = combp_AD_vs_CN$DMR[hits$queryHits],
  "Cell_type" = data.gr$CellType[hits$subjectHits]
) %>% unique %>% dplyr::group_by(DMR) %>% summarise("Cell_type" = paste(Cell_type,collapse = ";"))


#-----------------------------------------------------------------------------
# overlap with 97 prioritized cpgs + 10 prioritized dmrs
#-----------------------------------------------------------------------------
cpgs.prioritized <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Supp Table 3 prioritized-CpGs-crossTissue_brain_blood.xlsx",skip = 4
)
cpgs.prioritized  <- cpgs.prioritized$CpG %>% as.character
length(cpgs.prioritized)

cpg.gr <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations[cpgs.prioritized,] %>%
  GenomicRanges::makeGRangesFromDataFrame(start.field = "pos",end.field = "pos")
cpg.gr <- cpg.gr + 250
cpg.gr$cpg <- cpgs.prioritized
hits <- findOverlaps(cpg.gr,data.gr) %>% as.data.frame()
cpgs.prioritized.is.enahncer.nasser <- data.frame(
  "Cpg" = cpg.gr[hits$queryHits,]$cpg,
  "Cell_type" = data.gr$CellType[hits$subjectHits]
) %>% unique %>% dplyr::group_by(Cpg) %>% summarise("Cell_type" = paste(Cell_type,collapse = ";"))
out <- cpgs.prioritized.is.enahncer.nasser[match(cpgs.prioritized,cpgs.prioritized.is.enahncer.nasser$Cpg),]
out$isEnahncer <- ifelse(is.na(out$Cpg),"No","Yes")
writexl::write_xlsx(x = out %>% as.data.frame(),path = "code/annotate_enhancer/cpgs.prioritized.is.enahncer.nasser.xlsx")

dmrs.prioritized <- readxl::read_xlsx(
  "DRAFT-TABLES_FIGURES_4-17-2021/_Main Table 3 Top 10 prioritized-CpGs_and_DMRs-crossTissue_brain_blood.xlsx",skip = 16,sheet = 1
)
dmrs.prioritized  <- dmrs.prioritized$DMR
length(dmrs.prioritized)

dmr.gr <- dmrs.prioritized %>% MethReg::make_granges_from_names()
dmr.gr <- dmr.gr + 250
hits <- findOverlaps(dmr.gr,data.gr) %>% as.data.frame()
dmrs.prioritized.is.enahncer.nasser <- data.frame(
  "DMR" = dmrs.prioritized[hits$queryHits],
  "Cell_type" = data.gr$CellType[hits$subjectHits]
) %>% unique %>% dplyr::group_by(DMR) %>% summarise("Cell_type" = paste(Cell_type,collapse = ";"))
