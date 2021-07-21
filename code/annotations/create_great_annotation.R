library(rGREAT)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(S4Vectors)

epic.hg19 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
epic.hg19.gr <- epic.hg19 %>% makeGRangesFromDataFrame(
  start.field = "pos", end.field = "pos", keep.extra.columns = TRUE
)

regionsToGenes.list <- plyr::alply(
  seq(1,length(epic.hg19.gr),50000),
  .margins = 1,
  .fun = function(start){
    end <- (start + 50000 - 1)
    if(end > length(epic.hg19.gr)) end <- length(epic.hg19.gr)
    job <- submitGreatJob(epic.hg19.gr[start:end], species = "hg19")
    Sys.sleep(70)
    data.frame(plotRegionGeneAssociationGraphs(job))
  },.progress = "time")
regionsToGenes <- plyr::rbind.fill(regionsToGenes.list)

regionsToGenes$GREAT_annotation <- ifelse(
  regionsToGenes$distTSS > 0,
  paste0(regionsToGenes$gene, " (+", regionsToGenes$distTSS, ")"),
  paste0(regionsToGenes$gene, " (", regionsToGenes$distTSS, ")"))
regionsToGenes <- regionsToGenes[
  ,c("seqnames", "start", "end", "GREAT_annotation")
]

great <- regionsToGenes %>%
  group_by(seqnames, start, end) %>%
  mutate(GREAT_annotation = paste0(GREAT_annotation,collapse = ";")) %>% unique()
save(great,file = "great_EPIC_array_annotation.rda")
