suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)
library(ggtree)
library(ggtreeExtra)
library(tidytree)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }


########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

PLOTDIR='figures/exploratory/cell_type_hierarchy/plots/'


######################
## get ArchR project
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_Multiome_DLPFC")
proj = loadArchRProject(path = PROJDIR2)
proj = proj[! proj$predictedGroup_RNA2RNACo %in% c('TH', 'L5.POU3F1')]
table(proj$predictedGroup_RNA2RNACo)

celltypes = getCellColData(proj) %>% as.data.frame() %>%
  filter(!duplicated(predictedGroup_RNA2RNACo)) %>% 
  mutate(cell_class = case_when(grepl('^L[2-6]', predictedGroup_RNA2RNACo) ~ 'EXC', 
                                grepl('SST|LAMP|NDN|VIP|PVAL', predictedGroup_RNA2RNACo) ~ 'INH',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(predictedGroup_RNA2RNACo, cell_class) %>% as.data.frame() %>% 
  split(x = .$predictedGroup_RNA2RNACo,f = .$cell_class)
unlist(celltypes[c('EXC', 'INH', 'GLIA')])

###############################
## compute differential peaks
diffPeak_fn = here('figures/exploratory/cell_type_hierarchy/rdas/differentialPeaks_1vAll.rds')
if(FALSE){
  markersPeaks <- getMarkerFeatures( proj,useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_RNA2RNACo",  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"))
  saveRDS(markersPeaks, diffPeak_fn)
} else {
  markersPeaks = readRDS(diffPeak_fn)
}

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

markerDat = markerList %>% lapply(as.data.frame) %>%
  rbindlist(idcol = 'celltype') %>% 
  mutate(celltype = factor(celltype, unlist(celltypes[c('EXC', 'INH', 'GLIA')])))
markerDat2 = markerDat %>% group_by(idx) %>% filter(length(idx) > 1) %>%
  ungroup()

tab1 = full_join(
  getCellColData(proj) %>% as.data.frame() %>%
    mutate(celltype = factor(predictedGroup_RNA2RNACo, unlist(celltypes[c('EXC', 'INH', 'GLIA')]))) %>% count(celltype) %>% rename('numCells' = 'n'),
  markerDat %>% count(celltype) %>% rename('differentialEnhancer' = 'n')) %>%
  full_join(markerDat2 %>% count(celltype) %>% rename('uniqDiffEnhancer' = 'n'))

markerDat %>% writexl::write_xlsx(here('figures/exploratory/cell_type_hierarchy/tables/differentialEnhancerTable_multiomeATAC_DLPFC.xlsx'))

markerDat2 %>% writexl::write_xlsx(here('figures/exploratory/cell_type_hierarchy/tables/uniqDiffEnhancerTable_multiomeATAC_DLPFC.xlsx'))

tab1 %>% writexl::write_xlsx(here('figures/exploratory/cell_type_hierarchy/tables/candidateCelltypeEnhancerSummaryTable_multiomeATAC_DLPFC.xlsx'))


