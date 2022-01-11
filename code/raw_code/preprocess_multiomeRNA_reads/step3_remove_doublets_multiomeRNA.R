library(scds)
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)
library(here)
library(tidyverse)
library(Matrix)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

# read in srna-seq matrix files, these folders should have a barcodes.tsv, genes.tsv, and matrix.mtx files
CellRanger_fn = 'data/raw_data/CellRangerOutput' %>% 
  list.dirs(full.names = T, recursive = T) %>%
  stringr::str_subset(pattern =  'filtered_feature_bc_matrix$')
names(CellRanger_fn) = ss(CellRanger_fn, '/', 9)

sce_list <- lapply(CellRanger_fn, read10xCounts, version = '3')
sce_list <- lapply(sce_list, cxds_bcds_hybrid,list("retRes"=TRUE))
sce_list <- lapply(sce_list, function(sce) {
  mcols(sce) = mcols(sce)[,!names(mcols(sce)) %in% c('cxds_hvg_bool', 'cxds_hvg_ordr')]
  sce = sce[, colData(sce)$hybrid_score < 1.0 ]
  return(sce)
})

sce_list <- lapply(sce_list, function(sce){
  sce$Sample = gsub('_PFC', '_DLPFC', sce$Sample)
  sce$Sample = ss(sce$Sample, '/', 4)
  sce$Barcode = paste0(sce$Sample, '#', ss(sce$Barcode, '-',1))
  colnames(sce) = sce$Barcode  
  return(sce)
})

## load in raw counts of Jing's sce RNA profiles
sce_fn = file.path('data/tidy_data/rdas','JH_PFC_LabeledNuclei_20220104', 'all_nuclei_final.sce.rds')
sce_jh = readRDS(sce_fn)

sce_list <- lapply(sce_list, function(sce) {
  sce$cell_class = sce_jh$cell_class[match(sce$Barcode, colnames(sce_jh))]
  sce$cell_type2 = sce_jh$cell_type2[match(sce$Barcode, colnames(sce_jh))]
  return(sce)
})

sce_fn = here('data/tidy_data/rdas', 'raw_multiomeRNA_DLPFC_sce_list_N4.rds')
saveRDS(sce_list, sce_fn)

