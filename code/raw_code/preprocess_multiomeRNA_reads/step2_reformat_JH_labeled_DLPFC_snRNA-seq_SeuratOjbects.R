# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(here)

DATADIR='data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104'

## load in the sample sheet
pd = read_csv(here('data/raw_data/tables/SampleSheet_RNA_DLPFC.csv')) %>%
  rename_with(make.names) %>% 
  rename_all(funs(stringr::str_replace_all(., '^X\\.\\.', ''))) %>%
  mutate(Animal = paste0('Monkey_',stringr::str_sub(Nickname, start = 1, end = 1)),
         Sample = paste0(Nickname, '_', Region, '-',Batch)) %>%
  relocate(Sample, Region, Animal, Batch, .before = everything())

#######################################################
## 1) load SeuratObjects of labeled all nuclei, glia
obj = readRDS(file.path(DATADIR, 'nuclei.integrated_1200_7000_30_30_neuron2200_final_2_updated.rds'))
# obj2 = obj

## relabel cell phenotype data to match ATAC
JHbatch_to_BNPbatch = setNames(c(1,2,4:9), c(4:7, 0:3))
obj@meta.data = obj[[]] %>% 
  mutate(Batch = JHbatch_to_BNPbatch[batch], 
         Sample = paste0(monkey, '_DLPFC-', Batch)) %>%
  select(-c(batch, sample_names, monkey)) %>%
  left_join(pd) %>%
  relocate(all_of(names(pd)), .before = everything())

## relabel cell barcodes to match ATAC
obj = RenameCells(obj, new.names = paste0(obj$Sample, '#', ss(Cells(obj), '-', 1)))

## relabel glia to be consistent and relevel
gliaTypes = c('Astrocytes', 'Endothelial', 'Oligos', 'Oligo_Pre', 'Microglia')
gliaTypes2 = setNames(c('Astro', 'Endo', 'Oligo','OPC', 'Microglia'), gliaTypes)

obj$cell_class = ifelse(as.character(obj$cell_class) %in% gliaTypes,
             gliaTypes2[as.character(obj$cell_class)], as.character(obj$cell_class))
obj$cell_class = factor(obj$cell_class, levels = c('Excita_neurons', 'Interneurons', gliaTypes2))

#######################################
## 2) subset glia dataset of all nuclei
obj_glia = subset(obj, subset = cell_class %in% gliaTypes2 )
sce_glia = as.SingleCellExperiment(obj_glia, assay = 'RNA')
saveRDS(sce_glia, file.path(DATADIR, 'glia_final.sce.rds'))
SaveH5Seurat(obj_glia, file.path(DATADIR, 'glia_final.h5Seurat'), overwrite = TRUE)

#####################################################################
## 3) load SeuratObjects of labeled excitatory and inhibitory neurons
obj_exc = readRDS(file.path(DATADIR, 'excita_2_final_updated.rds'))
obj_inh = readRDS(file.path(DATADIR, 'interneurons_1_updated.rds'))
obj_neuron <- merge(obj_exc, y = obj_inh, add.cell.ids = c("EXC", "INH"), 
                    project = "DLPFC_snRNA")
# obj_neuron2 = obj_neuron

## fix neuron phenotype labels
cols2match = c('n_genes', 'nCount_RNA', 'percent_ribo')
indNeur = which(obj$cell_class %in% c( 'Interneurons','Excita_neurons'))
JH_to_BNP_samples = inner_join(obj_neuron[[]][,c(cols2match, 'sample_names')], 
                 obj[[]][indNeur,c(cols2match, 'Sample')]) %>%
  filter(!duplicated(sample_names)) %>% select(sample_names, Sample) %>% deframe()


## add BNP sample labels
obj_neuron@meta.data = obj_neuron[[]] %>% 
  mutate(Sample = JH_to_BNP_samples[sample_names]) %>%
  select(-c(batch, sample_names, monkey)) %>%
  left_join(pd) %>%
  relocate(all_of(names(pd)), .before = everything())

## relabel cell barcodes to match ATAC
obj_neuron = RenameCells(obj_neuron, new.names = paste0(obj_neuron$Sample, '#', ss(Cells(obj_neuron), '_|-', 2)))

## relevel the cell type labels, group the various sub groups together
obj_neuron$cell_type2 =
  obj_neuron$cell_type %>% gsub(' [1-3]', '', .) %>% gsub('^PV', 'PV.', .) %>%
  gsub('\\+', '', .) %>% gsub('^SST.*', 'SST', .) %>% make.names()

## save the neuron dataset
sce_neuron = as.SingleCellExperiment(obj_neuron, assay = 'RNA')
saveRDS(sce_neuron, file.path(DATADIR, 'neuron_final.sce.rds'))
SaveH5Seurat(obj_neuron, file.path(DATADIR, 'neuron_final.h5Seurat'), overwrite = TRUE)


#####################################################################
## add neuronal subtype labels to the full nuclei dataset
obj$cell_type2 = obj_neuron$cell_type2[match(Cells(obj), Cells(obj_neuron))]
obj$cell_type2 = ifelse(is.na(obj$cell_type2), as.character(obj$cell_class), obj$cell_type2) %>%
  factor(c(unique(obj_neuron$cell_type2) %>% sort(), gliaTypes2))

## save this dataset in SingleCellExperiment and h5Seurat
sce = as.SingleCellExperiment(obj, assay = 'RNA')
saveRDS(sce, file.path(DATADIR, 'all_nuclei_final.sce.rds'))
SaveH5Seurat(obj, file.path(DATADIR, 'all_nuclei_final.h5Seurat'), overwrite = TRUE)
