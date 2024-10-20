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
library(future)
library(data.table)

## for pseudo-bulk analyses
library(scuttle)
library(DESeq2)
library(limma)
library(BiocParallel)
library(IHW)

register(MulticoreParam(18))


DATADIR='data/tidy_data/celltype_specific_enhancers'

## load in the sample sheet
pd = read_csv(here('data/raw_data/tables/SampleSheet_RNA_DLPFC.csv')) %>%
  rename_with(make.names) %>% 
  rename_all(funs(stringr::str_replace_all(., '^X\\.\\.', ''))) %>%
  mutate(Animal = paste0('Monkey_',stringr::str_sub(Nickname, start = 1, end = 1)),
         Sample = paste0(Nickname, '_', Region, '-',Batch)) %>%
  relocate(Sample, Region, Animal, Batch, .before = everything())

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 16)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

###############################################
## 1) load SeuratObjects of labeled all nuclei
sce = readRDS(file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                             'all_nuclei_final.sce.rds'))

sce_agg <- aggregateAcrossCells(sce, id=colData(sce)[,c("cell_type2", "Sample")])
assays(sce_agg)$counts = round(assays(sce_agg)$counts)
with(colData(sce_agg), table(cell_type2))
rm(sce)

## calculate the detection rate of genes in this cell type
y <- assays(sce_agg)$counts
sce_agg$cdr <- scale(colMeans(y > 0))


dds <- DESeqDataSet(sce_agg, design = ~ 0 + cell_type2 + Sample + cdr)
dds = DESeq(dds, parallel = TRUE)

###############################################
## 2) read in the differential peak comparisons
save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_models_fn)
model_name = cell_type_df$model
indlist = setNames(seq_along(model_name), model_name)
fgd_labels = setNames(cell_type_df$label, model_name)
bgd_labels = setNames(cell_type_df$bgd.labels %>% sapply(strsplit,';'), model_name)
bgd_labels = lapply(bgd_labels, function(b) b[b %in% levels(dds$cell_type2)])

contrastsList = mapply(function(f, g){
  paste(f , '-(', paste(g, collapse = '+'), ')/', length(g) )
}, f = fgd_labels, g = bgd_labels)

##############################################################
## 2) compute differential marker genes across all comparisons
de.markers.list = lapply(contrastsList, function(cont){
  print(paste('Getting DEGs for', cont))
  contrasts =  c(makeContrasts(contrasts = cont, levels = levels(dds$cell_type2)), rep(0,8))
  results(dds, filterFun=ihw, contrast =contrasts)
})
names(de.markers.list) = model_name
save_DEGs_fn = here(DATADIR, 'rdas/cell_type_models.all_DEGs.DESeq2.rds')
#saveRDS(de.markers.list %>% lapply(as.data.frame), save_DEGs_fn)

##############################################################
## 4) filter to significant DEGs, 
alpha = 0.05
de.markers.list2 = de.markers.list %>% lapply(as.data.frame) %>%
  lapply(rownames_to_column, "Gene.Symbol") %>% rbindlist(idcol = 'model') %>%
  mutate(Gene.Symbol = as.character(Gene.Symbol)) %>%
  filter(padj < alpha, log2FoldChange > 2) %>%
  full_join(x = cell_type_df) %>%
  group_by(label, Gene.Symbol) %>% top_n(1, -log10(padj)) %>% ungroup() %>% 
  group_by(Gene.Symbol) %>% filter(n() <= 6) %>%  ungroup() %>%
  dplyr::select(-c('model', 'bgd.group', 'bgd.labels')) %>%
  arrange(Gene.Symbol) %>% split(f = .$label)
  
sapply(de.markers.list2, nrow)
sapply(de.markers.list2, '[[', 'Gene.Symbol' )  %>% unlist() %>% head()

save_DEGs_fn2 = here(DATADIR, 'rdas/cell_type_models.top_DEGs.DESeq2.rds')
print('done')
#saveRDS(de.markers.list2, save_DEGs_fn2)

