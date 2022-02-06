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
obj = LoadH5Seurat(file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                             'all_nuclei_final.h5Seurat'), assay = 'RNA')

# fix some NA cell type labels ~30% of cells
tbl = with(obj[[]], table(cell_type2, seurat_clusters))
tmp = apply(tbl, 2, function(x) sum(x>10))
clust2celltype = setNames(rownames(tbl)[apply(tbl, 2, which.max)], colnames(tbl))
clust2celltype = clust2celltype[apply(tbl, 2, function(x) sum(x>10)) == 1]
obj@meta.data$cell_type2 = 
  with(obj[[]], ifelse(!is.na(cell_type2),as.character(cell_type2), 
                       clust2celltype[as.character(seurat_clusters)]))
Idents(obj) = 'cell_type2' ## set for Seurat to know what to group cells for FindMarkers


###############################################
## 2) read in the differential peak comparisons
save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_models_fn)
model_name = cell_type_df$model
indlist = setNames(seq_along(model_name), model_name)
fgd_labels = setNames(cell_type_df$label, model_name)
bgd_labels = setNames(cell_type_df$bgd.labels %>% sapply(strsplit,';'), model_name)
bgd_labels = sapply(bgd_labels, function(x) x[x %in% levels(Idents(obj))])


##############################################################
## 2) compute differential marker genes across all comparisons
de.markers.list = lapply(indlist, function(idx){
  print(paste('Compute differential genes on:', names(indlist)[idx]))
  ret <- FindMarkers(obj, ident.1 =fgd_labels[idx], ident.2 = bgd_labels[[idx]], 
                     max.cells.per.ident = 1000, min.pct = 0.25, test.use = "wilcox")
  return(ret)
})

save_DEGs_fn = here(DATADIR, 'rdas/cell_type_models.all_DEGs.rds')
saveRDS(de.markers.list, save_DEGs_fn)

##############################################################
## 4) filter to significant DEGs, 
alpha = 0.05
de.markers.list2 = de.markers.list %>% 
  lapply(rownames_to_column, "Gene.Symbol") %>% rbindlist(idcol = 'model') %>%
  mutate(Gene.Symbol = as.character(Gene.Symbol), 
         p_val_adj = p.adjust(p_val, 'fdr')) %>%
  filter(p_val_adj < alpha, avg_log2FC > 0) %>%
  full_join(x = cell_type_df) %>%
  group_by(label, Gene.Symbol) %>% top_n(1, -log10(p_val_adj)) %>% ungroup() %>% 
  group_by(Gene.Symbol) %>% filter(n() <= 6) %>%  ungroup() %>%
  dplyr::select(-c('model', 'bgd.group', 'bgd.labels', 'pct.1', 'pct.2')) %>%
  arrange(Gene.Symbol) %>% split(f = .$label)
  
sapply(de.markers.list2, nrow)
sapply(de.markers.list2, '[[', 'Gene.Symbol' )  %>% unlist() %>% head()

save_DEGs_fn2 = here(DATADIR, 'rdas/cell_type_models.top_DEGs.rds')
saveRDS(de.markers.list2, save_DEGs_fn2)

