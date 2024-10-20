library(SeuratDisk)
library(Seurat)
library(metap)
library(here)
library(tidyverse)
library(future)

plan("multicore", workers = 12)
PLOTDIR='figures/exploratory/celltype_markers_L5PT'

dir.create(here(PLOTDIR, 'plots'))
dir.create(here(PLOTDIR, 'tables'))
dir.create(here(PLOTDIR, 'rdas'))

## read in human-mouse ortholog gene list
hg2mm_df = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv') %>% 
  rename_with(make.names)

############################################
## 1) find macaque conserved L5PT marker genes
obj_rm = LoadH5Seurat(here('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104/all_nuclei_final_withSCTransform.h5Seurat'))

obj_rm$cell_type2

head(obj_rm[[]])
table(obj_rm$Sample)
table(obj_rm$cell_type2)

## find L5.ET marker genes vs. all other cell types across each animal ID
## combine w/ the largest p-value
markers_rm = FindConservedMarkers(
  obj_rm, ident.1 = 'L5.POU3F1', grouping.var = 'Sample', 
  assay = "RNA", slot = "data", min.cells.group = 3, 
  meta.method = metap::maximump, verbose = TRUE)

## filter out data frame w/ extraneous columns
markers_rm2 = markers_rm %>% rownames_to_column('Gene.name') %>% 
  dplyr::select(-contains('DLPFC')) %>% 
  mutate(min_pbonf = p.adjust(minimump_p_val, 'bonferroni')) %>% 
  inner_join(x = hg2mm_df) %>% arrange(minimump_p_val)


et_markers_fn = here(PLOTDIR, 'tables', 'JH_PFC_LabeledNuclei_20220104_L5ET_neuron_markers.xlsx')
markers_rm2 %>% writexl::write_xlsx(et_markers_fn)


############################################
## 2) find mouse conserved L5PT marker genes

obj_mm = LoadH5Seurat('/projects/pfenninggroup/singleCell/allen_mouse_brain/data/SeuratProjects/aibs_mouse_ctx_10x_combined_subclass_label.downSample100k.h5Seurat')
Idents(obj_mm) = 'combined_subclass_label'

head(obj_mm[[]])
table(obj_mm$external_donor_name_id)
table(obj_mm$combined_subclass_label)

## find L5.ET marker genes vs. all other cell types across each animal ID
## combine w/ the largest p-value
markers_mm = FindConservedMarkers(
  obj_mm, ident.1 = 'L5.ET', grouping.var = 'external_donor_name_id', 
  assay = "RNA", slot = "data", min.cells.group = 3, 
  meta.method = metap::maximump, verbose = TRUE)

## filter out data frame w/ extraneous columns
markers_mm2 = markers_mm %>% rownames_to_column('Mouse.gene.name') %>% 
  dplyr::select(-matches('^[0-9]')) %>% 
  mutate(min_pbonf = p.adjust(minimump_p_val, 'bonferroni')) %>% 
  inner_join(x = hg2mm_df) %>% arrange(minimump_p_val)


et_markers_fn2 = here(PLOTDIR, 'tables', 'aibs_mouse_ctx_L5ET_neuron_markers.xlsx')
markers_mm2 %>% writexl::write_xlsx(et_markers_fn2)


###########################
## 3) combine the 2 lists
alpha = 0.05
markers_shared = inner_join(
  markers_rm2 %>% filter(max_pval < alpha) %>% 
    dplyr::rename('max_pval_rm' = 'max_pval', 'minimump_p_val_rm' = 'minimump_p_val', 
                   'min_pbonf_rm' = 'min_pbonf'), 
  markers_mm2 %>% filter(max_pval < alpha) %>% 
    dplyr::rename('max_pval_mm' = 'max_pval', 'minimump_p_val_mm' = 'minimump_p_val', 
                  'min_pbonf_mm' = 'min_pbonf')) %>% 
  arrange(desc( log(max_pval_mm)*log(max_pval_rm)))
  
  
et_markers_fn3 = here(PLOTDIR, 'tables', 'mouse-monkey_shared_L5ET_neuron_markers.xlsx')
markers_shared %>% writexl::write_xlsx(et_markers_fn3)


