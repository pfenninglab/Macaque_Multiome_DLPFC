suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(RColorBrewer)
library(tidytree)
library(aplot)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
PLOTDIR='figures/exploratory/celltype_specific_enhancers/plots/'
DATADIR='data/tidy_data/celltype_specific_enhancers'

in2mm = 25.4

#################################
## 1) grab the scATAC monkey data
proj = loadArchRProject(here('data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC'))

## sort the cell types
cellTypes = getCellColData(proj) %>% as.data.frame()%>%
  arrange(desc(grepl('^L[2-3].[A-Z]', Celltype2)), desc(grepl('^L[2-3]', Celltype2)), 
          desc(grepl('^L[4-5].[A-Z]', Celltype2)), desc(grepl('^L[4-5]', Celltype2)), 
          desc(grepl('^L6.[A-Z]', Celltype2)), desc(grepl('^L6.', Celltype2)), 
          desc(grepl('SST|PV', Celltype2)),
          desc(grepl('LAMP5|NDNF|VIP|PAX', Celltype2)), 
          Celltype2) %>%
  pull(Celltype2) %>% unique()

## set the colors for the joint cell types
exc_col = setNames(c(brewer.pal(10, 'Paired'), '#D56DA9'),
                   c('L2.CUX2.MEIS2', 'L4.5.TBX15', 'L4.ALPL', 'L3.CUX2.RORB',
                     'L6.ITGA8','L4.TYR', 'L6.NKD1', 'L6.SYT6', 'L5.PCP4', 
                     'L5.6.NR4A2', 'L5.POU3F1'))
inh_col = setNames(brewer.pal(8, 'Dark2'), cellTypes[12:17])
glia_col = setNames(brewer.pal(9, 'Set3'), cellTypes[18:23])
cellTypes_cols = c(exc_col, inh_col, glia_col) [cellTypes]


#####################################
## 1) read in the cell type hierarchy
tree_dat = here(DATADIR, 'rdas/cell_type_hierarchy_multiome_DLPFC.rds') %>% readRDS(save_tree_fn)
enhancer_df = here(DATADIR, 'rdas/rheMac10_DLPFC.candidate_celltype_enhancers.rds') %>% 
  readRDS()%>% 
  lapply(function(x) x %>% dplyr::select(-starts_with('score_'), -label) )%>% 
  data.table::rbindlist(idcol = 'label',fill=TRUE) 

enhancer_df_agg = enhancer_df %>% filter(compositeScore>0) %>% group_by(label) %>% 
  summarise(numEnhancer = n()) %>% 
  mutate(label = factor(label, names(cellTypes_cols) ))


pdf(here(PLOTDIR, 'cell_type_hierarchy_multiome_DLPFC_sfn2022.pdf'), 
    height = 240/in2mm, width = 240/in2mm)
ggtree(tree_dat) + 
  geom_tiplab(align=TRUE, aes(color = label)) +  hexpand(.3)+
  geom_facet(panel = "# of candidate enhancers", data = enhancer_df_agg, geom = geom_col, 
               width = 1, aes( x = numEnhancer, fill = label), 
               orientation = 'y') +
  scale_color_manual(values = cellTypes_cols, guide = 'none')+
  scale_fill_manual(values = cellTypes_cols, guide = 'none') +
  theme_bw(base_size = 45) + 
  theme_tree2(legend.position=c(.05, .85)) + 
  theme(strip.text.x = element_text(size = 26))
dev.off()

