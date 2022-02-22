ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(tidyverse))
library(here)

addArchRThreads(threads = 4)

######################
## load ArchR projects
DATADIR='data/tidy_data'
obj =LoadH5Seurat(file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                                 'all_nuclei_final.h5Seurat'))
obj_neuron =LoadH5Seurat(file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                                 'neuron_final.h5Seurat'))


cell_class = obj[[]] %>%
  filter(!duplicated(cell_type2)) %>% mutate(cell_type2 = as.character(cell_type2)) %>%
  filter(cell_type2 != 'TH') %>% arrange(cell_type2) %>%
  mutate(cell_class = case_when(grepl('ITGA|NR4A|TBX|ALPL', cell_type2) ~ 'EXC.Lower_ET_IT',
                                grepl('PCP4|SYT6|NKD1|POU3F1', cell_type2) ~  'EXC.Lower_CT_NP',
                                grepl('^L[2-4]', cell_type2) ~ 'EXC.Upper_IT', 
                                grepl('LAMP|NDN|VIP', cell_type2) ~ 'INH.CGE',
                                grepl('SST|PV|TH', cell_type2) ~ 'INH.MGE',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(cell_type2, cell_class) %>% as.data.frame() %>% 
  split(x = .$cell_type2,f = .$cell_class)
cell_class = cell_class[c('EXC.Upper_IT', 'EXC.Lower_ET_IT','EXC.Lower_CT_NP','INH.CGE', 'INH.MGE',"GLIA")]
cell_class_col = setNames(RColorBrewer::brewer.pal(6, 'Dark2'), names(cell_class))
cell_class_col2 = mapply(rep, cell_class_col, lengths(cell_class)) %>% unlist()
names(cell_class_col2) = unlist(cell_class)
cell_class_col2 = cell_class_col2[!is.na(names(cell_class_col2))]


obj@meta.data$cell_type2 = factor(obj@meta.data$cell_type2, names(cell_class_col2))
obj_neuron@meta.data$cell_type2 = factor(obj_neuron@meta.data$cell_type2, names(cell_class_col2)[1:17])
Idents(obj) = 'cell_type2'
Idents(obj_neuron) = 'cell_type2'

############################################
## get markerPeaks scored and ranked by CNNs
DATADIR='data/tidy_data/celltype_specific_enhancers'
save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
markerGeneList = readRDS(save_top_enh_fn) %>% lapply('[[', 'markerGene') %>%
  lapply(function(x) x[!is.na(x)])
markerGeneList = markerGeneList[lengths(markerGeneList) > 0]

FIGDIR='figures/exploratory/celltype_specific_enhancers'
dir.create(here(FIGDIR, 'plots', 'geneplots'), showWarnings = F)

for(name in names(markerGeneList)){
  plot_fn = paste0('DLPFC_markerGene_SNAIL_candidate_ranked_AllTracks.', name, '.20220214.pdf')
  out_pdf = here(FIGDIR, 'plots', 'geneplots',plot_fn)

  ggList = VlnPlot(obj, features = markerGeneList[[name]], assay = 'RNA', 
                   cols = cell_class_col2, combine = FALSE, pt.size = FALSE) 

  pdf(out_pdf, height = 4, width = 8)
  for(gg in ggList){
    print(gg+ theme(legend.position = 'none') )}
  dev.off()
}




