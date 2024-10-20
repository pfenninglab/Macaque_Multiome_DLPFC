ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(tidyverse))
library(future)
library(scater)
library(here)

addArchRThreads(threads = 4)

######################
## load ArchR projects
DATADIR='data/tidy_data'
h5File = file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                   'all_nuclei_final_withSCTransform.h5Seurat')
obj =LoadH5Seurat(h5File)

Idents(obj) = 'cell_type2'
DefaultAssay(obj) = 'SCT'

celltypes = obj[[]] %>%
  filter(!duplicated(cell_type2)) %>% mutate(cell_type2 = as.character(cell_type2)) %>%
  filter(cell_type2 != 'TH') %>% arrange(cell_type2) %>%
  mutate(cell_class = case_when(grepl('ITGA|NR4A|TBX|ALPL', cell_type2) ~ 'EXC.Lower_ET_IT',
                                grepl('PCP4|SYT6|NKD1|POU3F1', cell_type2) ~  'EXC.Lower_CT_NP',
                                grepl('^L[2-4]', cell_type2) ~ 'EXC.Upper_IT', 
                                grepl('LAMP|NDN|VIP', cell_type2) ~ 'INH.CGE',
                                grepl('SST|PV|TH', cell_type2) ~ 'INH.MGE',
                                TRUE ~ "GLIA"), 
         cell_class = factor(cell_class, c('EXC.Upper_IT', 'EXC.Lower_CT_NP', 
                                           'EXC.Lower_ET_IT', 'INH.CGE', 'INH.MGE', 'GLIA'))) %>%
  dplyr::select(cell_type2, cell_class) %>% as.data.frame() %>% 
  arrange(cell_class, cell_type2) %>% pull(cell_type2)

celltypes_col = paletteDiscrete(celltypes)


pdf(here('figures/exploratory/celltype_specific_enhancers', 'fox_expression_monkey_pfc.pdf'))
VlnPlot(object = pbmc_small, features = 'PC_1')
dev.off()