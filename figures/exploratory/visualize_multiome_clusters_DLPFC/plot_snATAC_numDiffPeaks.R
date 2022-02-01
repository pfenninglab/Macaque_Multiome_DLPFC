suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)
library(aplot)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(tidytree)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
PLOTDIR='figures/exploratory/visualize_multiome_clusters_DLPFC/plots/'
DATADIR='data/tidy_data/celltype_specific_enhancers'

#######################
## get ArchR project
PROJDIR=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
proj = loadArchRProject(path = PROJDIR)

cell_class = getCellColData(proj) %>% as.data.frame() %>%
  filter(!duplicated(Celltype2)) %>% 
  mutate(cell_class = case_when(grepl('ITGA|NR4A|TBX|ALPL', Celltype2) ~ 'EXC.Lower_ET_IT',
                                grepl('PCP4|SYT6|NKD1|POU3F1', Celltype2) ~  'EXC.Lower_CT_NP',
                                grepl('^L[2-4]', Celltype2) ~ 'EXC.Upper_IT', 
                                grepl('LAMP|NDN|VIP', Celltype2) ~ 'INH.CGE',
                                grepl('SST|PV', Celltype2) ~ 'INH.MGE',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(Celltype2, cell_class) %>% as.data.frame() %>% 
  split(x = .$Celltype2,f = .$cell_class)
cell_class = cell_class[c('EXC.Upper_IT', 'EXC.Lower_ET_IT','EXC.Lower_CT_NP','INH.CGE', 'INH.MGE',"GLIA")]
cell_class_col = setNames(RColorBrewer::brewer.pal(6, 'Dark2'), names(cell_class))

cell_class2 = list(GLIA = cell_class[['GLIA']],
                   EXC = cell_class[grep('EXC', names(cell_class))] %>% unlist(),
                   INH = cell_class[grep('INH', names(cell_class))] %>% unlist())
cell_types = unlist(cell_class)
cell_types2 = setNames(gsub('[0-9]$', '', names(cell_types)), cell_types)


### read in the tree cell type hierarchy
save_tree_fn = here(DATADIR, 'rdas/cell_type_hierarchy_multiome_DLPFC.rds')
tree_dat = readRDS(save_tree_fn)



save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
models_df = readRDS(save_models_fn)
save_peaks_fn = here(DATADIR, 'rdas/cell_type_models_diffPeakList_DLPFC.rds')
diffPeakList = readRDS(save_peaks_fn)
diffPeakList2 = lapply(diffPeakList, function(gr) with(mcols(gr), gr[Log2FC > 0]))
candidateList = lapply(diffPeakList2, names)
candidateList = split(candidateList[models_df$model], models_df$label)
candidateList = lapply(candidateList, function(ll){
  ll = ll[lengths(ll) > 1000]
  Reduce('intersect', ll)
})

clade <- MRCA(tree_dat, cell_types[!grepl('GLIA', names(cell_types))])
tree_dat2 <- tree_subset(tree_dat, clade, levels_back = 0)
tree_dat2 = groupOTU(tree_dat2,cell_class )

df_sum = data.frame(Celltype2 = names(candidateList), count = lengths(candidateList))

dir.create(here(PLOTDIR), showWarnings = F, recursive = T)
pdf(here(PLOTDIR, 'numDiffPeaks_hierarchy_barplot_DLPFC.Neurons.pdf'), height = 4, width = 7)
g = ggtree(tree_dat2, aes(color=group)) + 
  geom_tiplab(align=TRUE, color = 'black') + 
  xlim(c(0,1.5)) + 
  scale_color_manual(values = cell_class_col[-6], name = 'Grouping') + 
  guides(color = guide_legend(ncol = 1, title.position = 'top')) +
  theme(legend.position = 'bottom',
        legend.spacing.x = unit(.2, 'cm'), 
        legend.key.size = unit(.2, "cm"))+
  theme(legend.position='none')

p2 <- ggplot(df_sum, aes(x = Celltype2, y =count)) + 
  geom_text(aes(label=count, y= count), hjust = 0,angle=-0) +
  geom_col() + coord_flip() + theme_tree2() + 
  ylim(c(0, max(df_sum$count) + 1000)) +
  theme(legend.position='none')
p2 %>% insert_left(g, width=.4)
dev.off()