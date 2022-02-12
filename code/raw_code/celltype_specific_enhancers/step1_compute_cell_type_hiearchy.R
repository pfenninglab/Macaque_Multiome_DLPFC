
suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(tidytree)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
PLOTDIR='figures/exploratory/celltype_specific_enhancers/plots/'
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

#############################################################################
# ## get the data-driven cell type-hiearchy from the join RNA-ATAC embedding
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_multiomeATAC")
proj2 = loadArchRProject(path = PROJDIR2)
table(proj2$Celltype2)

## embedding in the dimensionality reduction, use harmony corrected version
dimRed  = getReducedDims(proj2, 'HarmonyI_Combined') %>% as.data.frame()
colnames(dimRed) = paste0(colnames(dimRed), seq(ncol(dimRed)))

## get the average high dim distance per cell type
dimRed = cbind(dimRed, getCellColData(proj2) %>% as.data.frame()) %>% 
  group_by(Celltype2) %>%
  summarise_if(is.numeric, mean) %>% 
  column_to_rownames('Celltype2') %>%
  select(starts_with('LS')) %>% as.matrix()

## make into a distance matrix and then a tree
Matrix <- as.matrix(t(dimRed))
distMat <- as.dist(1 - cor(Matrix, method = "pearson"))

set.seed(101)
my_tree <- ape::nj(distMat)
my_tree = root(my_tree,outgroup =  'Microglia')
tree_dat = as.treedata(my_tree) %>% groupOTU(cell_class, "cell_class")

pdf(here(PLOTDIR, 'cell_type_hierarchy_multiome_DLPFC.pdf'), height = 4, width = 2)
ggtree(tree_dat, aes(color=cell_class)) + 
  geom_tiplab(align=TRUE, color = 'black') + 
  scale_color_manual(values = cell_class_col, name = 'Grouping') + 
  xlim(c(0,1.5)) + 
  guides(color = guide_legend(ncol = 1, title.position = 'top')) +
  theme(legend.position = 'bottom',
        legend.spacing.x = unit(.2, 'cm'), 
        legend.key.size = unit(.2, "cm"))
dev.off()

## save the tree
dir.create(here(DATADIR, 'rdas'), recursive = T, showWarnings = F)
dir.create(here(DATADIR, 'tables'), recursive = T, showWarnings = F)
save_tree_fn = here(DATADIR, 'rdas/cell_type_hierarchy_multiome_DLPFC.rds')
saveRDS(tree_dat, save_tree_fn)


groupings = c(cell_class, cell_class2)
df = cbind(do.call("rbind", replicate(6, as_tibble(tree_dat) , simplify = FALSE)),
           data.frame(bgd1 = rep(names(cell_class), each = nrow(as_tibble(tree_dat))))) %>%
  mutate(tmp = cell_class %>% as.character() %>% ss('\\.'),
         tmp2 = bgd1 %>% ss('\\.')) %>%
  filter(!is.na(label), tmp == tmp2, tmp2 != 'GLIA') %>% dplyr::select(-c(tmp, tmp2))
  
cell_type_df= cbind(do.call("rbind", replicate(length(cell_class2), df, simplify = FALSE)),
        data.frame(bgd2 = rep(names(cell_class2), each = nrow(df)))) %>%
  select(label, starts_with('bgd')) %>%
  pivot_longer(-label, values_to = 'bgd.group') %>%
  mutate(bgd.labels = map2_chr(bgd.group, label, function(bgd.group, label){
    tmp = groupings[bgd.group] %>% unlist() %>% sort()
    tmp = tmp[tmp != label]
    tmp = paste(unlist(tmp), collapse = ';')
    return(tmp)})) %>%
  arrange(label, name, bgd.group) %>% dplyr::select(-name)%>%
  filter(!duplicated(paste0(label,bgd.labels))) %>%
  mutate(model = paste0(label, 'vs', bgd.group)) %>%
  relocate(model, .before = everything())

save_model_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
saveRDS(cell_type_df, save_model_fn)

