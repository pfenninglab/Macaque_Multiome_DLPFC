### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
suppressMessages(library(tidyverse))
suppressMessages(library(ArchR))
library(RColorBrewer)
library(here)

PLOTDIR= here('figures/exploratory/Robjects_for_bill','plots')

#################################################
## create Arrow file and comput duplicate scores
PROJDIR= here('data/tidy_data/ArchRProjects','ArchR_DLPFC_scATAC')
proj = loadArchRProject(path = PROJDIR)


umap_coords = getEmbedding(proj, 'UMAPH200_ATAC') %>%
  as.data.frame() %>% rownames_to_column('cellNames') %>%
  rename_with(~gsub(pattern = 'HarmonyI200_ATAC#', replacement = '',.)) %>%
  rename_with(~gsub(pattern = 'UMAP_Dimension_', replacement = 'UMAP_Dim_',.)) 


df = getCellColData(proj)
df = df %>% as.data.frame() %>% 
  mutate(Class = case_when(grepl('EXC', Celltype1) ~'EXC',
                           grepl('IN', Celltype1) ~'INH',
                           TRUE ~'Glia'), 
         Class = factor(Class, c('EXC', 'INH', 'Glia'))) %>%
  arrange(Class, Celltype2) %>%
  rownames_to_column('cellNames') %>% inner_join(umap_coords) %>%
  mutate(Celltype2 = factor(Celltype2, unique(Celltype2)))


###############################################
# make a color palette for all the cell types
cluster_cols = c(brewer.pal(11, 'Paired'), # 10 excitatory Celltype2
                 brewer.pal(6, 'Set3'),  # 6 inhibitory Celltype2
                 brewer.pal(6, 'Dark2'))   # 6 glia Celltype2
names(cluster_cols) = levels(df$Celltype2) # set names of colors as the 

## make the plot with ggplot2 
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
pdf(here(PLOTDIR, 'Macaque_DLPFC_scATAC_UMAP.pdf'), height = 5, width = 4)
## use filled point
ggplot(df, aes(x = UMAP_Dim_1, y = UMAP_Dim_2,  fill = Celltype2)) + 
  geom_point(pch = 21) + 
  scale_fill_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 8) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use solid point different colors
# set alpha for transparency to see overlapping points
ggplot(df, aes(x = UMAP_Dim_1, y = UMAP_Dim_2, color = Celltype2)) + 
  geom_point(pch = 20, alpha = .6) + 
  scale_color_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 8) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use different solid points for class
## number mapping for shapes here: 
# http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
ggplot(df, aes(x = UMAP_Dim_1, y = UMAP_Dim_2, fill = Celltype2)) + 
  geom_point(aes(shape = Class), alpha = .6) +
  scale_shape_manual(values = 21:23, guide = 'none') + 
  scale_fill_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 8) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

dev.off()






## make the plot with ggplot2 

###############################################
# make a color palette for all the cell types
cluster_cols2 = cluster_cols[1:11]
table(df$Celltype2)

pdf(here(PLOTDIR, 'Macaque_DLPFC_scATAC_UMAP_EXC.pdf'), height = 5, width = 4)
## use filled point
ggplot(df%>% filter(Class == 'EXC'), aes(x = UMAP_Dim_1, y = UMAP_Dim_2,  fill = Celltype2)) + 
  geom_point(pch = 21) + 
  scale_fill_manual(values = cluster_cols2, name = '') + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use solid point different colors
# set alpha for transparency to see overlapping points
ggplot(df %>% filter(Class == 'EXC'), aes(x = UMAP_Dim_1, y = UMAP_Dim_2, color = Celltype2)) + 
  geom_point(pch = 20, alpha = .6) + 
  scale_color_manual(values = cluster_cols2, name = '') + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use different solid points for class
## number mapping for shapes here: 
# http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
ggplot(df%>% filter(Class == 'EXC'), aes(x = UMAP_Dim_1, y = UMAP_Dim_2, fill = Celltype2)) + 
  geom_point(aes(shape = Class), alpha = .6) +
  scale_shape_manual(values = 21:23, guide = 'none') + 
  scale_fill_manual(values = cluster_cols2, name = '') + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

dev.off()






umap_coords = getEmbedding(proj, 'UMAPI230_ATAC') %>%
  as.data.frame() %>% rownames_to_column('cellNames') %>%
  rename_with(~gsub(pattern = 'IterativeLSI230_ATAC#', replacement = '',.))

df = getCellColData(proj) %>% as.data.frame() %>%
  rename('predictedGroup_RNA2ATACCo' = 'Celltype2' )%>% 
  relocate('Celltype2', .after = 'Sample') %>%
  rownames_to_column('cellNames')  %>% 
  inner_join(umap_coords) %>%
  # order by Class (EXC, INH, GLIA), then by Celltype2
  arrange(!grepl('EXC', Celltype1), !grepl('INH', Celltype1), Celltype2) %>%
  mutate(Celltype2 = factor(Celltype2, unique(Celltype2))) %>%
  filter(Celltype2 != 'Drop')

# make a color palette for all the cell types
cluster_cols = c(brewer.pal(10, 'Paired'), # 10 excitatory Celltype2
                 brewer.pal(6, 'Set3'),  # 6 inhibitory Celltype2
                 brewer.pal(5, 'Dark2'))   # 5 glia Celltype2
names(cluster_cols) = levels(df$Celltype2) # set names of colors as the 

df = df %>% mutate(col = cluster_cols[Celltype2]) # assign color to dataframe

save_fn = here('figures/exploratory/Robjects_for_bill/rdas/multiomeATAC_DLPFC_celldata_UMAPcoords.rds')
saveRDS(df, file = save_fn)

