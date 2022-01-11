### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
library(tidyverse)
library(RColorBrewer)

#################################################
## load in the RDS file (it's a regular data frame)
save_fn = file.path('rdas/multiomeATAC_DLPFC_celldata_UMAPcoords.rds')
df=  readRDS(df, file = save_fn)

df = df %>% mutate(Class = case_when(grepl('EXC', Celltype1) ~'EXC',
                                     grepl('IN', Celltype1) ~'INH',
                                     TRUE ~'Glia'), 
                   Class = factor(Class, c('EXC', 'INH', 'Glia')))
head(df)

table(df$Clusters) # how many cells of each type
table(df$Class) # how many cells of each class

###############################################
# make a color palette for all the cell types
cluster_cols = c(brewer.pal(10, 'Paired'), # 10 excitatory clusters
                 brewer.pal(6, 'Set3'),  # 6 inhibitory clusters
                 brewer.pal(5, 'Dark2'))   # 5 glia clusters
names(cluster_cols) = levels(df$Clusters) # set names of colors as the 

## make the plot with ggplot2 
dir.create('plots', showWarnings = F)
pdf('plots/multiomeATAC_DLPFC_UMAP_celltypes.pdf')
## use filled point
ggplot(df, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, fill = Clusters)) + 
  geom_point(pch = 21) + 
  scale_fill_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 12) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use solid point different colors
# set alpha for transparency to see overlapping points
ggplot(df, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, color = Clusters)) + 
  geom_point(pch = 20, alpha = .6) + 
  scale_color_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 12) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

## use different solid points for class
## number mapping for shapes here: 
# http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
ggplot(df, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, fill = Clusters)) + 
  geom_point(aes(shape = Class), alpha = .6) +
  scale_shape_manual(values = 21:23, guide = 'none') + 
  scale_fill_manual(values = cluster_cols, name = '') + 
  theme_bw(base_size = 12) + 
  guides(fill = guide_legend(size = 2)) +
  theme(legend.position = 'bottom')

dev.off()



