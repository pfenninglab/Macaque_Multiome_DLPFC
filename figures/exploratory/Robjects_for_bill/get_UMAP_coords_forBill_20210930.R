### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
suppressMessages(library(tidyverse))
suppressMessages(library(ArchR))
library(RColorBrewer)
library(here)

#################################################
## create Arrow file and comput duplicate scores
PROJDIR= file.path('data/tidy_data/ArchRProjects','ArchR_multiomeATAC_DLPFC')
proj = loadArchRProject(path = PROJDIR)

umap_coords = getEmbedding(proj, 'UMAPI230_ATAC') %>%
  as.data.frame() %>% rownames_to_column('cellNames') %>%
  rename_with(~gsub(pattern = 'IterativeLSI230_ATAC#', replacement = '',.))

df = getCellColData(proj) %>% as.data.frame() %>%
  rename('predictedGroup_RNA2ATACCo' = 'Clusters' )%>% 
  relocate('Clusters', .after = 'Sample') %>%
  rownames_to_column('cellNames')  %>% 
  inner_join(umap_coords) %>%
  # order by Class (EXC, INH, GLIA), then by Clusters
  arrange(!grepl('EXC', Celltype1), !grepl('INH', Celltype1), Clusters) %>%
  mutate(Clusters = factor(Clusters, unique(Clusters))) %>%
  filter(Clusters != 'Drop')

# make a color palette for all the cell types
cluster_cols = c(brewer.pal(10, 'Paired'), # 10 excitatory clusters
                 brewer.pal(6, 'Set3'),  # 6 inhibitory clusters
                 brewer.pal(5, 'Dark2'))   # 5 glia clusters
names(cluster_cols) = levels(df$Clusters) # set names of colors as the 

df = df %>% mutate(col = cluster_cols[Clusters]) # assign color to dataframe

save_fn = here('figures/exploratory/Robjects_for_bill/rdas/multiomeATAC_DLPFC_celldata_UMAPcoords.rds')
saveRDS(df, file = save_fn)

