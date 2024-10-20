suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(rtracklayer)
library(future)
library(harmony)
library(pheatmap)

## for plotting 
library(ggsankey)
library(RColorBrewer)

PLOTDIR='figures/exploratory/jing_trackPlots/plots'
in2mm = 25.4

#################################
## 1) grab the scATAC monkey data
proj = loadArchRProject('/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC')

## sort the cell types
cellTypes = getCellColData(proj) %>% as.data.frame()%>%
  pivot_longer(cols = c(Celltype2, allen_labels), 
               names_to = 'celltype', values_to = 'values') %>% 
  arrange(desc(grepl('^L[2-3].[A-Z]', values)), desc(grepl('^L[2-3]', values)), 
          desc(grepl('^L[4-5].[A-Z]', values)), desc(grepl('^L[4-5]', values)), 
          desc(grepl('^L6.[A-Z]', values)), desc(grepl('^L6.', values)), 
          desc(grepl('SST|PV', values)),
          desc(grepl('LAMP5|NDNF|VIP|PAX', values)), 
          values) %>%
  pull(values) %>% unique()

cellTypes = cellTypes[c(1:5, 7,6, 8, 10:15, 9, 16:17, 19:20, 18, 21:37 )]

## set the colors for the joint cell types
exc_col = setNames(c(brewer.pal(10, 'Paired'), '#D56DA9'),
                   c('L2.CUX2.MEIS2', 'L4.5.TBX15', 'L4.ALPL', 'L3.CUX2.RORB',
                     'L6.ITGA8','L4.TYR', 'L6.NKD1', 'L6.SYT6', 'L5.PCP4', 
                     'L5.6.NR4A2', 'L5.POU3F1'))
inh_col = setNames(brewer.pal(8, 'Dark2'), cellTypes[21:28])
glia_col = setNames(brewer.pal(9, 'Pastel1'), cellTypes[29:37])

other_cols = cellTypes[! cellTypes %in% (c(exc_col, inh_col, glia_col) %>% names())]
other_cols = setNames(brewer.pal(length(other_cols), 'Set3'), other_cols)
cellTypes_cols = c(exc_col, other_cols, inh_col, glia_col) [cellTypes]
excTypes_cols = cellTypes_cols[1:20]
inhTypes_cols = cellTypes_cols[21:28]

###########################################
## 2) make the sankey diagram of allen to JH clusters, excitatory neurons
df = getCellColData(proj) %>% as.data.frame() %>% as_tibble() 

cM = df %>% filter( Celltype2 %in% names(excTypes_cols), 
                    allen_labels %in% names(excTypes_cols)) %>% 
  with(confusionMatrix(Celltype2, allen_labels))

rowmatch = match(names(excTypes_cols), rownames(cM))
rowmatch = rowmatch[complete.cases(rowmatch)]
colmatch = match(names(excTypes_cols), colnames(cM))
colmatch = colmatch[complete.cases(colmatch)]

cM = cM[rowmatch,  colmatch]
cM <- cM / Matrix::rowSums(cM)

pdf(here(PLOTDIR, 'macaque_PFC_scATAC_JH_AllenLabels.heatmap.exc.pdf'), onefile = F, 
    height = 120/in2mm, width = 120/in2mm)
pheatmap::pheatmap( 
  mat = as.matrix(cM), 
  color =colorRampPalette( c('white', 'black'))(100), 
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()



#########################################################################
## 3) make the sankey diagram of allen to JH clusters, excitatory neurons
cM = df %>% filter( Celltype2 %in% names(inhTypes_cols), 
                    allen_labels %in% names(inhTypes_cols)) %>% 
  with(confusionMatrix(Celltype2, allen_labels))

cM = cM[names(inhTypes_cols)[!grepl('PVALB|PAX6', names(inhTypes_cols))], 
        names(inhTypes_cols)[!grepl('\\.|NDNF', names(inhTypes_cols))]]
cM <- cM / Matrix::rowSums(cM)

pdf(here(PLOTDIR, 'macaque_PFC_scATAC_JH_AllenLabels.heatmap.inh.pdf'), onefile = F, 
    height = 120/in2mm, width = 120/in2mm)
pheatmap::pheatmap( 
  mat = as.matrix(cM), 
  color =colorRampPalette( c('white', 'black'))(100), 
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()




