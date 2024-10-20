suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(rtracklayer)
library(future)
library(harmony)
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(GenomicRanges))

## for plotting 
library(ggsankey)
library(RColorBrewer)

DATADIR='data/tidy_data/celltype_specific_enhancers'
PLOTDIR='figures/exploratory/jing_trackPlots/plots'
in2mm = 25.4

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

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



save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
models_df = readRDS(save_models_fn)


alpha = 0.05
save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models.all_diffPeaks.DESeq2.rds')
diffPeakList =  readRDS(save_diffPeaks_fn) %>%
  lapply(rownames_to_column, "peakName") %>% rbindlist(idcol = 'model') %>%
  filter(padj < alpha, log2FoldChange > 0) %>% 
  full_join(models_df) %>% ungroup() %>% split(f = .$model )
  candidateList = lapply(diffPeakList, '[[', 'peakName')
  candidateList = split(candidateList[models_df$model], models_df$label)

  candidateList = lapply(candidateList, function(ll){
  data.frame(peak = Reduce('intersect', ll))
})

candidateList = candidateList[sapply(candidateList, nrow) > 10]
candidateEnhancers = rbindlist(candidateList, idcol = 'label') %>%
  group_by(peak) %>% filter(n()==1) %>% ungroup()
table(candidateEnhancers$label)


save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
markerPeaks = readRDS(save_top_enh_fn) %>% rbindlist(fill = TRUE) %>% 
  dplyr::select(c(label, peak, AvgScore, AvgRank))%>%  
  mutate(label= ss(label, '_'))
table(markerPeaks$label)
summary(markerPeaks$AvgScore)

table(candidateEnhancers$peak %in% markerPeaks$peak)
table(markerPeaks$peak %in% candidateEnhancers$peak)

df = left_join(candidateEnhancers, markerPeaks) %>% 
  mutate(AvgScore = ifelse(is.na(AvgScore), 0, AvgScore)) %>% 
  group_by(label) %>% 
  arrange(is.na(AvgRank), desc(AvgScore)) %>% 
  mutate(AvgRank = seq(n())) %>% 
  ungroup() %>% filter(grepl('^L[1-6]', label)) %>% 
  filter(label %in% names(excTypes_cols)) %>% 
  filter(label %in% unique(markerPeaks$label)) %>% 
  mutate(group = ifelse(AvgScore ==0, 'NA', label), 
         label = factor(label, names(excTypes_cols)), 
         label = droplevels(label))

pdf(here(PLOTDIR, 'macaque_PFC_scATAC_candidateEnhancer.rankplot.pdf'), onefile = F, 
    height = 9, width = 30)
ggplot(df, aes(x = AvgRank, y = AvgScore)) + 
  geom_line(size = 2, color = 'black') +
  geom_point(pch = 20, size = 10 , aes(color = group)) +
  facet_wrap(~label, scales = 'free_x', nrow = 2) +
  theme_classic(base_size = 26) + scale_x_continuous(trans='log10')  +
  ylab("ML Enhancer Ranking Score") + xlab("Candidate enhancer rank") +
  scale_color_manual(values=excTypes_cols[levels(df$label)])+
  theme(legend.position = 'none')
dev.off()


