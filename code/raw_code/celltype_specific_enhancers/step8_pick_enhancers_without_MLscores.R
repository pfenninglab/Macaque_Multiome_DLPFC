ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(ArchR)
library(Rfast)

DATADIR = 'data/tidy_data/celltype_specific_enhancers'
PLOTDIR = 'figures/exploratory/celltype_specific_enhancers'

## function to compute geometric mean
gm = function(...){ expm1(mean(log1p(...), na.rm = T))}

################################################
## 1) read in the differential peaks
save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
models_df = readRDS(save_models_fn)

save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models.all_diffPeaks.DESeq2.rds')

alpha = 0.05
diffPeakList =  readRDS(save_diffPeaks_fn) %>%
  lapply(rownames_to_column, "peakName") %>% rbindlist(idcol = 'model') %>%
  filter(padj < alpha, log2FoldChange > 0) %>% 
  full_join(models_df) %>% group_by(peakName) %>%
  filter(length(unique(label)) <= 6) %>% ungroup() %>%
  split(f = .$model )

candidateList = lapply(diffPeakList, '[[', 'peakName')
candidateList = split(candidateList[models_df$model], models_df$label)
candidateList = lapply(candidateList, function(ll){
  data.frame(peak = Reduce('intersect', ll))
})

candidate_peaks = lapply(diffPeakList, '[[', 'peakName') %>% unlist()

save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
candidate_enh_list_withML = readRDS(save_top_enh_fn)
candidate_enh_list_withML = candidate_enh_list_withML[sapply(candidate_enh_list_withML, nrow)> 0]

candidateList = candidateList[names(candidateList) %ni% names(candidate_enh_list_withML)]
candidateEnhancers = rbindlist(candidateList, idcol = 'label') %>%
  group_by(peak) %>% filter(n()<3) %>% ungroup()
table(candidateEnhancers$label)



################################################
## 3) retrieve Peak2GeneLinks and MarkerGenes

## cell type-specific marker genes
save_DEGs_fn = here(DATADIR, 'rdas/cell_type_models.all_DEGs.DESeq2.rds')
de.markers.list = readRDS(save_DEGs_fn)
de.markers.list2 = de.markers.list %>% lapply(as.data.frame) %>%
  lapply(rownames_to_column, "Gene.Symbol") %>% rbindlist(idcol = 'model') %>%
  mutate(Gene.Symbol = as.character(Gene.Symbol)) %>%
  filter(padj < alpha, log2FoldChange > 0) %>%
  full_join(x = models_df) %>%
  group_by(label, Gene.Symbol) %>% top_n(1, -log10(padj)) %>% ungroup() %>% 
  group_by(Gene.Symbol) %>% filter(length(unique(label)) <=3) %>%  ungroup() %>%
  dplyr::select(-c('model', 'bgd.group', 'bgd.labels')) %>%
  arrange(Gene.Symbol) %>% split(f = .$label)

de.markers = de.markers.list2 %>% rbindlist(idcol = 'celltype') 
markerGenes = de.markers %>% dplyr::select(celltype, Gene.Symbol) %>% deframe()
markerGenes2 = de.markers %>% dplyr::select(Gene.Symbol, celltype) %>% deframe()
markerGenes_log2FC = setNames(de.markers$log2FoldChange, de.markers$Gene.Symbol)
markerGenes_log2FC = sort(markerGenes_log2FC, decreasing= TRUE)

## load the ArchR project w/ the samples and get peak counts matrix
proj = loadArchRProject(path = file.path('data/tidy_data/ArchRProjects',
                                         'ArchR_DLPFC_multiomeATAC'))

p2g <- getPeak2GeneLinks(proj, corCutOff = 0.1, resolution = 1, returnLoops = FALSE)
metadata(p2g)$peakSet$name = with(data.frame(metadata(p2g)$peakSet), 
                                  paste0('rheMac10:',seqnames, ':', start, '-', end, ':', 250))
p2g$gene =metadata(p2g)$geneSet$name[p2g$idxRNA]
p2g$peak =metadata(p2g)$peakSet$name[p2g$idxATAC]

## subset to just marker genes, add the cell type each gene is a marker for
p2g = p2g[p2g$peak %in% candidate_peaks, ]
p2g = p2g[order(log(p2g$FDR)),]
p2g$celltype = markerGenes2[p2g$gene]
table(p2g$celltype)
p2g_lookup = split(p2g$gene, p2g$peak) %>% sapply(function(x) paste(sort(x), collapse = ','))
p2g_cor = setNames(p2g$Correlation, p2g$peak)


########################################
## 4) retrieve cell type unique motifs
save_motif_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.noncoding_peak.HOCOMOCO11_celltypeMotifs.se.rds')
motif_unique_sums_se = readRDS(save_motif_fn)
motif_unique_sums_se = motif_unique_sums_se[candidate_peaks, ]
motif_counts = assays(motif_unique_sums_se)[[1]]


###############################################################################
## 5) add multimodal evidence of cell type specificity to ML model predictions
candidate_enh_list = lapply(split(candidateEnhancers, candidateEnhancers$label), function(df){
  celltype = unique(df$label)
  df2 = df %>% 
    mutate(peak2gene = p2g_lookup[peak], peak2gene.cor = p2g_cor[peak],
           peak2gene.cor = ifelse(is.na(peak2gene.cor), 0, peak2gene.cor),
           markerGene = map_chr(strsplit(peak2gene, ','), function(x)
             x[min(which(unlist(x) %in% names(markerGenes_log2FC)), na.rm = T)]),
           markerGene.log2fc = markerGenes_log2FC[markerGene],
           markerGene.log2fc = ifelse(is.infinite(markerGene.log2fc), 0, markerGene.log2fc),
           MotifZscore = motif_counts[peak,celltype])
  tmp = apply(df2%>% dplyr::select(peak2gene.cor:MotifZscore, -c(peak2gene, markerGene)) , 1, function(x) mean(x, na.rm = T) )
  df2 = df2 %>% mutate(compositeScore = tmp) %>%
    arrange(desc(compositeScore),  is.na(peak2gene), peak2gene.cor, !is.na(markerGene.log2fc), 
            markerGene.log2fc, is.na(MotifZscore)) %>%
    relocate(compositeScore, peak2gene:MotifZscore, .before = MotifZscore) %>%
    filter(!is.na(markerGene)) %>%
    mutate(AvgRank = order(-compositeScore))
  df2 = df2 %>% mutate(label =ifelse(is.na(peak2gene),  paste(label, AvgRank, sep = "_"), 
                                     paste(label, AvgRank, peak2gene, sep = "_"))) 
  return(df2)
})
candidate_enh_list = candidate_enh_list[sapply(candidate_enh_list, nrow)> 0]
sapply(candidate_enh_list, nrow)

save_top_enh_excel = here(DATADIR, 'tables', 'rheMac10_DLPFC.candidate_celltype_enhancers_withoutML.xlsx')
writexl::write_xlsx(candidate_enh_list, save_top_enh_excel)

save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers_withoutML.rds')
saveRDS(candidate_enh_list, save_top_enh_fn)






