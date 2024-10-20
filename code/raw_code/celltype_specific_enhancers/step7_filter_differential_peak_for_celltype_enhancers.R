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
  filter(length(unique(label)) <= 3) %>% ungroup() %>%
  split(f = .$model )

candidateList = lapply(diffPeakList, '[[', 'peakName')
candidateList = split(candidateList[models_df$model], models_df$label)
candidateList = lapply(candidateList, function(ll){
  data.frame(peak = Reduce('intersect', ll))
})

candidateEnhancers = rbindlist(candidateList, idcol = 'label') %>%
  group_by(peak) %>% filter(n()==1) %>% ungroup()
table(candidateEnhancers$label)


################################################
## 2) read in the predictions for all the models
eval_df = readRDS(here(DATADIR, 'rdas/cell_type_bestCNNandSVM_DLPFC.rds'))
pred_dir = here(DATADIR, 'predictions')
grep_pat = paste(eval_df$model_name, collapse = '|')

## read in predictions from the CNN, take the logit, akin to SVM score range
candidate_enh_pred_cnn_fn = list.files(pred_dir, full.names = T, recursive = T, 
                                       pattern = 'noncoding_peak.predictions.txt.gz') %>% str_subset(grep_pat)
names(candidate_enh_pred_cnn_fn) = candidate_enh_pred_cnn_fn %>% str_extract(grep_pat)
candidate_enh_pred_cnn_df = lapply(candidate_enh_pred_cnn_fn, fread) %>%
  lapply(function(df) {
    # some of the CNN predictions is 1, makes the logit +Inf
    df$y_pred_logit[df$y_pred_logit== Inf] = max(df$y_pred_logit[is.finite(df$y_pred_logit)])
    return(df[df$Name %in%candidateEnhancers$peak,])
  }) %>%
  rbindlist(idcol = 'model_name') %>% 
  dplyr::select(c('Name', 'y_pred_logit', 'model_name')) %>%
  dplyr::rename('score' = 'y_pred_logit', 'peak' = 'Name')

## add the per-model ranks
candidate_enh_pred_cnn_df = candidate_enh_pred_cnn_df %>% arrange(desc(score)) %>%
  group_by(model_name) %>% mutate(rank = order(score, decreasing = T)) %>% ungroup()


## read in the predictions from the SVM
candidate_enh_pred_svm_fn = list.files(pred_dir,full.names = T, recursive = T, 
                                       pattern = '_candidatesEnhancers.txt') %>% str_subset(grep_pat)
names(candidate_enh_pred_svm_fn) = candidate_enh_pred_svm_fn %>% str_extract(grep_pat)
candidate_enh_pred_svm_df = lapply(candidate_enh_pred_svm_fn, fread, col.names = c('peak', 'score')) %>%
  lapply(function(df)  {
    # some of the CNN predictions is 1, makes the logit +Inf
    df$score[df$score== Inf] = max(df$score[is.finite(df$score)])
    return(df[df$peak %in% candidateEnhancers$peak,])
  }) %>%  rbindlist(idcol = 'model_name')

## add the per-model ranks
candidate_enh_pred_svm_df = candidate_enh_pred_svm_df %>% arrange(desc(score)) %>%
  group_by(model_name) %>% mutate(rank = order(score, decreasing = T)) %>% ungroup()


## combine the models together, averaging scores
candidate_enh_pred_df = bind_rows(candidate_enh_pred_cnn_df, candidate_enh_pred_svm_df) %>%
  full_join(x = eval_df, y = .) %>% inner_join(y = candidateEnhancers) %>% 
  as.data.table() %>% dplyr::select(-c(bgd.group:model_type))

candidate_enh_pred_list = split(candidate_enh_pred_df, candidate_enh_pred_df$label)
candidate_enh_pred_wide = lapply(candidate_enh_pred_list, pivot_wider, names_from = 'model', 
              values_from = c('score', 'rank'), values_fn = mean) %>%
  lapply(function(df){
    offTarget = df %>% dplyr::select(starts_with('score')) %>%
      apply(1,FUN =  function(x) any(x < 0))
    df = df[!offTarget, ]
    avgScore = df %>% dplyr::select(starts_with('score')) %>%
      apply(1,FUN =  mean)
    avgRank = df %>% dplyr::select(starts_with('rank')) %>%
      apply(1,FUN =  mean) %>% order()
    df %>% mutate(AvgRank = avgRank, AvgScore = avgScore) %>%
      filter(AvgScore > 0) %>% arrange(AvgRank) %>%
      relocate(AvgRank, AvgScore, .after = peak) %>%
      dplyr::select(-starts_with('rank'))
  })
candidate_peaks = sapply(candidate_enh_pred_list, '[[', 'peak') %>% unlist()
candidate_enh_pred_wide %>% sapply(nrow)


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
save_motif_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.noncoding_peak.HOCOMOCO11_celltypeMotifs.se')
motif_unique_sums_se = readRDS(save_motif_fn)
motif_unique_sums_se = motif_unique_sums_se[candidate_peaks, ]
motif_counts = assays(motif_unique_sums_se)[[1]]


###############################################################################
## 5) add multimodal evidence of cell type specificity to ML model predictions
candidate_enh_list = lapply(candidate_enh_pred_wide, function(df){
  celltype = unique(df$label)
  df2 = df %>% 
    mutate(peak2gene = p2g_lookup[peak], peak2gene.cor = p2g_cor[peak],
           peak2gene.cor = ifelse(is.na(peak2gene.cor), 0, peak2gene.cor),
           markerGene = map_chr(strsplit(peak2gene, ','), function(x)
             x[min(which(unlist(x) %in% names(markerGenes_log2FC)), na.rm = T)]),
           markerGene.log2fc = markerGenes_log2FC[markerGene],
           markerGene.log2fc = ifelse(is.infinite(markerGene.log2fc), 0, markerGene.log2fc),
           MotifZscore = motif_counts[peak,celltype])
  tmp = apply(df2%>% dplyr::select(AvgScore:MotifZscore, -c(peak2gene, markerGene)) , 1, gm)
  df2 = df2 %>% mutate(compositeScore = tmp) %>%
    arrange(desc(compositeScore), AvgRank,  is.na(peak2gene), peak2gene.cor, !is.na(markerGene.log2fc), 
            markerGene.log2fc, is.na(MotifZscore)) %>%
    relocate(compositeScore, peak2gene:MotifZscore, .before = AvgScore) %>%
    mutate(AvgRank = order(-compositeScore))
  df2 = df2 %>% mutate(label =ifelse(is.na(peak2gene),  paste(label, AvgRank, sep = "_"), 
                              paste(label, AvgRank, peak2gene, sep = "_"))) 
  return(df2)
})


save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
saveRDS(candidate_enh_list, save_top_enh_fn)

save_top_enh_excel = here(DATADIR, 'tables', 'rheMac10_DLPFC.candidate_celltype_enhancers.xlsx')
writexl::write_xlsx(candidate_enh_list, save_top_enh_excel)

head(candidate_enh_list[[1]] %>% data.frame())
sapply(candidate_enh_list, nrow)







