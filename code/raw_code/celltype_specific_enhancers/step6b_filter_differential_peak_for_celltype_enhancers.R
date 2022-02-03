ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(ArchR)

DATADIR = 'data/tidy_data/celltype_specific_enhancers'
PLOTDIR = 'figures/exploratory/celltype_specific_enhancers'

################################################
## 1) read in the differential peaks
save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
models_df = readRDS(save_models_fn)

save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models_diffPeakList_DLPFC.rds')
diffPeakList = readRDS(save_diffPeaks_fn) %>% lapply(function(gr) gr[gr$Log2FC > 0])
candidateList = lapply(diffPeakList, names)
candidateList = split(candidateList[models_df$model], models_df$label)
candidateList = lapply(candidateList, function(ll){
  ll = ll[lengths(ll) > 1000]
  data.frame(peak = Reduce('intersect', ll))
})

candidateEnhancers = rbindlist(candidateList, idcol = 'label') %>%
  group_by(peak) %>% filter(n()==1) %>% ungroup()

################################################
## 2) read in the predictions for all the models
eval_df = readRDS(here(DATADIR, 'rdas/cell_type_bestCNNandSVM_DLPFC.rds'))
pred_dir = here(DATADIR, 'predictions')
grep_pat = paste(eval_df$model_name, collapse = '|')

## read in predictions from the SVM
candidate_enh_pred_cnn_fn = list.files(pred_dir, full.names = T, recursive = T, 
                                       pattern = 'noncoding_peak.predictions.txt.gz') %>% str_subset(grep_pat)
names(candidate_enh_pred_cnn_fn) = candidate_enh_pred_cnn_fn %>% str_extract(grep_pat)
candidate_enh_pred_cnn_df = lapply(candidate_enh_pred_cnn_fn, fread) %>%
  lapply(function(df) df[df$Name %in%candidateEnhancers$peak,]) %>%
  rbindlist(idcol = 'model_name') %>% 
  dplyr::select(c('Name', 'y_pred_logit', 'model_name')) %>%
  dplyr::rename('score' = 'y_pred_logit', 'peak' = 'Name')

## read in the predictions from the CNN, take the logit, akin to SVM score range
candidate_enh_pred_svm_fn = list.files(pred_dir,full.names = T, recursive = T, 
                                       pattern = '_candidatesEnhancers.txt') %>% str_subset(grep_pat)
names(candidate_enh_pred_svm_fn) = candidate_enh_pred_svm_fn %>% str_extract(grep_pat)
candidate_enh_pred_svm_df = lapply(candidate_enh_pred_svm_fn, fread, col.names = c('peak', 'score')) %>%
  lapply(function(df) df[df$peak %in%candidateEnhancers$peak,]) %>%
  rbindlist(idcol = 'model_name')

candidate_enh_pred_df = bind_rows(candidate_enh_pred_cnn_df, candidate_enh_pred_svm_df) %>%
  full_join(x = eval_df, y = .) %>% inner_join(y = candidateEnhancers) %>% 
  as.data.table() %>% dplyr::select(-c(bgd.group:model_type))

candidate_enh_pred_list = split(candidate_enh_pred_df, candidate_enh_pred_df$label)
candidate_enh_pred_wide = lapply(candidate_enh_pred_list, pivot_wider, names_from = 'model', 
              values_from = 'score', values_fn = mean) %>%
  lapply(function(df){
    df %>% mutate(AvgScore = rowMeans(across(where(is.numeric)))) %>%
      filter(AvgScore > 0) %>% arrange(desc(AvgScore)) %>%
      mutate(label = paste(label, seq(n()), sep = '_')) %>%
      relocate(AvgScore, .after = peak)
  })

