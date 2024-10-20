ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
library(tidyverse)
library(arrow)
library(here)
library(data.table)

DATADIR = 'data/tidy_data/celltype_specific_enhancers'
PLOTDIR = 'figures/exploratory/celltype_specific_enhancers'

#####################################################
## 0) read in all cell type comparisons
save_model_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_model_fn)

#####################################################
## 1) find all SVM model evaluations 
eval_svm_fn = list.files(here(DATADIR, 'predictions'), pattern =  '_eval.feather', 
                     full.names = T, recursive = T) %>% str_subset('gkmpredict')

eval_svm_df = lapply(eval_svm_fn, read_feather) %>% 
  data.table::rbindlist() %>% 
  filter(auROC > .65, accuracy > .65) %>%
  mutate(model_path = model,model = basename(model) %>% ss('_fold'),
         model_type = 'SVM') %>% 
  rename(eval = label) %>%
  arrange(model, desc(fhalf_score)) %>% 
  filter(!duplicated(model)) %>% 
  inner_join(x = cell_type_df) %>% 
  group_by(label) %>% 
  filter(n() > 2) %>% ungroup() %>%
  dplyr::select(-c(prefix, pred_pos:tp))

eval_svm_df %>% write_tsv(here(DATADIR, 'tables/cell_type_bestSVMs_DLPFC.tsv'))
eval_svm_df %>% saveRDS(here(DATADIR, 'rdas/cell_type_bestSVMs_DLPFC.rds'))

pdf(here(PLOTDIR, 'plots', 'gkm_svm_validation_evaluations.pdf'), width = 8, height = 4)
ggplot(eval_svm_df, aes(x = label, y= auPRC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_svm_df, aes(x = label, y= auROC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_svm_df, aes(x = label, y= f1_score)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


######################################
## 2) find all CNN model evaluations 
eval_cnn_fn = list.files(here(DATADIR, 'predictions'), pattern =  'performance.feather', 
                         full.names = T, recursive = T) %>% str_subset('cnn') %>%
  str_subset('.h5_', negate = T)

eval_cnn_df = lapply(eval_cnn_fn, read_feather) %>% 
  data.table::rbindlist(fill = TRUE) %>% 
  filter(auROC > .65, accuracy > .65) %>%
  mutate(model_path = model,model = basename(model) %>% ss('_b64e30.h5'),
         model_type = 'CNN', eval = 'valid') %>% 
  arrange(model, desc(fhalf_score)) %>% 
  filter(!duplicated(model)) %>% 
  inner_join(x = cell_type_df) %>% 
  group_by(label) %>% 
  filter(n() > 2) %>% ungroup() %>%
  dplyr::select(all_of(names(eval_svm_df)))


eval_cnn_df %>% write_tsv(here(DATADIR, 'tables/cell_type_bestCNNs_DLPFC.tsv'))
eval_cnn_df %>% saveRDS(here(DATADIR, 'rdas/cell_type_bestCNNs_DLPFC.rds'))

pdf(here(PLOTDIR, 'plots', 'cnn_validation_evaluations.pdf'), width = 8, height = 4)
ggplot(eval_cnn_df, aes(x = label, y= auPRC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_cnn_df, aes(x = label, y= auROC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_cnn_df, aes(x = label, y= f1_score)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


######################################
## 3) join the best models together 
eval_df = bind_rows(eval_svm_df, eval_cnn_df)
eval_df$model_name = basename(eval_df$model_name) %>% ss('.model.txt|.h5')
eval_df %>% write_tsv(here(DATADIR, 'tables/cell_type_bestCNNandSVM_DLPFC.tsv'))
eval_df %>% saveRDS(here(DATADIR, 'rdas/cell_type_bestCNNandSVM_DLPFC.rds'))

eval_col = c( "accuracy", "auROC", "auPRC", "f1_score", 'fhalf_score')
eval_pivot_df = eval_df %>% dplyr::select(-c(model, model_path, model_name)) %>%
  pivot_longer(cols = all_of(eval_col), names_to = 'metric', values_to = 'value') %>%
  pivot_wider(names_from = 'model_type', values_from = 'value', values_fill = 0)

pdf(here(PLOTDIR, 'plots', 'compare_CNNvsSVM_validation_evaluations.pdf'), width = 8, height = 4)
ggplot(eval_pivot_df, aes(x = SVM, y= CNN)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red')+
  geom_point(pch = 21, aes(fill = label)) +w
  theme_bw() +
  facet_wrap(~metric, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

  
  





