ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
library(tidyverse)
library(arrow)
library(here)

DATADIR = 'data/tidy_data/celltype_specific_enhancers'
PLOTDIR = 'figures/exploratory/celltype_specific_enhancers'

#####################################################
## 0) read in all cell type comparisons
save_model_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_model_fn)

#####################################################
## 1) find all SVM model evaluations 
eval_fn = list.files(here(DATADIR, 'predictions'), pattern =  '_eval.feather', 
                     full.names = T, recursive = T)
eval_df = lapply(eval_fn, read_feather) %>% 
  data.table::rbindlist() %>% 
  filter(auROC > 0.51, accuracy > 0.51) %>%
  mutate(model_path = model,model = basename(model) %>% ss('_fold')) %>% 
  rename(eval = label) %>%
  arrange(model, desc(fhalf_score)) %>% 
  filter(!duplicated(model)) %>% 
  inner_join(x = cell_type_df) %>% 
  group_by(label) %>% 
  filter(n() > 2) %>% ungroup() %>%
  dplyr::select(-c(prefix, pred_pos:tp))

eval_df %>% write_tsv(here(DATADIR, 'tables/cell_type_bestSVMs_DLPFC.tsv'))
eval_df %>% saveRDS(here(DATADIR, 'rdas/cell_type_bestSVMs_DLPFC.rds'))

table(eval_df$label)
cell_type_df %>% filter(! model %in% eval_df$model )

pdf(here(PLOTDIR, 'plots', 'gkm_svm_validation_evaluations.pdf'), width = 8, height = 4)
ggplot(eval_df, aes(x = label, y= auPRC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_df, aes(x = label, y= auROC)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(eval_df, aes(x = label, y= f1_score)) + 
  geom_boxplot() + geom_jitter(pch = 21) +
  theme_bw() + ylim(c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


