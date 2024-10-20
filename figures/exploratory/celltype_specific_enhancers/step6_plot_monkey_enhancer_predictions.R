library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(data.table)
library(readxl)
library(here)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

## import the tested enhancers
path2 = here('data/raw_data/tables/Enhancer Sequence and Summary Data 20231805.xlsx')
test_df <- read_excel(path2, sheet = 'PFC 1st Enhancers') %>% 
  dplyr::select(peak, barcode)

## import the initial ranked enhancers
path = here('data/tidy_data/celltype_specific_enhancers/tables/rheMac10_DLPFC.candidate_celltype_enhancers.xlsx')
sheets <- excel_sheets(path)
enh_list <- lapply(set_names(sheets), function(x) read_excel(path, sheet = x))

shared_columns = Reduce(intersect, lapply(enh_list, names))
enh_df = lapply(enh_list, function(df) 
  df %>% pivot_longer(cols = starts_with('score'), 
                      names_to = 'background', values_to = 'pred_logit')
  ) %>% 
  rbindlist(idcol = 'celltype', fill=TRUE)

enh_df = enh_df %>% 
  mutate(background = str_replace_all(background, 'score_', ''), 
         background = str_replace_all(background, 'L3.CUX2.RORB', 'L3PN '), 
         background = str_replace_all(background, 'vs', 'vs\n'), 
         background = str_replace_all(background, '_|\\.', ' '),
         background = str_replace_all(background, 'Lower |Upper ', ''),
         celltype = str_replace_all(celltype, 'L3.CUX2.RORB', 'L3PN')
  ) %>% inner_join(test_df)

enh_df2 = enh_df %>% group_by(peak, celltype) %>% 
  mutate(background = 'Avg\nPrediction', pred_logit = mean(pred_logit)) %>% 
  ungroup() %>% arrange(background) %>% 
  distinct(peak, celltype, pred_logit, .keep_all = T)


for(model in unique(enh_df$celltype)){
  fn = here('figures/exploratory/celltype_specific_enhancers', 'plots')
  fn = file.path(fn, paste0('macaque_only_model.', model,'.pdf'))
  p1 = ggplot(bind_rows(enh_df, enh_df2) %>% filter(celltype %in% model), 
              aes(x = celltype, y = pred_logit)) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_boxplot(aes(fill = celltype), outlier.shape = NA, alpha = 0.7) + 
    geom_text(position=position_jitter(), size = 2,
              aes(label = barcode)) +
    scale_fill_brewer(palette = 'Set1', guide = 'none') + 
    scale_color_brewer(palette = 'Set1', guide = 'none') + 
    facet_grid(~background, space = 'free', scales = 'free_x') + 
    theme_bw(base_size = 8) + 
    ylab('Prediction Score') + xlab('Tested Enhancers') 
  
  pdf(fn, width = 4.25, height = 1.5)
  print(p1)  
  dev.off()
}



