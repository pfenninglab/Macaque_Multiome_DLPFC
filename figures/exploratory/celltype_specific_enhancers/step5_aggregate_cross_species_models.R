library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(data.table)
library(readxl)
library(here)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

## read in the fasta sequences 
fa_500 = import(here('data/tidy_data/celltype_specific_enhancers/fasta/rheMac10_DLPFC.noncoding_peak.500.fa'), format = 'fasta')

## import the initial ranked enhancers
path = here('data/tidy_data/celltype_specific_enhancers/tables/rheMac10_DLPFC.candidate_celltype_enhancers.xlsx')
sheets <- excel_sheets(path)
enh_list <- lapply(set_names(sheets), function(x) read_excel(path, sheet = x))

path2 = here('data/tidy_data/celltype_specific_enhancers/tables/rheMac10_DLPFC.candidate_celltype_enhancers_withoutML.xlsx')
sheets2 <- excel_sheets(path2)
enh_list2 <- lapply(set_names(sheets2), function(x) read_excel(path2, sheet = x))

shared_columns = Reduce(intersect, lapply(c(enh_list, enh_list2), names))
enh_df = lapply(c(enh_list, enh_list2), function(df) 
  df %>% dplyr::select(all_of(shared_columns))) %>% 
  rbindlist(idcol = 'celltype')

## import the tested enhancers
path2 = here('data/raw_data/tables/Enhancer Sequence and Summary Data 20231805.xlsx')
test_df <- read_excel(path2, sheet = 'PFC 1st Enhancers')



## import the predictions from the multi-species model
fn = list.files(path = '/projects/pfenninggroup/machineLearningForComputationalBiology/Cortex_Cell-TACIT/data/tidy_data/Cortex-SNAIL/predictions/monkey_snail_v1',
                full.names = T)
names(fn) = basename(fn)

pred = lapply(fn, fread, header = F, col.names = 'pred') %>% 
  lapply(function(df){
    df %>% mutate(peak = names(fa_500)) %>% filter(peak %in% enh_df$peak)
  }) %>% rbindlist(idcol = 'model')

pred2 = pred %>% 
  mutate(fold = ss(model, '_', 6), model = ss(model, '_', 5), 
         celltype2 = ss(model, 'vs', 1), 
         background = ss(model, 'vs', 2),
         background = paste(celltype2, 'vs', background),
         pred_logit = log(pred/(1-pred))) %>% 
  group_by(celltype2, background, peak) %>% 
  summarise(pred_logit = mean(pred_logit),
            pred = 1/(1 + exp(-pred_logit))) %>% 
  inner_join(test_df) %>% 
  mutate(celltype = ss(celltype, '\\.') %>% paste0('PN')) %>% 
  filter(background!= 'L2.3.IT vs ITexc')

with(pred2, table(celltype, background))

pred3 = pred2 %>% group_by(peak, celltype2) %>% 
  mutate(background = 'Avg Prediction', pred_logit = mean(pred_logit)) %>% 
  ungroup() %>% arrange(background) %>% 
  distinct(peak, celltype2, pred_logit, .keep_all = T)


for(model in unique(pred3$celltype2)){
  fn = here('figures/exploratory/celltype_specific_enhancers', 'plots')
  fn = file.path(fn, paste0('cross_species_model.', model,'.pdf'))
  p1 = ggplot(bind_rows(pred2, pred3) %>% filter(celltype2 %in% model), 
              aes(x = celltype, y = pred_logit)) +
    geom_hline(yintercept = 0, color = 'black') +
    geom_boxplot(aes(fill = celltype), outlier.shape = NA, alpha = 0.7) + 
    geom_text(position=position_jitter(), size = 2,
              aes(label = barcode, color = celltype)) +
    scale_fill_brewer(palette = 'Set1', guide = 'none') + 
    scale_color_brewer(palette = 'Set1', guide = 'none') + 
    facet_grid(~background, space = 'free', scales = 'free_x') + 
    theme_bw(base_size = 8) + 
    ylab('Prediction Score') + xlab('Tested Enhancers') 
  
  pdf(fn, width = 4.25, height = 1.5)
  print(p1)  
  dev.off()
}



