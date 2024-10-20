ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(here)

source(here('code/final_code/hal_scripts/narrowPeakFunctions.R'))
source(here('code/final_code/hal_scripts/gen_enh_ortholog_sets.R'))

######################
## load ArchR projects
DATADIR='data/tidy_data'
filter_xlsx_path = here(DATADIR,'celltype_specific_enhancers/tables',
                        'rheMac10_DLPFC_L3CUX3RORB_L5.POU3F1.enhancers_for_screening_20220610.xlsx')
enhancer_list = filter_xlsx_path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = filter_xlsx_path)

cols = Reduce('intersect', lapply(enhancer_list, names))

enhancer_df = lapply(enhancer_list, function(df) 
  df %>% dplyr::select(all_of(cols))
) %>% rbindlist(idcol = 'celltype')

enhancer_gr = enhancer_df %>% mutate(
  name = peak,
  peak = gsub('rheMac10:|:250$', '', peak)
) %>% dplyr::select(label, peak)%>% deframe() %>%
  GRanges()
enhancer_gr$name = enhancer_df$peak

enhancer_fa = here(DATADIR,'celltype_specific_enhancers/fasta','rheMac10_DLPFC_L3CUX3RORB_L5.POU3F1.candidate_celltype_enhancers.toOrder.20220613.fa')
enhancer_seq =  writeGRangesToFasta(gr = enhancer_gr, file = enhancer_fa, genome = 'rheMac10')


enhancer_df = enhancer_df %>%
  mutate(enhancer_seq = data.frame(enhancer_seq)$enhancer_seq)
filter_xlsx_path2 = here(DATADIR,'celltype_specific_enhancers/tables',
                         'rheMac10_DLPFC_L3CUX3RORB_L5.POU3F1.candidate_celltype_enhancers.toOrder.20220613.xlsx')
split(enhancer_df, enhancer_df$celltype) %>%
  writexl::write_xlsx(filter_xlsx_path2)

filter_xlsx_path3 = here(DATADIR,'celltype_specific_enhancers/tables',
                         'rheMac10_DLPFC_L3CUX3RORB_L5.POU3F1.candidate_celltype_enhancers.toOrderLong.20220613.xlsx')
writexl::write_xlsx(enhancer_df, filter_xlsx_path3)


