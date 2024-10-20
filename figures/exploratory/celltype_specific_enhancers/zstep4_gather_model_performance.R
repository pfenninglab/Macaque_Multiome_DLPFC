ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(tidyverse)
library(readxl)
library(here)

addArchRThreads(threads = 4)

######################
## load ArchR projects
DATADIR=

df = here('data/tidy_data/celltype_specific_enhancers/tables', 
          'cell_type_bestCNNandSVM_DLPFC.tsv') %>% read_tsv()

df %>% group_by(label) %>% 
  summarise(mean_auROC = mean(auROC), 
            mean_auPRC = mean(auPRC)) %>% 
  summarise(range(mean_auROC), 
           range(mean_auPRC))
  

with(df, range(auROC))
with(df, range(auPRC))


df2 = here('figures/exploratory/celltype_specific_enhancers/tables',
           'candidateCelltypeEnhancerSummaryTable_multiomeATAC_DLPFC_20220418.xlsx') %>% read_xlsx()

df2 %>% dplyr::select(numCandidateEnhancers, numMLselectedEnhancers) %>% 
  apply(2, mean, na.rm = T)

df2 %>% dplyr::select(numCandidateEnhancers, numMLselectedEnhancers) %>% 
  apply(2, sd, na.rm = T)/sqrt(nrow(df2))


