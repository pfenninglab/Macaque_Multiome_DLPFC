suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)
library(data.table)
library(ggtree)
library(ggtreeExtra)
library(tidytree)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }


########################################
## load genomeAnnotation, geneAnnotation
PLOTDIR='figures/exploratory/cell_type_hierarchy/plots/'
DATADIR = 'data/tidy_data/celltype_specific_enhancers'
DATADIR2 = 'figures/exploratory/celltype_specific_enhancers/tables/'

save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
models_df = readRDS(save_models_fn)

###############################
## compute differential peaks
alpha = 0.05
save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models.all_diffPeaks.DESeq2.rds')
diffPeakList =  readRDS(save_diffPeaks_fn) %>%
  lapply(rownames_to_column, "peakName") %>% rbindlist(idcol = 'model') %>%
  filter(padj < alpha, log2FoldChange > 0) %>% 
  full_join(models_df) %>% group_by(peakName) %>%
  filter(length(unique(label)) <= 3) %>% ungroup() %>%
  split(f = .$model )

candidateList = lapply(diffPeakList, '[[', 'peakName')
candidateList = split(candidateList[models_df$model], models_df$label)
candidateList = lapply(candidateList, function(ll){
  data.frame(peak = Reduce('union', ll))
})

candidateEnhancersList = rbindlist(candidateList, idcol = 'label') %>%
  group_by(peak) %>% filter(n()==1) %>% ungroup() %>%
  split(f = .$label)

tmp = sapply(candidateEnhancersList, nrow)
mean(tmp)
sd(tmp)/sqrt(length(candidateEnhancersList))

######################
## get ArchR project
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
proj = loadArchRProject(path = PROJDIR2)
celltypesList = split(proj$Celltype2, proj$Celltype2)
celltypesList = celltypesList[names(candidateEnhancersList)]


## read in table of prioritized enhancers
save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
candidate_enh_list = readRDS(save_top_enh_fn)
numCandidatePeaks = sapply(candidate_enh_list, nrow)
numCandidatePeaks = numCandidatePeaks[names(candidateEnhancersList)]

df = data.frame(celltype = names(candidateEnhancersList),
                numCells = lengths(celltypesList), 
                numCandidateEnhancers = sapply(candidateEnhancersList, nrow),
                numMLselectedEnhancers = numCandidatePeaks)

save_table_fn = here(DATADIR2,'candidateCelltypeEnhancerSummaryTable_multiomeATAC_DLPFC_20220418.xlsx')
writexl::write_xlsx(df, save_table_fn)

