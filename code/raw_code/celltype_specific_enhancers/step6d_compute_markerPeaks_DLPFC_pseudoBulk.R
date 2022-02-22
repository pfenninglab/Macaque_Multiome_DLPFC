# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(ArchR)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(here)
library(future)
library(data.table)

## for pseudo-bulk analyses
library(scuttle)
library(DESeq2)
library(limma)
library(BiocParallel)
library(IHW)


DATADIR='data/tidy_data/celltype_specific_enhancers'
source(here('code/final_code/hal_scripts/narrowPeakFunctions.R'))
source(here('code/final_code/hal_scripts/gen_enh_ortholog_sets.R'))

## load in the sample sheet
pd = read_csv(here('data/raw_data/tables/SampleSheet_RNA_DLPFC.csv')) %>%
  rename_with(make.names) %>% 
  rename_all(funs(stringr::str_replace_all(., '^X\\.\\.', ''))) %>%
  mutate(Animal = paste0('Monkey_',stringr::str_sub(Nickname, start = 1, end = 1)),
         Sample = paste0(Nickname, '_', Region, '-',Batch)) %>%
  relocate(Sample, Region, Animal, Batch, .before = everything())


#######################################################
# 0) Seurat uses the future package for parallelization
register(MulticoreParam(24))
addArchRThreads(threads = 24)


###############################################
## 1) load SeuratObjects of labeled all nuclei
pbulk_fn = file.path('data/tidy_data/rdas', 'pseudobulk_candidate_enhancers_snATAC_DLPFC.sce.rds')
if(file.exists(pbulk_fn)) {
  sce_agg = readRDS(pbulk_fn)
} else {
  proj = loadArchRProject(path = file.path('data/tidy_data/ArchRProjects',
                                           'ArchR_DLPFC_scATAC'))
  proj$tmp = paste(proj$Sample, proj$Celltype2, sep = '#')
  
  sce_agg <- getGroupSE( proj, useMatrix = 'PeakMatrix', groupBy = "tmp",
                         divideN = FALSE, scaleTo = NULL, verbose = TRUE)
  assays(sce_agg)$PeakMatrix = round(assays(sce_agg)$PeakMatrix)
  colData(sce_agg)$Celltype2 = ss(rownames(colData(sce_agg)), '#', 2)
  colData(sce_agg)$Sample = ss(rownames(colData(sce_agg)), '#', 1)
  with(colData(sce_agg), table(Celltype2))
  
  ## subset to just candidate enhancer peaks
  peak_gr = getPeakSet(proj)
  peak_gr = peak_gr %>% addSummitCenter() %>% nameNarrowPeakRanges(genome = 'rheMac10')
  indKeep = which(peak_gr$peakType %in% c('Distal', 'Intronic') & abs(peak_gr$distToTSS) > 2000)
  peak_gr = peak_gr[indKeep]
  rownames(sce_agg) = with(rowData(sce_agg) %>% data.frame, 
                           paste0('rheMac10:', seqnames, ':', start, '-', end, ':250'))
  rowData(sce_agg) = data.frame(name = rownames(sce_agg))
  sce_agg = sce_agg[mcols(peak_gr)$name,]
  saveRDS(sce_agg, pbulk_fn)
}


## calculate the detection rate of genes in this cell type
y <- assays(sce_agg)$PeakMatrix
sce_agg$cdr <- scale(colMeans(y > 0))
sce_agg = sce_agg
dds <- DESeqDataSet(sce_agg, design = ~ 0 + Celltype2 + Sample + cdr)
dds = DESeq(dds, parallel = TRUE)


###############################################
## 2) read in the differential peak comparisons
save_models_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_models_fn)
model_name = cell_type_df$model
indlist = setNames(seq_along(model_name), model_name)
fgd_labels = setNames(cell_type_df$label, model_name)
bgd_labels = setNames(cell_type_df$bgd.labels %>% sapply(strsplit,';'), model_name)
bgd_labels = lapply(bgd_labels, function(b) b[b %in% levels(dds$Celltype2)])

contrastsList = mapply(function(f, g){
  paste(f , '-(', paste(g, collapse = '+'), ')/', length(g) )
}, f = fgd_labels, g = bgd_labels)

##############################################################
## 3) compute differential marker genes across all comparisons
de.markers.list = lapply(contrastsList, function(cont){
  print(paste('Getting diffPeaks for', cont))
  contrasts =  c(makeContrasts(contrasts = cont, levels = levels(dds$Celltype2)), rep(0,5))
  results(dds, filterFun=ihw, contrast =contrasts)
})
names(de.markers.list) = model_name
save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models.all_diffPeaks.DESeq2.rds')
saveRDS(de.markers.list %>% lapply(as.data.frame), save_diffPeaks_fn)

##############################################################
## 4) filter to significant diffPeaks, 
alpha = 0.05
de.markers.list2 = de.markers.list %>% lapply(as.data.frame) %>%
  lapply(rownames_to_column, "peakName") %>% rbindlist(idcol = 'model') %>%
  mutate(peakName = as.character(peakName)) %>%
  filter(padj < alpha, log2FoldChange > 2) %>%
  full_join(x = cell_type_df) %>%
  group_by(label, peakName) %>% top_n(1, -log10(padj)) %>% ungroup() %>% 
  ## diff peak unique to this cell type
  group_by(peakName) %>% filter(n() == 1) %>%  ungroup() %>%
  dplyr::select(-c('model', 'bgd.group', 'bgd.labels')) %>%
  arrange(peakName) %>% split(f = .$label)
  
sapply(de.markers.list2, nrow)
sapply(de.markers.list2, '[[', 'peakName' )  %>% unlist() %>% head()

save_diffPeaks_fn2 = here(DATADIR, 'rdas/cell_type_models.top_diffPeaks.DESeq2.rds')
saveRDS(de.markers.list2, save_diffPeaks_fn2)

