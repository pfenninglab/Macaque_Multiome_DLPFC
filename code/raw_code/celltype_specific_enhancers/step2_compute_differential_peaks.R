suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(rtracklayer)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
PLOTDIR='figures/exploratory/celltype_specific_enhancers/plots/'
DATADIR='data/tidy_data/celltype_specific_enhancers'
GENOME = 'rheMac10'

#######################
## get ArchR project
PROJDIR=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
proj = loadArchRProject(path = PROJDIR)
peak_gr = getPeakSet(proj)
names(peak_gr) = paste0(GENOME, ':', seqnames(peak_gr), ':', start(peak_gr), '-', end(peak_gr), ':250')
indKeep = which(peak_gr$peakType %in% c('Distal', 'Intron') & abs(peak_gr$distToTSS) > 5000)

save_model_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_model_fn)

model_name = cell_type_df$model
indlist = setNames(seq_along(model_name), model_name)
fgd_labels = setNames(cell_type_df$label, model_name)
bgd_labels = setNames(cell_type_df$bgd.labels %>% sapply(strsplit,';'), model_name)

save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models_diffPeakList_DLPFC.rds')
if(file.exists(save_diffPeaks_fn)){
  diffPeak_list = readRDS(save_diffPeaks_fn) 
} else {
  diffPeak_list = lapply(indlist, function(i){
    print(paste('Getting diff peaks for:', model_name[i]))
    markers = getMarkerFeatures(
      ArchRProj = proj,
      groupBy = "Celltype2",
      useGroups = fgd_labels[i],
      bgdGroups = bgd_labels[[i]],
      useMatrix = "PeakMatrix",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      maxCells = 1000, scaleTo = 10^4, k = 100,
      bufferRatio = 0.8,
      binarize = TRUE)
    rownames(markers) = names(peak_gr)
    markers <- getMarkers(markers, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.26")[[1]]
    return(GRanges(markers))
  })
  saveRDS(diffPeak_list, save_diffPeaks_fn)
}

lengths(diffPeak_list) %>% sort()
cell_type_df = cell_type_df %>% 
  mutate(numDiffPeaks.All = lengths(diffPeak_list), 
         numDiffPeaks.Up = sapply(diffPeak_list, function(gr) sum(mcols(gr)$Log2FC > 0)),
         numDiffPeaks.Dn = sapply(diffPeak_list, function(gr) sum(mcols(gr)$Log2FC < 0)), 
         posToNegRatio = numDiffPeaks.Up/numDiffPeaks.Dn,
         canTrainModel = ifelse(numDiffPeaks.All > 800, 'yes', 'no'),
         wparam = signif(posToNegRatio, digits = 2)) %>% 
  arrange(desc(canTrainModel == 'yes'), model)


save_model_fn = here(DATADIR, 'tables/cell_type_models_to_train_DLPFC.csv')
cell_type_df %>% write_csv(save_model_fn)

save_model_fn2 = here(DATADIR, 'tables/cell_type_models_to_train_DLPFC.xlsx')
cell_type_df %>% writexl::write_xlsx(save_model_fn2)

