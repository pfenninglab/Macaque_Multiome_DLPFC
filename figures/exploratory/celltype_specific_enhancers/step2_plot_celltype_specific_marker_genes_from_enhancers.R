ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(tidyverse))
library(future)
library(scater)
library(here)

addArchRThreads(threads = 4)

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 16)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')


######################
## load ArchR projects
DATADIR='data/tidy_data'
h5File = file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                   'all_nuclei_final_withSCTransform.h5Seurat')
if(!file.exists(h5File)){
  obj =LoadH5Seurat(file.path('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104',
                                   'all_nuclei_final.h5Seurat'), 
                    assay = 'RNA')
  
  obj <- SCTransform(obj, method = "glmGamPoi", verbose = TRUE, 
                     vars.to.regress = c('percent_ribo', 'percent.mt'))
  obj %>% SaveH5Seurat(h5File)
} else {
  obj =LoadH5Seurat(h5File)
  DefaultAssay(obj) = 'SCT'
}

Idents(obj) = 'cell_type2'
DefaultAssay(obj) = 'SCT'

celltypes = obj[[]] %>%
  filter(!duplicated(cell_type2)) %>% mutate(cell_type2 = as.character(cell_type2)) %>%
  filter(cell_type2 != 'TH') %>% arrange(cell_type2) %>%
  mutate(cell_class = case_when(grepl('ITGA|NR4A|TBX|ALPL', cell_type2) ~ 'EXC.Lower_ET_IT',
                                grepl('PCP4|SYT6|NKD1|POU3F1', cell_type2) ~  'EXC.Lower_CT_NP',
                                grepl('^L[2-4]', cell_type2) ~ 'EXC.Upper_IT', 
                                grepl('LAMP|NDN|VIP', cell_type2) ~ 'INH.CGE',
                                grepl('SST|PV|TH', cell_type2) ~ 'INH.MGE',
                                TRUE ~ "GLIA"), 
         cell_class = factor(cell_class, c('EXC.Upper_IT', 'EXC.Lower_CT_NP', 
                                           'EXC.Lower_ET_IT', 'INH.CGE', 'INH.MGE', 'GLIA'))) %>%
  dplyr::select(cell_type2, cell_class) %>% as.data.frame() %>% 
  arrange(cell_class, cell_type2) %>% pull(cell_type2)

celltypes_col = paletteDiscrete(celltypes)


obj_agg = AggregateExpression(obj, group.by = c('Sample', 'cell_type2'), return.seurat = TRUE)
obj_agg@meta.data$celltype = ss(colnames(obj_agg), '-[0-9]_', 2) %>%
  factor(levels = celltypes)
obj_agg = SCTransform(obj_agg)

###############################################################
## get markerGenes around the candidate enhancers ranked by CNNs
FIGDIR='figures/exploratory/celltype_specific_enhancers'
DATADIR='data/tidy_data/celltype_specific_enhancers'
save_top_enh_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.candidate_celltype_enhancers.rds')
markerGeneList = readRDS(save_top_enh_fn) %>% lapply(function(df){
  df %>% dplyr::select(label, markerGene) %>% deframe()
}) %>% lapply(function(x) {
  x = x[!is.na(x)] %>% strsplit(',') %>% unlist()
  x = x[x %in% rownames(obj)]
})
markerGeneList = markerGeneList[lengths(markerGeneList) > 0]


for(celltype in names(markerGeneList)){
  plot_fn = paste0('DLPFC_markerGene_aroung_SNAIL_candidates.', celltype, '.20220214.pdf')
  out_pdf = here(FIGDIR, 'plots', 'geneplots',plot_fn)
  
  pdf(out_pdf, height = 6, width = 16)
  for(i in seq_along(markerGeneList[[celltype]])){
    gg = VlnPlot(obj_agg, features = markerGeneList[[celltype]][i],
                 group.by = 'celltype',cols = celltypes_col, pt.size = FALSE)
    print(gg+ theme(legend.position = 'none')  + 
            ggtitle(paste0('celltype: ', celltype, 
                           ', gene: ', markerGeneList[[celltype]][i], 
                           ', peak:', names(markerGeneList[[celltype]])[i]
            )))
  }
  dev.off()
}


###############################################################
## get peak2genes around the candidate enhancers ranked by CNNs
peak2geneList = readRDS(save_top_enh_fn) %>% lapply(function(df){
  df %>% dplyr::select(label, peak2gene) %>% deframe()
}) %>% lapply(function(x) {
    x = x[!is.na(x)] %>% strsplit(',') %>% unlist()
    x = x[x %in% rownames(obj)]
    })
peak2geneList = peak2geneList[lengths(peak2geneList) > 0]

dir.create(here(FIGDIR, 'plots', 'geneplots'), showWarnings = F)

for(celltype in names(peak2geneList)){
  plot_fn = paste0('DLPFC_peak2Gene_aroung_SNAIL_candidates.', celltype, '.20220214.pdf')
  out_pdf = here(FIGDIR, 'plots', 'geneplots',plot_fn)

  pdf(out_pdf, height = 6, width = 16)
  for(i in seq_along(peak2geneList[[celltype]])){
    gg = plotExpression(sce_agg, x = 'celltype', 
                        features = peak2geneList[[celltype]][i], 
                        log2_values = TRUE, colour_by = 'celltype' ) +
      scale_fill_manual(celltypes_col)
    print(gg+ theme(legend.position = 'none')  + 
            ggtitle(paste0('celltype: ', celltype, 
                           ', gene: ', peak2geneList[[celltype]][i], 
                           ', peak:', names(peak2geneList[[celltype]])[i]
            )))
  }
  dev.off()
}




