# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(ArchR)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(tidyverse)
library(here)
library(future)
library(data.table)
library(BSgenome.Mmulatta.UCSC.rheMac10)
DATADIR='data'

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 12 cores
nThread = min(parallel::detectCores()-2, 12)
plan("multicore", workers = nThread)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

## grab the rheMac10 annotations of the genome
annotation_fn=file.path( "/projects/pfenninggroup/singleCell/Macaque_snATAC-seq", 
                         'macaque_snATAC-seq/Macaque_Multiome_Striatum/data',
                         'tidy_data/rdas/EnsDb.Hsapiens.v86_liftOver_rheMac10.gr.rds')
if(file.exists(annotation_fn)){
  annotation_rheMac10 = readRDS(annotation_fn)
} else {
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  seqlevels(annotation) = paste0('chr', seqlevels(annotation))
  
  chainFile =here("/home/bnphan/resources/liftOver_chainz", 'hg38ToRheMac10.over.chain')
  ch = import.chain(chainFile)
  annotation_rheMac10 = liftOver(annotation, ch = ch)
  annotation_rheMac10 = unlist(annotation_rheMac10)
  saveRDS(annotation_rheMac10, annotation_fn)
}
######################################
## grab the snATAC cells that pass QC
cells_fn=here(DATADIR,'tidy_data/rdas/ArchR_snATAC_DLPFC_cells.rds')

if(FALSE){
  # load the all nuclei seurat object
  obj = here(DATADIR, 'tidy_data/rdas/JH_PFC_LabeledNuclei_20220104', 
             'all_nuclei_final.h5Seurat') %>% LoadH5Seurat()
  obj = subset(obj, subset = Chemistry.library == 'Multiome kit')
  table(obj$Sample)
  
  # grab the cells from the ArchR project
  ARCHDIR1=here(DATADIR,'tidy_data/ArchRProjects/ArchR_DLPFC_scATAC')
  proj_multiome = loadArchRProject(ARCHDIR1)
  cells = getCellColData(proj_multiome) %>% as.data.frame() %>%
    rownames_to_column('Barcode') %>% mutate(Barcode2 = ss(Barcode, '#', 2)) %>%
    ## filter out to snATAC cells w/ an RNA profile
    filter(Barcode %in% Cells(obj)) %>% split(f = .$Sample) %>%
    lapply(function(x) x %>% dplyr::select(Barcode, Barcode2) %>% deframe())
  # cells = cells[names(fragments_fn)]
  lengths(cells) 
  
  saveRDS(cells, file = cells_fn)
  
  # grab the peak matrix from the ArchR project
  ARCHDIR2=here(DATADIR,'tidy_data/ArchRProjects','ArchR_DLPFC_multiomeATAC')
  if(dir.exists(ARCHDIR2)) {
    proj_multiome = loadArchRProject(ARCHDIR1)
  } else {
    proj_multiome = subsetArchRProject( proj_multiome, outputDirectory = ARCHDIR2,
                                        cells = unlist(sapply(cells, names)))
  }
  
  ## subset the Seurat RNA profiles
  obj = subset(obj, cells = Cells(obj)[Cells(obj) %in% getCellNames(proj_multiome)])
  obj <- SplitObject(obj, split.by = "Sample") # 
  obj = obj[names(cells)]
  h5Seurat_fn = here(DATADIR, 'tidy_data/SeuratMultiome', 
                     paste0('all_nuclei.', names(obj), '.h5Seurat'))
  mapply(SaveH5Seurat, obj = obj, file = h5Seurat_fn, overwrite = TRUE)
} 

## load the ArchR project w/ the samples and get peak counts matrix
ARCHDIR2=here(DATADIR,'tidy_data/ArchRProjects/ArchR_DLPFC_multiomeATAC')
proj_multiome = loadArchRProject(ARCHDIR2)

## get peak counts matrix w/ old peaks
feats_multiome = getMatrixFromProject(proj_multiome, useMatrix = "PeakMatrix", binarize = FALSE)

## grab the peakset from the snATAC analyses
peaks_fn=here(DATADIR,'tidy_data/rdas/ArchR_DLPFC_peakSet.gr.rds')
if(!file.exists(peaks_fn)) {
  peaks_gr = rowRanges(feats_multiome)
  saveRDS(peaks_gr, file = peaks_fn)
} else {
  peaks_gr= peaks_fn %>% readRDS()
}



##############################################
## 1) load the h5Seurat objects of the snRNA-seq

for( i in 1:4) {
  cells = readRDS(file = cells_fn)[[i]]
  cellNames = setNames(names(cells), cells)
  sample_name = ss(names(cells), '#') %>% unique()
  
  ## see if Seurat object already has ATAC profile
  h5Seurat_fn = here(DATADIR, 'tidy_data/SeuratMultiome', 
                     paste0('all_nuclei.', sample_name,'.h5Seurat'))
  
  hfile <- Connect(h5Seurat_fn)
  alreadyMade = "ATAC" %in% names(hfile$index())
  hfile$close_all()
  if(alreadyMade) next
  obj = LoadH5Seurat(h5Seurat_fn)
  
  #######################################
  ## 2) initiate the Signac ChromAssay 
  ## make fragments for each fragments file
  fragments_fn = here(DATADIR, 'raw_data/bed', paste0( sample_name,'.corrected.tsv.gz'))
  names(fragments_fn) = basename(fragments_fn) %>% ss('\\.')
  fragments_fn = fragments_fn[sample_name]
  fragments <- CreateFragmentObject(path = fragments_fn, cells = cells)
  
  ###################################################
  ## 2) merge the ATAC counts w/ the Seurat RNA object
  feats = feats_multiome[,cellNames]
  
  ## initiate ATAC assay
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = assays(feats)$PeakMatrix, fragments = fragments, ranges = rowRanges(feats), 
    genome = seqinfo(BSgenome.Mmulatta.UCSC.rheMac10), annotation = annotation_rheMac10)
  
  ## SCTransform RNA object
  DefaultAssay(obj) <- "RNA"
  obj <- obj %>% 
    SCTransform(method = "glmGamPoi", vars.to.regress = c("percent.mt", 'percent_ribo')) %>% 
    RunPCA(assay = 'SCT', verbose = F)
  SaveH5Seurat(obj, h5Seurat_fn, overwrite = T)
}