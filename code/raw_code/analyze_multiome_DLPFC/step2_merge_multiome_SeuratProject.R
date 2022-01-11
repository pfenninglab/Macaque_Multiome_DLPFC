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
library(harmony)
library(BSgenome.Mmulatta.UCSC.rheMac10)
DATADIR='data'

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 12 cores
nThread = min(parallel::detectCores()-2, 12)
plan("multicore", workers = nThread)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

## load the ArchR project w/ the samples and get peak counts matrix
ARCHDIR2=here(DATADIR,'tidy_data/ArchRProjects/ArchR_DLPFC_multiomeATAC')
proj_multiome = loadArchRProject(ARCHDIR2)

###############################
## grab the individual samples 
cells_fn=here(DATADIR,'tidy_data/rdas/ArchR_snATAC_DLPFC_cells.rds')
h5Seurat_fn = here(DATADIR, 'tidy_data/SeuratMultiome') %>% 
  list.files(pattern = '.h5Seurat', full.names = TRUE)
names(h5Seurat_fn) = basename(h5Seurat_fn) %>% ss('\\.', 2)
h5Seurat_fn = h5Seurat_fn[grepl('all_nuclei', h5Seurat_fn)]
obj_list = h5Seurat_fn %>% lapply(LoadH5Seurat) 

## merge the SCTransform RNA features
obj_list = obj_list %>% lapply(function(obj) {
  DefaultAssay(obj) <- "SCT"; return(obj)})
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)

## make sure seqinfo of genome and annotations the same
obj_list = obj_list %>% lapply(function(obj) {
  DefaultAssay(obj) <- "ATAC"; return(obj)})
annotations <- GetAssayData(obj_list[[1]], slot = "annotation")
granges <- GetAssayData(obj_list[[1]], slot = "ranges")
seqinfo(annotations) = seqinfo(obj_list[[1]])[seqlevels(annotations)]
seqinfo(granges) = seqinfo(obj_list[[1]])[seqlevels(granges)]
obj_list = lapply(obj_list, SetAssayData, slot = "annotation", new.data = annotations)
obj_list = lapply(obj_list, SetAssayData, slot = "ranges", new.data = granges)

## merge together and perform TF-IDF lsi scaling / SVD
combined <- merge(x = obj_list[[1]], y = obj_list[-c(1)])
combined <- combined %>% FindTopFeatures(min.cutoff = 'q1') %>% RunTFIDF() %>%  RunSVD()
combined = combined %>% RegionStats(genome = BSgenome.Mmulatta.UCSC.rheMac10)
combined = combined %>% LinkPeaks(peak.assay = "ATAC", expression.assay = "SCT",
                                  genes.use = c('SLC17A6', 'GAD1'))

## recompute RNA PCA
DefaultAssay(combined) = 'SCT'
combined <- combined %>% RunPCA(assay = 'SCT', verbose = F, features = features)

## compute Harmony on the RNA PCA and ATA LSI embeddings
combined <- RunHarmony(combined, group.by.vars = c("Sample"), project.dim = F,
                       assay.use = "SCT", reduction.save = "harmonyPCA")
combined <- RunHarmony(combined, group.by.vars = c("Sample"), project.dim = F,
                       assay.use = "ATAC", reduction.save = "harmonyLSI")

# build a joint neighbor graph using both assays
combined <- combined %>% FindMultiModalNeighbors(reduction.list = list("harmonyPCA", "harmonyLSI"),
                               dims.list = list(1:50, 2:40), verbose = TRUE, 
                               modality.weight.name = "RNA.weight")
combined <- combined %>% RunUMAP(nn.name = "weighted.nn", assay = "RNA", verbose = F)


## save the project
save_multiome_fn = here(DATADIR, 'tidy_data/SeuratMultiome',
                        paste0('mulitome_DLPFC.N',length(obj_list),'.h5Seurat'))
SaveH5Seurat(combined, save_multiome_fn, overwrite = TRUE)

# create a new UMAP using the integrated embeddings
Idents(obj_list[[1]]) = 'cell_type2'; DefaultAssay(obj_list[[1]]) = 'ATAC'
idents.plot <- levels(droplevels(factor(combined[[]]$cell_type2)))
idents.plot = idents.plot[grepl('^L[1-6]|LAMP|SST|PV|NDN|VIP', idents.plot)]

p1 <- CoveragePlot(obj_list[[1]], region = "SLC17A6", features = "SLC17A6",
                   expression.assay = "SCT", idents = idents.plot,
                   extend.upstream = 800, extend.downstream = 800)
p2 <- CoveragePlot(obj_list[[1]], region = "GAD1", features = "GAD1",
                   expression.assay = "SCT",idents = idents.plot,
                   extend.upstream = 800, extend.downstream = 800)

pdf('tmp.pdf', height = 6*4, width = 8*4)
DimPlot(obj_list[[1]], label = TRUE, repel = TRUE, group.by = 'cell_type2' , reduction = "umap") + NoLegend()
patchwork::wrap_plots(p1, p2, ncol = 1)
dev.off()





