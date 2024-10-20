suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(rtracklayer)
library(SingleCellExperiment)
library(future)
library(harmony)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

addArchRThreads(6)
plan("multiprocess", workers = 6)

PLOTDIR='figures/exploratory/jing_trackPlots/plots'

############################################
## load in ArchR project w/ just scATAC data
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
proj = loadArchRProject(path = PROJDIR2)


## sort the cell types
cellTypes = getCellColData(proj) %>%as.data.frame()%>%
  arrange(desc(grepl('^L[2-6]', Celltype2)), 
          Celltype2 %ni% c('LAMP5', 'PV.BC', 'PV.ChC', 'SST', 'NDNF', 'VIP'),
          Celltype2) %>%
  pull(Celltype2) %>% unique()

cellTypes_pal = paletteDiscrete(cellTypes)

## make the plots for Jing
pdf(here(PLOTDIR, 'DLPFC_scATAC_UMAP_all.pdf'), onefile = F)
plotEmbedding( ArchRProj = proj, embedding = "UMAP_Peak60", name = "Celltype2",
  pal = cellTypes_pal, size = 0.1, rastr = FALSE, plotAs = "points")
dev.off()

pdf(here(PLOTDIR, 'DLPFC_scATAC_UMAP_exc.pdf'), onefile = F)
plotEmbedding( ArchRProj = proj[grepl('^L[2-6]', proj$Celltype2),], 
               embedding = "UMAP_Peak60", name = "Celltype2",
               pal = cellTypes_pal, size = 0.1, rastr = FALSE, plotAs = "points")
dev.off()


pdf(here(PLOTDIR, 'DLPFC_scATAC_UMAP_inh.pdf'), onefile = F)
plotEmbedding( ArchRProj = proj[ proj$Celltype2 %in% c('LAMP5', 'PV.BC', 'PV.ChC', 'SST', 'NDNF', 'VIP'),], 
               embedding = "UMAP_Peak60", name = "Celltype2",
               pal = cellTypes_pal, size = 0.1, rastr = FALSE, plotAs = "points")
dev.off()



############################################
## load in Seurat object w/ just scRNA data
obj = LoadH5Seurat(here('data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104', 
                   'all_nuclei_final_withSCTransform.h5Seurat')) 

obj = obj[,!is.na(obj$cell_type2)]
Idents(obj) = 'cell_type2'
sum(is.na(obj$cell_type2))

obj = obj %>% RunPCA(verbose = FALSE) %>% 
  ## remove sample-wise batch effects in the dimension reduction
  RunHarmony("Sample", plot_convergence = F, assay.use = "SCT")
  ## make a UMAP clustering
obj = obj %>% RunUMAP(reduction = "Seurat..ProjectDim.SCT.harmony", dims = 1:30)

## make the plots for Jing
pdf(here(PLOTDIR, 'DLPFC_scRNA_UMAP_all.pdf'), onefile = F)
DimPlot(obj, cols = cellTypes_pal[cellTypes]) + 
  theme(legend.position = 'bottom')
dev.off()

pdf(here(PLOTDIR, 'DLPFC_scRNA_UMAP_exc.pdf'), onefile = F)
DimPlot(obj, cols = cellTypes_pal[cellTypes], 
        cells =colnames(obj)[grepl('^L[2-6]', obj$cell_type2)]) + 
  theme(legend.position = 'bottom')
dev.off()

pdf(here(PLOTDIR, 'DLPFC_scRNA_UMAP_inh.pdf'), onefile = F)
DimPlot(obj, cols = cellTypes_pal[cellTypes], 
        cells =colnames(obj)[obj$cell_type2 %in% c('LAMP5', 'PV.BC', 'PV.ChC', 'SST', 'NDNF', 'VIP')]) + 
  theme(legend.position = 'bottom')
dev.off()


# ## get the GeneScoreMatrix
# se_geneScore = getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
# rowMax(assays(se_geneScore)[['GeneScoreMatrix']]) %>% summary()
# rownames(se_geneScore) = rowData(se_geneScore)$name
# counts = round(assays(se_geneScore)[['GeneScoreMatrix']])
# counts[1:5,1:5] ## make sure has genes and cell barcodes
# 
# ## get the cell metadata from ArchR
# meta.data = getCellColData(proj) %>% as.data.frame()
# colnames(counts) %in% rownames(meta.data) %>% table()
# rownames(meta.data) %in% colnames(counts) %>% table()
# 
# ## create the Seurat object version of this 
# obj_geneScore = CreateSeuratObject( counts = counts, meta.data = meta.data, 
#                                     assay = 'ATAC',
#                                     project = "PFC_scATAC_geneScore")

