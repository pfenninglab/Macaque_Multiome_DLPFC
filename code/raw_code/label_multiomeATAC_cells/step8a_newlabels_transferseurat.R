suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(Seurat))
library(SummarizedExperiment)
suppressMessages(library(SingleCellExperiment))
library(tidyverse)
library(tidyseurat)
library(ggplot2)
library(sctransform)
library(SeuratDisk)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 24)

# load in annotations
annotations = readRDS('/projects/pfenninggroup/allen_human_brain/data/human_annotations.rds')

## load in the project
PROJDIR='../../../data/tidy_data/ArchRProjects'
proj = loadArchRProject(file.path(PROJDIR,'ArchR_DLPFC_multiomeATAC'), showLogo = FALSE)
proj_data = getCellColData(proj) %>% as.data.frame()
row_names = rownames(proj_data)

## load in the Seurat object for this project
obj = LoadH5Seurat('/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/rdas/JH_PFC_LabeledNuclei_20220104/all_nuclei_final_withSCTransform.h5Seurat')

## process the annotations
obj <- obj %>%
     NormalizeData(verbose = FALSE) %>%
     FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
     ScaleData(verbose = FALSE) %>%
     RunPCA(npcs = 30, verbose = FALSE)

 ## process the Seurat object
annotations <- annotations %>%
     NormalizeData(verbose = FALSE) %>%
     FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
     ScaleData(verbose = FALSE) %>%
     RunPCA(npcs = 30, verbose = FALSE)

anchors <- FindTransferAnchors(reference = annotations, query = obj,
                               dims = 1:30, reference.reduction = "pca", normalization.method = 'SCT')

## transfer the annotations
predictions <- TransferData(anchorset = anchors, refdata = annotations$subclass.id,
                            dims = 1:30)
  # find out the weaker predictions and mark them as unknown
  ind_unknown = which(predictions$prediction.score.max < 0.5)
  predictions$predicted.id[ind_unknown] = 'Unknown'
predictions <- predictions[1]
table(predictions$predicted.id)
obj <- AddMetaData(obj, metadata = predictions[1])
table(obj$predicted.id)

## save the RNA Seurat object from the Multiome object
save_h5 = '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/SeuratMultiome/allen_human_annotated.h5Seurat'
SaveH5Seurat(obj, save_h5, overwrite = T)

## save just the metadata with cell barcodes
df = obj[[]] %>% as_tibble() %>%
  mutate(cellNames = Cells(obj))
head(df)
save_metadata_fn = '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/rdas/allen_human_annotated.meta.rds'
saveRDS(df, save_metadata_fn)

##DISREGARD EVERYTHING BELOW THIS FOR NOW##


# #observing that in the column we want to use as the reference for label transfer, only a small subset
# #of the loaded annotations are in the project, which is fine, just something to think about
# table(proj_data$predictedCell_RNA2ATACCo %in% colnames(annotations))
# table(colnames(annotations) %in% proj_data$predictedCell_RNA2ATACCo)
# 
# #making the seurat object a tibble to make manipulation easier
# annotations = annotations %>% as_tibble()
# 
# #transfer labels based on cell barcodes
# proj_data$sample_name = proj_data$predictedCell_RNA2ATACCo
# proj_data = left_join(proj_data, annotations, by='sample_name')
# 
# #drop the extra columns that we don't need, and rename columns 
# proj_data <- proj_data[c(-86:-147)]
# proj_data <- proj_data[c(-40:-84)]
# names(proj_data)[names(proj_data) == 'subclass.id.y'] <- 'allen_labels'
# table(proj_data$allen_labels, proj_data$subclass.id.x)
# rownames(proj_data) <- row_names
# head(proj_data, 2)
# 
# proj <- addCellColData(ArchRProj = proj, data = proj_data$allen_labels, cells = rownames(proj_data), name = "allen_labels", force = TRUE)
# 
# saveArchRProject(proj, outputDirectory = file.path(PROJDIR,'relabeled_DLPFC'))
# 
