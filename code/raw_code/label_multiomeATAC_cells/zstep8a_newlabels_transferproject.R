### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

library(Seurat)
library(SeuratDisk) ## remotes::install_github("mojaveazure/seurat-disk")
library(here)
library(tidyverse)
library(future) ## for parallel processing w/ Seurat
library(ArchR)

## load in the project
PROJDIR='../../../data/tidy_data/ArchRProjects'
proj = loadArchRProject(file.path(PROJDIR,'ArchR_DLPFC_multiomeATAC'), showLogo = FALSE)
proj_data = getCellColData(proj) %>% as.data.frame()
row_names = rownames(proj_data)

## load in the cell type annotations from allen human
annotations = readRDS('/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/rdas/allen_human_annotated.meta.rds')
annotations = annotations %>% as_tibble()

#transfer labels based on cell barcodes
proj_data$cellNames = rownames(proj_data)
proj_data = left_join(proj_data, annotations, by='cellNames')
#drop the extra columns that we don't need, and rename columns 
proj_data <- proj_data[c(-45:-70)]
names(proj_data)[names(proj_data) == 'predicted.id'] <- 'allen_labels'
table(proj_data$allen_labels, proj_data$subclass.id)
rownames(proj_data) <- row_names
head(proj_data, 2)
 
proj <- addCellColData(ArchRProj = proj, data = proj_data$allen_labels, cells = rownames(proj_data), name = "allen_labels", force = TRUE)
 
saveArchRProject(proj, outputDirectory = file.path(PROJDIR,'ArchR_DLPFC_multiomeATAC'))



