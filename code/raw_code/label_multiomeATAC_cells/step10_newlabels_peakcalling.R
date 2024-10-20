# devtools::install_github("GreenleafLab/ArchR", 
# 	ref="release_1.0.2", repos = BiocManager::repositories())
library(ArchR)
library(parallel)
library(tidyverse)
library(here)
set.seed(1)

## add general functions at the top of the R scripts
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

options(repr.plot.width=15, repr.plot.height=8.5)

library(here)
library(Seurat)
library(SeuratDisk) ## remotes::install_github("mojaveazure/seurat-disk")
library(future) ## for parallel processing w/ Seurat

plan("multicore", workers = 2)
options(future.globals.maxSize = 40 * 1024 ^ 3) ## need to have requested more than 50Gb of RAM
library('BSgenome.Mmulatta.UCSC.rheMac10')

### load rheMac10 ArchR genome ###
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

#load in the project to use
proj = loadArchRProject('/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ArchRProjects/signal_corrected')

#column in proj to group by
to_group_by = 'allen_labels'

# proj <- addGroupCoverages(ArchRProj = proj, groupBy = to_group_by, force=T)
# 
# proj <- saveArchRProject(proj)
# 
# pathToMacs2 <- findMacs2()
# 
# proj <- addReproduciblePeakSet(
#   ArchRProj = proj, 
#   groupBy = to_group_by, 
#   reproducibility = "(n+1)/2",
#   pathToMacs2 = pathToMacs2,
#   genomeSize = 2.7e9,
#   force = TRUE
# )
# 
# proj <- saveArchRProject(proj)
# 
# proj <- addPeakMatrix(proj)
# 
# proj <- saveArchRProject(proj)
# 
# proj <- addCoAccessibility(
#   ArchRProj = proj,
#   reducedDims = "Harmony_Peak60"
# )
# 
# proj <- saveArchRProject(proj)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = to_group_by,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

#vector of genes to plot

genes = c("Bdnf", "Fos", "Egr3", "FosB", "NPAS4", "SRF", "NARP", "DHCR7", "SHANK3", "CACNB2", "Foxp1", "Slc4a10")
  
  p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = to_group_by,
    geneSymbol = genes,
    features =  markerList,
    upstream = 70000,
    downstream = 70000
  )
  
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  
p1 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = to_group_by, 
  geneSymbol = genes, 
  upstream = 70000,
  downstream = 70000,
  loops = getCoAccessibility(proj)
)

plotPDF(plotList = p1,
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


proj <- saveArchRProject(proj)



