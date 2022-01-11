suppressMessages(library(ArchR))
library(here); library(tidyverse)
library(rtracklayer)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/exploratory/jing_trackPlots/plots'

######################################
## load in ArchR project w/ multiome
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_Multiome_DLPFC")
proj = loadArchRProject(path = PROJDIR2)

## sort the cell types
cellTypes = getCellColData(proj) %>%as.data.frame()%>%
  arrange(desc(grepl('^L[2-6]', predictedGroup_RNA2RNACo)), 
          predictedGroup_RNA2RNACo %ni% c('LAMP5', 'PVALB.BC', 'PVALB.ChC', 'SST', 'NDNF', 'VIP'),
          predictedGroup_RNA2RNACo) %>%
  pull(predictedGroup_RNA2RNACo) %>% unique()


# proj$predictedGroup_RNA2RNACo = factor(proj$predictedGroup_RNA2RNACo, cellTypes)

PeakSet_L23 <- getMarkerFeatures(
  proj, useMatrix = "PeakMatrix",  groupBy = "predictedGroup_RNA2RNACo",  
  testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c("L2.CUX2.MEIS2", "L3.CUX2.RORB"))
markerSet_L23 <- getMarkers(PeakSet_L23, cutOff = "FDR <= 0.1 & Log2FC >1", returnGR = TRUE)

markerGenes  <- c( "CUX2", 'TYR', 'TBX15', 'POU3F1', 
                   'PCP4', 'OPRK1', 'ITGA8', 'NKD1', 'SYT6', 'NR4A2')

pdf(file.path(PLOTDIR, 'fig1_trackPlots_marker_genes_exc.pdf'), height = 5, width = 8)
pList <- plotBrowserTrack( ArchRProj = proj, groupBy = "predictedGroup_RNA2RNACo", 
  geneSymbol = markerGenes, upstream = 50000,downstream = 50000, minCells = 100,
  features = markerSet_L23, useGroups = cellTypes, loops = getPeak2GeneLinks(proj, corCutOff = 0.25))
for (p in pList) { grid::grid.draw(p); grid::grid.newpage()}
dev.off()


## interneuron subtypes 
PeakSet_INH <- getMarkerFeatures(
  proj, useMatrix = "PeakMatrix",  groupBy = "predictedGroup_RNA2RNACo",  
  testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c('LAMP5', 'PVALB.BC', 'SST', 'VIP'))
markerSet_INH <- getMarkers(PeakSet_INH, cutOff = "FDR <= 0.1 & Log2FC >1", returnGR = TRUE)

markerGenes2  <- c( "LAMP5", "PVALB", 'SST', 'VIP', 'NDNF')

pdf(file.path(PLOTDIR, 'fig1_trackPlots_marker_genes_inh.pdf'), height = 5, width = 8)
pList <- plotBrowserTrack( ArchRProj = proj, groupBy = "predictedGroup_RNA2RNACo", 
  geneSymbol = markerGenes2, upstream = 50000, downstream = 50000, minCells = 100,
  features = markerSet_INH, useGroups = cellTypes, loops = getPeak2GeneLinks(proj, corCutOff = 0.25))
for (p in pList) {grid::grid.draw(p); grid::grid.newpage()}
dev.off()


### chandelier PV cells
PeakSet_PV <- getMarkerFeatures(
  proj, useMatrix = "PeakMatrix",  groupBy = "predictedGroup_RNA2RNACo",  
  testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c('PVALB.BC'))
markerSet_PV <- getMarkers(PeakSet_PV, cutOff = "FDR <= 0.1 & Log2FC >.2", returnGR = TRUE)

markerGenes3  <- c( "PVALB", 'UNC5B', 'TRPS1', 'NFIB','RORA', 'SCN1A')

pdf(file.path(PLOTDIR, 'fig1_trackPlots_marker_genes_pvalb.pdf'), height = 5, width = 8)
pList <- plotBrowserTrack( ArchRProj = proj, groupBy = "predictedGroup_RNA2RNACo", 
                           useGroups  = cellTypes, features = markerSet_PV,
                           geneSymbol = markerGenes3, upstream = 50000, downstream = 50000, minCells = 100,
                           loops = getPeak2GeneLinks(proj, corCutOff = 0.25))
for (p in pList) { grid::grid.draw(p); grid::grid.newpage()}
dev.off()

