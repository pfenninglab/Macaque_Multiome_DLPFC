suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

## read in gene annotations
annot = import('/home/bnphan/resources/genomes/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3')
names(annot) = gsub('-[1-9]$', '', annot$Name)
names(annot) = gsub('_', '-', names(annot))

## load in raw counts of unlabeled multiomeRNA
sce_fn = file.path('data/tidy_data/rdas', 'raw_multiomeRNA_DLPFC_sce_list_N4.rds')
sce_list = readRDS(sce_fn)
sce = do.call('cbind', sce_list)

mcols(sce)$Symbol = gsub('gene-|-[1-9]$', '', mcols(sce)$Symbol)
sce = sce[!duplicated(mcols(sce)$Symbol)]
sce = sce[mcols(sce)$Symbol %in% names(annot), ]
rowRanges(sce) = annot[match(mcols(sce)$Symbol, names(annot))]

table(sce_list[[1]]$cell_class)
table(sce$cell_type2)

########################################################
## combine multiomeRNA with the ArchR multiomeATAC cells
PROJDIR1=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_multiomeATAC")
if(!dir.exists(PROJDIR2)){
  proj1 = loadArchRProject(path = PROJDIR1)
  proj = saveArchRProject(proj1,outputDirectory = PROJDIR2 )
} else{
  proj = loadArchRProject(path = PROJDIR2)
}


# 12299 cells have both ATAC and RNA profiles
table( proj$cellNames %in% colnames(sce))
table( colnames(sce) %in% proj$cellNames )

## Add the multiomeRNA 
## London_DLPFC-1=7604,Memphis_DLPFC-1=8705
proj <- addGeneExpressionMatrix(input = proj, seRNA = sce, force = TRUE)

sce = sce[,colnames(sce) %in% proj$cellNames ]
proj = addCellColData( ArchRProj = proj, data = sce$cell_type2 %>% as.character(),
  name = 'Celltype2_RNA', cells = colnames(sce), force = TRUE)

proj <- proj[!is.na(proj$Gex_nUMI) & !is.na(proj$Celltype2_RNA)]
proj = saveArchRProject(ArchRProj = proj)

with(getCellColData(proj), table(Celltype2, Celltype2_RNA))
table(proj$Celltype2_RNA)

###############################################
## add iterative LSI for the gene expression
pd = getCellColData(proj)
dimRed = names(attributes(proj)$reducedDims)
embedNames = names(attributes(proj)$embeddings)

iterLSIName = paste0("IterativeLSI_RNA")
print(iterLSIName)
if (iterLSIName %ni% dimRed){
  pdf()
  proj <- addIterativeLSI(
    ArchRProj = proj, 
    clusterParams = list( resolution = 0.2, 
      sampleCells = 10000, n.start = 10),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI", force = T, 
    varFeatures = 2500, firstSelection = "variable",
    binarize = FALSE, name = iterLSIName)
  dev.off()
  proj = saveArchRProject(ArchRProj = proj)
}

# add Harmony batch correction
HarmonyName = paste0("HarmonyI_RNA")
if (HarmonyName %ni% dimRed ){
  print(HarmonyName)
  proj <- addHarmony(proj, reducedDims = iterLSIName, 
                     max.iter.harmony = 15, name = HarmonyName, 
                     groupBy = c("Sample", 'Animal'), force = T)
  }

# add umap
UMAPName2 = paste0("UMAPH_RNA")
if (UMAPName2 %ni% embedNames){
  print(UMAPName2)
  proj <- addUMAP(proj, reducedDims = HarmonyName, 
                  name = UMAPName2, nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine", force = T)
  }

# add clusters
ClustersName2 = paste0("ClustersH_RNA")
if (ClustersName2 %ni% names(pd)){
  print(ClustersName2)
  proj <- addClusters(proj, reducedDims = HarmonyName, method = "Seurat", 
                      name = ClustersName2, resolution = 1, force = T)
}
proj = saveArchRProject(ArchRProj = proj,)

####################################################
## add combined ATAC + RNA dimensionality reductions
proj <- addCombinedDims(proj, reducedDims = c("IterativeLSI200_ATAC", iterLSIName), 
                        name =  "LSI_Combined")

# add Harmony batch correction
proj <- addHarmony(proj, reducedDims = 'LSI_Combined', 
                   max.iter.harmony = 20, name = 'HarmonyI_Combined', 
                   groupBy = c("Sample", 'Animal'), force = T)

# add umap
proj <- addUMAP(proj, reducedDims = 'HarmonyI_Combined', 
                name = 'UMAPH_Combined', nNeighbors = 30, minDist = 0.5, 
                metric = "cosine", force = T)

# add clusters
proj <- addClusters(proj, reducedDims = 'HarmonyI_Combined', method = "Seurat", 
                      name = "ClustersH_Combined", resolution = 4, force = T)

proj = addImputeWeights(proj, reducedDims = "LSI_Combined")
proj = saveArchRProject(ArchRProj = proj)

# add peak2gene links matrix
proj <- addPeak2GeneLinks( ArchRProj = proj,dimsToUse = 1:30,
  reducedDims = "HarmonyI_Combined", useMatrix = 'GeneExpressionMatrix',
  scaleDims = TRUE, corCutOff = 0.75, k = 100, knnIteration = 500, 
  overlapCutoff = 0.8,  maxDist = 1e+05, scaleTo = 10^4, log2Norm = TRUE)
proj = saveArchRProject(ArchRProj = proj)


