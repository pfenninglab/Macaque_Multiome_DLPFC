ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

########################################
## load genomeAnnotation, geneAnnotation
suppressMessages(library(ArchR))
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

#################################################
## create Arrow file and comput duplicate scores
for(PROJDIR in c('ArchR_DLPFC_multiomeATAC','ArchR_DLPFC_scATAC')){
  proj = loadArchRProject(path = file.path('data/tidy_data/ArchRProjects',PROJDIR))
  varFeat = 60
  pd = getCellColData(proj)
  dimRed = names(attributes(proj)$reducedDims)
  embedNames = names(attributes(proj)$embeddings)
  
  # add iterative LSI
  iterLSIName = paste0("IterativeLSI_Peak",varFeat)
    pdf()
    proj <- addIterativeLSI( proj, useMatrix = "PeakMatrix", 
                             name = iterLSIName,
                             LSIMethod = 2, 
                             iterations = 4, # increase this if noticing subtle batch effects
                             scaleTo = 15000, # median unique fragment per cell
                             selectionMethod = 'var',
                             clusterParams = list( # See Seurat::FindClusters
                               resolution = .2, # lower this if noticing subtle batch effects
                               sampleCells = 10000,  n.start = 10), 
                             varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
                             dimsToUse = 1:30, force = TRUE)
    dev.off()
    proj = saveArchRProject(ArchRProj = proj)
  
  # add Harmony batch correction
  HarmonyName = paste0("Harmony_Peak",varFeat)
    print(HarmonyName)
    proj <- addHarmony( proj, reducedDims = iterLSIName, 
                        max.iter.harmony = 15, name = HarmonyName, 
                        groupBy = c("Sample"), force = TRUE)
  
  # add umap for harmony
  UMAPName2 = paste0("UMAP_Peak",varFeat)
    print(UMAPName2)
    proj <- addUMAP(proj, reducedDims = HarmonyName,  name = UMAPName2, 
                    nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
  
  # add clusters for Harmony
  ClustersName2 = paste0("Clusters_Peak",varFeat)
    print(ClustersName2)
    proj <- addClusters( proj, reducedDims = HarmonyName, 
                         method = "Seurat", name = ClustersName2, 
                         filterBias = TRUE,resolution = 2, force = TRUE)
  proj = saveArchRProject(proj)
}
