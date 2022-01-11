suppressMessages(library(ArchR))
library(here)
library(tidyverse)
library(MASS)
library(broom)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

#############################################################################
### make that arrow file from bgzipped fragments file (tsv.gz or bed.gz) ####
ArrowFiles = here('data/raw_data/arrow') %>% list.files(pattern = 'arrow$', full.names = T)
PROJDIR=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")

if(!dir.exists(PROJDIR)){
  proj = ArchRProject( ArrowFiles = ArrowFiles,
                       outputDirectory = PROJDIR, 
                       copyArrows = TRUE,
                       geneAnnotation = geneAnnotation, #must be the custom rheMac10 version
                       genomeAnnotation = genomeAnnotation #must be the custom rheMac10 version
  )

  ## read in sample sheet
  pd = read_csv(here('data/raw_data/tables/SampleSheet_ATAC_DLPFC.csv')) %>%
    rename_with(make.names) %>% 
    rename_all(funs(stringr::str_replace_all(., '^X\\.\\.', ''))) %>%
    mutate(Animal = paste0('Monkey_',stringr::str_sub(Nickname, start = 1, end = 1)),
           Sample = paste0(Nickname, '_', Region, '-',Batch)) %>%
    relocate(Sample, Region, Animal, Batch, .before = everything())
  
  ## Add metadata columns to ArchR Project
  pd = pd[match(proj$Sample, pd$Sample),]
  for(col in names(pd)){
    if( col != 'Sample'){
      proj = addCellColData( ArchRProj = proj, name = col, cells = getCellNames(proj),
        data = pd %>% pull(all_of(col)), force = TRUE)
    }}
  
  ## increase unique fragments cutoff to 10^3.4, remove cluster of low QC cell, based on QC plots
  idxSample <- BiocGenerics::which(proj$nFrags > 10^3.4)
  cellsSample <- proj$cellNames[idxSample]
  proj = subsetCells(ArchRProj = proj, cellNames = cellsSample)
  
  ## filter doublets
  proj = filterDoublets( proj, cutEnrich = 0.5, cutScore = -log10(.05), filterRatio = 1)
  
  ## fit probability regression for outlier cell detection
  df = getCellColData(proj) %>% as.data.frame() %>% mutate(tmp = seq(n())) %>%
    rownames_to_column("CellBarcode")
  
  resid_cutoff = 2
  cellsSampleKeep = split(df, df$Sample) %>% 
    map(., ~glm.nb(nFrags ~ TSSEnrichment + PromoterRatio + DoubletEnrichment, data = .)) %>% 
    map2(.x = ., .y = split(df, f = df$Sample), .f = ~augment_columns(x = .x, data = .y)) %>% 
    bind_rows() %>% mutate(outlierScore =  abs(.std.resid), 
                           isOutlier = outlierScore > resid_cutoff) %>%
    dplyr::select(-starts_with('\\.')) %>% arrange(tmp) %>%
    filter(!isOutlier) %>% pull(CellBarcode)
  
  proj = subsetCells(ArchRProj = proj, cellNames = cellsSampleKeep)
  
  ## save the filtered project
  proj = saveArchRProject(proj)
} else {
  proj = loadArchRProject(path = PROJDIR)
}

# add iterative LSI
varFeat = 200
pd = getCellColData(proj)
dimRed = names(attributes(proj)$reducedDims)
embedNames = names(attributes(proj)$embeddings)

iterLSIName = paste0("IterativeLSI",varFeat,'_ATAC')
print(iterLSIName)
if (iterLSIName %ni% dimRed){
  pdf()
  proj <- addIterativeLSI( proj, useMatrix = "TileMatrix", 
                           name = iterLSIName,
                           LSIMethod = 2, 
                           iterations = 4, # increase this if noticing subtle batch effects
                           scaleTo = 20000, # median unique fragment per cell
                           selectionMethod = 'var',
                           clusterParams = list( # See Seurat::FindClusters
                             resolution = .2, # lower this if noticing subtle batch effects
                             sampleCells = 10000,  n.start = 10), 
                           varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
                           dimsToUse = 1:30, force = FALSE)
  dev.off()
  proj = saveArchRProject(ArchRProj = proj)}

# add umap
UMAPName = paste0("UMAPI",varFeat,'_ATAC')
if (UMAPName %ni% embedNames){
  print(UMAPName)
  proj <- addUMAP(proj, reducedDims = iterLSIName, 
                  name = UMAPName, nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine", force = F)}

# add clusters
ClustersName = paste0("ClustersI",varFeat,'_ATAC')
if (ClustersName %ni% names(pd)){
  print(ClustersName)
  proj <- addClusters(proj, reducedDims = iterLSIName, method = "Seurat", 
                      algorithm = 2,
                      filterBias = TRUE, name = ClustersName, resolution = 1, force = T)
  }

# add Harmony batch correction
HarmonyName = paste0("HarmonyI",varFeat,'_ATAC')
if (HarmonyName %ni% dimRed ){
  print(HarmonyName)
  proj <- addHarmony(proj, reducedDims = iterLSIName, 
                     max.iter.harmony = 15, name = HarmonyName, 
                     groupBy = c("Animal", 'Region', 'Batch'), force = F)}

# add umap
UMAPName2 = paste0("UMAPH",varFeat,'_ATAC')
if (UMAPName2 %ni% embedNames){
  print(UMAPName2)
  proj <- addUMAP(proj, reducedDims = HarmonyName, 
                  name = UMAPName2, nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine", force = F)}

# add clusters
ClustersName2 = paste0("ClustersH",varFeat,'_ATAC')
if (ClustersName2 %ni% names(pd)){
  print(ClustersName2)
  proj <- addClusters(proj, reducedDims = HarmonyName, method = "Seurat", 
                      algorithm = 2, 
                      name = ClustersName2, resolution = 1, force = T)
}

proj = addImputeWeights(proj, reducedDims = HarmonyName)
proj = saveArchRProject(ArchRProj = proj)