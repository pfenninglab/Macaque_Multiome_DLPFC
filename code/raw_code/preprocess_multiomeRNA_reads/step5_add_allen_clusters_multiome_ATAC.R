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

PROJDIR2=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_multiomeATAC")
proj = loadArchRProject(path = PROJDIR2)

## according to the mapping from the scATAC co-clustering w/ allen snRNA-seq
monkeyToAllen= c(Astro ='Astro', Endo = 'Endo', L2.CUX2.MEIS2 = 'L2.3', L3.CUX2.RORB = 'L4.5', L4.5.TBX15 = 'L4.5', L4.ALPL = 'L4.5', L4.TYR = 'L3', L5.6.NR4A2 = 'L6.Car3', L5.PCP4 = 'L5.6.NP', L5.POU3F1 = 'L5.ET', L6.ITGA8 = 'L6', L6.NKD1 = 'L6b', L6.SYT6 = 'L6.CT', LAMP5 = 'Lamp5', Microglia = 'Microglia', Mural = 'VLMC', NDNF = 'Lamp5', Oligo = 'Oligo', OPC = 'Opc', PV.BC = 'Pvalb', PV.ChC= 'Pvalb', SST = 'Sst',  VIP = 'Vip')

## remap to subclass.id
proj$subclass.id <- mapLabels(proj$Celltype2, newLabels = monkeyToAllen, oldLabels = names(monkeyToAllen))

## save the transferred labels
proj = saveArchRProject(proj)

table(proj$subclass.id, proj$Celltype2)