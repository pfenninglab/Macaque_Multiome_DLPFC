### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
suppressMessages(library(ArchR))
library(here)
library(BSgenome.Mmulatta.UCSC.rheMac10)

DATADIR='data/tidy_data/celltype_specific_enhancers'

source(here('code/final_code/hal_scripts/narrowPeakFunctions.R'))
source(here('code/final_code/hal_scripts/gen_enh_ortholog_sets.R'))

## chromosomal splits
GENOME = 'rheMac10'
testSet = c('chr1','chr2')
folds = list(fold1 = c('chr8', 'chr9'))

#############################################################################
## export the fasta sequence of all noncoding/nonpromoter enhancer candidates
peak_enh_fn = here(DATADIR, 'fasta', paste(GENOME,'DLPFC.noncoding_peak.fa', sep = '_'))
if(!file.exists(peak_enh_fn)){
  PROJDIR=here("data/tidy_data/ArchRProjects/ArchR_DLPFC_scATAC")
  proj = loadArchRProject(path = PROJDIR)
  peak_gr = getPeakSet(proj)
  indKeep = which(peak_gr$peakType %in% c('Distal', 'Intron') & abs(peak_gr$distToTSS) > 2000)
  peak_gr = peak_gr[indKeep]
  peak_gr = peak_gr %>% addSummitCenter() %>% nameNarrowPeakRanges(genome = GENOME)
  writeGRangesToFasta (gr = peak_gr, file = peak_enh_fn, genome = GENOME)
  save_fn = here(DATADIR, 'rdas', paste(GENOME,'DLPFC.noncoding_peak.gr.rds', sep = '_'))
  saveRDS(peak_gr, file = save_fn )
}

############################################
### get the  cell type differential peaks
save_diffPeaks_fn = here(DATADIR, 'rdas/cell_type_models_diffPeakList_DLPFC.rds')
rhesus_peakList = readRDS(save_diffPeaks_fn)

# label the summit and peak name using rheMac10 coordinates
rhesus_peakList = lapply(rhesus_peakList, addSummitCenter)
rhesus_peakList = lapply(rhesus_peakList, nameNarrowPeakRanges, genome = GENOME)
rhesus_peakList = lapply(rhesus_peakList, sort)
rhesus_peakList = GRangesList(rhesus_peakList)

# trim the peaks using rheMac10 chr lengths
seqinfo(rhesus_peakList) = seqinfo(BSgenome.Mmulatta.UCSC.rheMac10)
rhesus_peakList = lapply(rhesus_peakList, trim) %>% GRangesList()

# split peak by differential log2 fold change direction
rhesus_posList = lapply(rhesus_peakList, function(gr) gr[gr$Log2FC > 0])
rhesus_negList = lapply(rhesus_peakList, function(gr) gr[gr$Log2FC < 0])

# split peaks up into test, train, and validation by folds (here just 1 fold)
rhesus_posList_split = lapply(folds, function(fold){
  ret = lapply(rhesus_posList, splitPeakSet, testSet = testSet, validSet = fold)
  return(ret)
})

rhesus_negList_split = lapply(folds, function(fold){
  ret = lapply(rhesus_negList, splitPeakSet, testSet = testSet, validSet = fold)
  return(ret)
})

# write the peaks to fasta for each differential peak set
split = names(rhesus_posList_split[[1]][[1]])
system(paste('mkdir -p',  here(DATADIR, 'fasta')))
for(cell in names(rhesus_peakList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = here(DATADIR, 'fasta', 
                             paste(GENOME, cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = rhesus_posList_split[[fold]][[cell]],  
                      file = pos.fasta_fn, genome = GENOME)
    
    # write the negatives
    neg.fasta_fn = here(DATADIR, 'fasta', 
                        paste(GENOME, cell, fold, split, 'negative.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = rhesus_negList_split[[fold]][[cell]],  
                      file = neg.fasta_fn, genome = GENOME)
}}

# save all the relevant objects for later
system(paste('mkdir -p',  here(DATADIR, 'rdas')))
save_fn = here(DATADIR, 'rdas', paste('Macaque_DLPFC', GENOME,'.rds', sep = '.'))
saveRDS(rhesus_peakList, file = save_fn )

save_fn2 = here(DATADIR, 'rdas', paste('Macaque_DLPFC', GENOME,'_negativeSplits.rds', sep = '.'))
saveRDS(rhesus_negList_split, file = save_fn2 )

save_fn3 = here(DATADIR, 'rdas', paste('Macaque_DLPFC', GENOME,'_positiveSplits.rds', sep = '.'))
saveRDS(rhesus_posList_split, file = save_fn3 )