### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

# for ease of not in set operation
`%ni%` <- Negate(`%in%`)

options(stringsAsFactors = F)
library(rtracklayer)
library(tidyverse)
library(ArchR)
library(here)
source('/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/narrowPeakFunctions.R')

genome = 'rheMac10'


# compute geometric mean of a vector
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

################
# Try this out
PROJDIR = "/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/ArchRProjects/signal_corrected"
proj = PROJDIR %>% loadArchRProject()


peak_rds_fn = list.files(path = here('data', 'tidy_data','ArchRProjects', 'signal_corrected','PeakCalls'), 
                         full.names = T, pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

celltypes = unique(getCellColData(proj)$allen_labels)

writeToPeak <- function(gr, genome, fastaFile, labelFile, ycol = 'signalValueGM'){
  require(BSgenome.Mmulatta.UCSC.rheMac10)
  require(Biostrings)
  bsgenome = BSgenome.Mmulatta.UCSC.rheMac10
  names(gr) = paste0(genome, ':', seqnames(gr), ':', start(gr), '-', end(gr))
  
  ## write sequence to fasta file
  # if(grepl('.gz$', fastaFile)){
  #   writeXStringSet(seq,file=fastaFile, compress= TRUE)
  # }else {
  #   writeXStringSet(seq,file=fastaFile)
  # }
  
  # write peaks to narrowPeak file
  outRanges = write_GRangesToNarrowPeak(gr = gr, file = fastaFile, genome = genome)
  
  
  ## write numeric labels for sequences to txt file
  write.table(gr %>% as.data.frame() %>% dplyr::select(all_of(ycol)), 
              file = labelFile, quote = FALSE, 
              sep = "\t", row.names = TRUE, col.names = FALSE)
}

for (cell in celltypes){
  # get the reproducible peaks across cell types
  repSummits_fn = list.files(file.path(PROJDIR, 'PeakCalls', 'ReplicateCalls'), 
                             pattern = cell, full.names = TRUE)
  repSummits_gr = repSummits_fn %>% lapply(readRDS) %>% 
    lapply(as.data.frame) %>% rbindlist() %>% GRanges()
  
  # union of all peaks peak file
  peaks_gr = getPeakSet(proj)
  names(peaks_gr) = with(as.data.frame(peaks_gr, row.names = NULL), 
                         paste0(seqnames, ":", start, '-',end ))
  
  # filter out promoters, exonic, and filter out remotely promoter-like peaks
  table(peaks_gr$peakType)
  peaks_gr = peaks_gr[peaks_gr$peakType %in% c('Distal', 'Intronic') & 
                        peaks_gr$distToTSS > 20000 ]
  length(peaks_gr) 
  
  oo = findOverlaps(query = peaks_gr, subject = repSummits_gr)
  peaks_gr2 = cbind(peaks_gr[queryHits(oo)] %>% as.data.frame(row.names = NULL),
                    tmpsignalValue = repSummits_gr[subjectHits(oo)]$signalValue) %>%
    mutate(tmp = paste0(seqnames, ":", start, '-',end )) %>% group_by(tmp) %>%
    mutate(signalValueGM = gm_mean(tmpsignalValue)) %>%
    ungroup() %>% distinct(tmp, .keep_all = TRUE)  %>% 
    dplyr::select(-c(tmp, tmpsignalValue)) %>% 
    filter(distToTSS > 20000, peakType %in% c('Distal', 'Intronic')) %>%
    GRanges()
  names(peaks_gr2) = with(as.data.frame(peaks_gr2, row.names = NULL), 
                          paste0(seqnames, ":", start, '-',end ))
  
  # the peaks w/o any signal values
  peaks_gr3 = peaks_gr[names(peaks_gr) %ni% names(peaks_gr2)]
  mcols(peaks_gr3)$signalValueGM = 0
  # add these to the set
  peaks_gr4 = c(peaks_gr2, peaks_gr3)
  
  summary(peaks_gr2$signalValueGM)
  summary(peaks_gr4$signalValueGM)
  
  OUTDIR=here('data/tidy_data/signal_corrected_peaks')
  
  narrowPeak_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, cell, sep = '_'), '_signalValueGM.narrowPeak.gz'))
  label_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, cell, sep = '_'), '_signalValueGM.txt.gz'))
  
  out = writeToPeak(gr = peaks_gr2, fastaFile = narrowPeak_fn, 
               labelFile = label_fn, genome = genome, ycol = 'signalValueGM')
  
  
  narrowPeak_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, cell, sep = '_'), '_signalValueMoreZeroes.narrowPeak.gz'))
  label_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, cell, sep = '_'), '_signalValueMoreZeroes.txt.gz'))
  out = writeToPeak(gr = peaks_gr4, fastaFile = narrowPeak_fn, labelFile = label_fn, genome = genome, ycol = 'signalValueGM')
  
}














#########

celltype = unique(getCellColData(proj)$allen_labels)
cell = 'SST'
for (cell in celltypes){
  # replicate peak calls file
  repSummits_fn = list.files(file.path(PROJDIR, 'PeakCalls', 'ReplicateCalls'), 
                               pattern = cell, full.names = TRUE)
  repSummits_gr = repSummits_fn %>% lapply(readRDS) %>% 
    lapply(as.data.frame) %>% rbindlist() %>% GRanges()
  
  # compute geometric mean of "signalValue" column for a reproducible peak
  oo = findOverlaps(query = peaks_gr, subject = repSummits_gr)
  peaks_gr2 = cbind(peaks_gr[queryHits(oo)] %>% as.data.frame(row.names = NULL),
                    tmpsignalValue = repSummits_gr[subjectHits(oo)]$signalValue) %>%
    mutate(tmp = paste0(seqnames, ":", start, '-',end )) %>% group_by(tmp) %>%
    mutate(signalValueGM = gm_mean(tmpsignalValue)) %>%
    ungroup() %>% distinct(tmp, .keep_all = TRUE)  %>% 
    dplyr::select(-c(tmp, tmpsignalValue)) %>% 
    filter(distToTSS > 20000, peakType %in% c('Distal', 'Intronic')) %>%
    GRanges()
  names(peaks_gr2) = with(as.data.frame(peaks_gr2, row.names = NULL), 
                         paste0(seqnames, ":", start, '-',end ))
  
  # the peaks w/o any signal values
  peaks_gr3 = peaks_gr[names(peaks_gr) %ni% names(peaks_gr2)]
  mcols(peaks_gr3)$signalValueGM = 0
  # add these to the set
  peaks_gr4 = c(peaks_gr2, peaks_gr3)
  
  summary(peaks_gr2$signalValueGM)
  summary(peaks_gr4$signalValueGM)
  
  #####################
  ## split this peakset
  #fold = 'fold1'
  OUTDIR=here('data/tidy_data/signal_corrected_peaks')

  ## peaks detectable in cell types with non-zero signalValue
  peaks_grList = splitPeakSet(peaks_gr2, testSet = testSet, validSet = folds[[fold]], useCol = 'seqnames')
  lengths(peaks_grList)
  
  ## write summarised experiment to fasta
  narrowPeak_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, 
                                         names(peaks_grList), cell, sep = '_'), '_signalValueGM.narrowPeak.gz'))
  label_fn = here(OURDIR, paste0(paste('Macaque_DLPFC', genome, 
                                         names(peaks_grList), cell, sep = '_'), '_signalValueGM.txt.gz'))
  out = mapply(writeToFasta, gr = peaks_grList, narrowPeakFile = narrowPeak_fn, 
               labelFile = label_fn, genome = genome, ycol = 'signalValueGM')
  
  ## peaks with a lot more 0's
  peaks_grList = splitPeakSet(peaks_gr4, testSet = testSet, validSet = folds[[fold]], useCol = 'seqnames')
  lengths(peaks_grList)
  
  ## write summarised experiment to fasta
  narrowPeak_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, 
                                         names(peaks_grList), cell, sep = '_'), '_signalValueMoreZeroes.narrowPeak.gz'))
  label_fn = here(OUTDIR, paste0(paste('Macaque_DLPFC', genome, 
                                         names(peaks_grList), cell, sep = '_'), '_signalValueMoreZeroes.txt.gz'))
  out = mapply(writeToFasta, gr = peaks_grList, fastaFile = narrowPeak_fn, 
               labelFile = label_fn, genome = genome, ycol = 'signalValueGM')
}