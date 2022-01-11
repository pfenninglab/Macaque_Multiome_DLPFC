ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
library(here)

source('/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/narrowPeakFunctions.R')


##############################
### read in ArchR project ####
DATADIR='data/tidy_data'
CODEDIR='code/raw_code/label_multiomeATAC_cells'
LABEL='Macaque_PFC'; GENOME = 'rheMac10'; 


###############################################
# 1) grab peaks from the N=5 peak calls ##
PEAKDIR=here('data/raw_data','peak', 'old_peaks')
narrowPeak_mmul10_fn = list.files(here(PEAKDIR), full.names = T, 
                                  pattern = '.GenBankRheMac8.narrowPeak.gz')
names(narrowPeak_mmul10_fn) = basename(narrowPeak_mmul10_fn) %>%
  gsub('multiomeATAC_DLPFC.|.GenBankRheMac8.narrowPeak.gz','', .)

# read peaks to narrowPeak file
peakList = lapply(narrowPeak_mmul10_fn, import)
names(peakList) = ifelse(grepl('Astrocytes', names(peakList)), 'Astro', 
                  ifelse(grepl('Endothelial', names(peakList)) ,  'Endo',
                  ifelse(grepl('Oligos', names(peakList)) ,  'Oligo',
                  ifelse(grepl('Oligo_Pre', names(peakList)) , 'OPC',names(peakList)))))
names(peakList) = gsub('PVALB', 'PV', names(peakList))                        

###############################################
# 2) grab peaks from the N=5 peak calls ##
PEAKDIR2=here('data/raw_data','peak')
narrowPeak_mmul10_fn = list.files(here(PEAKDIR2), full.names = T, 
                                  pattern = '.GenBankRheMac8.narrowPeak.gz')
names(narrowPeak_mmul10_fn) = basename(narrowPeak_mmul10_fn) %>%
  gsub('Macaque_PFC.|.GenBankRheMac8.narrowPeak.gz','', .)
# read peaks to narrowPeak file
peakList2 = lapply(narrowPeak_mmul10_fn, import)

##########################################
# 2) grab peaks from the N=5 peak calls ##
df = list(numPeak_N2 = data.frame(numCells =  lengths(peakList), celltype = names(peakList)),  
          numPeak_N5 = data.frame(numCells =  lengths(peakList2), celltype = names(peakList2))) %>%
  data.table::rbindlist(idcol = 'sample.size') %>%
  pivot_wider(!sample.size, names_from = 'sample.size', values_from = numCells, values_fill = 0) %>%
  as.data.frame() %>% arrange(celltype)
  

  
  
  