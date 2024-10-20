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

##############################
### read in ArchR project ####
PROJDIR='data/tidy_data'
LABEL='Macaque_PFC'; GENOME = 'rheMac10'; 

#######################################################
# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = here(DATADIR,'ArchRProjects', 'ArchR_DLPFC_scATAC','PeakCalls'), 
                           full.names = T, pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# label the summit and peak name using rheMac10 coordinates
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

###############################################
# create directory and narrowPeak file names ##
PEAKDIR2=here('data/raw_data','peak_rheMac10')
system(paste('mkdir -p',PEAKDIR2))
narrowPeak_mmul10_fn = here(PEAKDIR2, paste0(LABEL, '.', names(peakList), 
                                          '.rheMac10.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList, file = narrowPeak_mmul10_fn, genome = GENOME)

##################################################################
## liftover peaks to rheMac8 (the version in Cactus hal file) ##
chainFile =here("/home/bnphan/resources/liftOver_chainz", 'rheMac10ToRheMac8.over.chain')
peakList_rheMac8 = lapply(peakList, liftOver_narrowPeak, chainFile = chainFile)

# map rheMac8 UCSC chr to GenBank chromosome names (for halper) ##
CHRDIR = '/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/Zoonomia_data/tables'
chrmap_fn = here(CHRDIR, 'GCF_000772875.2_Mmul_8.0.1_assembly_report.txt.gz')
map = read_tsv(chrmap_fn, skip = 30) %>%
  rename_with(gsub, pattern ='# ',replacement =  '') %>%
  rename_with(make.names) %>%
  filter(Assigned.Molecule != 'na') %>%
  filter(GenBank.Accn != 'na') %>% 
  select(GenBank.Accn, UCSC.style.name) %>%
  column_to_rownames(var = "UCSC.style.name")

peakList_rheMac8_renamed = lapply(peakList_rheMac8, function(gr){
  ret = gr %>% as.data.frame() %>%
    filter(seqnames %in% rownames(map)) %>%
    mutate(seqnames = map[as.character(seqnames),'GenBank.Accn']) %>% 
    GRanges()
  return(ret)
}) %>% GRangesList()
seqnames(peakList_rheMac8_renamed[[2]])

###############################################
# create directory and narrowPeak file names ##
PEAKDIR=here('data/raw_data','peak')
dir.create(PEAKDIR, showWarnings = F, recursive = T)
narrowPeak_fn = here(PEAKDIR, paste0(LABEL, '.', names(peakList_rheMac8), 
                                          '.GenBankRheMac8.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList_rheMac8_renamed, 
                   file = narrowPeak_fn, genome = GENOME)

sapply(outRanges, function(gr){
  tmp = gr %>% as.data.frame() %>%
    mutate(col = paste0(seqnames, start, end, peak, name)) %>% 
    filter(duplicated(col))
  return(nrow(tmp))
})



# take care of w/ Carson's 

# #####################################
# # halLiftOver and HALPER the peaks ##
# halmapper_script = '/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh'
# system('mkdir -p logs')
# sbatch = 'sbatch -p pfen1 -w compute-1-40 --mem 10G'
# SOURCE_SPECIES = 'Macaca_mulatta'
# TARGET_SPECIES = c('Homo_sapiens')
# target_species = paste('-t', TARGET_SPECIES)
# source_species = paste('-s', SOURCE_SPECIES)
# outdir = paste('-o', here('data/raw_data', 'halper'))
# peak_files = paste('-b',narrowPeak_fn)
# 
# # paste the parameter calls together
# thecall = paste(sbatch, halmapper_script, 
#                 source_species, 
#                 rep(target_species, each= length(peak_files)), 
#                 outdir, 
#                 rep(peak_files, times = length(target_species)))
# cat(thecall, file= here(CODEDIR, 'step5b_run_halmaper_Peanut.sh'), sep = '\n')
# system(paste('chmod u+x', here(CODEDIR,'step5b_run_halmaper_Peanut.sh')))
# # system('./step5b_run_halmaper_Peanut.sh')
# 

