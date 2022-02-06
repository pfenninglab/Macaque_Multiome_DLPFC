### conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
library(here)
library(tidyverse)
library(data.table)
library(Biostrings)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(TFutils)
library(motifmatchr)
library(TFBSTools)
library(universalmotif)
library(SummarizedExperiment)
library(GenomicRanges)

DATADIR='data/tidy_data/celltype_specific_enhancers'

hocomoco.mono = hocomoco.mono %>% dplyr::rename('MOTIF_ID' = 'Model')
save_model_fn = here(DATADIR, 'rdas/cell_type_models_to_train_DLPFC.rds')
cell_type_df = readRDS(save_model_fn)

#############################################################################
## gather the memechip_motifs
motif_fn = list.files(here(DATADIR, 'meme'), recursive = T, pattern = 'summary.tsv', 
                      full.names = T)
names(motif_fn) = dirname(motif_fn) %>% basename()

## read in motifs
motif_df = lapply(motif_fn, fread) %>% rbindlist(idcol = 'model') %>%
  mutate(MOTIF_ID = ifelse(MOST_SIMILAR_MOTIF=='', MOTIF_ID, MOST_SIMILAR_MOTIF)) %>%
  left_join(hocomoco.mono) %>% left_join(cell_type_df)

save_meme_fn = here(DATADIR, 'rdas', 'meme-chip_celltype_differentialMotif_DLPFC.rds')
saveRDS(motif_df, save_meme_fn)


#############################################################################
## count motifs for cell types in candidate enhancers
peak_enh_fn = here(DATADIR, 'fasta', 'rheMac10_DLPFC.noncoding_peak.fa')
enh_candidate_fa = readDNAStringSet(peak_enh_fn)
enh_candidate_gr = names(enh_candidate_fa) %>% ss('rheMac10:|:250$',2) %>% GRanges()

# Get motif matches for example motifs in peaks 
hocomoco_fn = '/home/bnphan/resources/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme'
PWMs = read_meme(hocomoco_fn) %>% convert_motifs('TFBSTools-PWMatrix')
names(PWMs) = sapply(PWMs, function(x) attributes(x)$name)
# idxKeep = sapply(PWMs, function(x) attributes(x)$name %in% unlist(motifid_unique))
PWMlist = do.call('PWMatrixList', PWMs)

# scan candidate enhancers for cell type-unique motifs
motif_ix <- matchMotifs(PWMlist, enh_candidate_fa) 
motif_counts = motifMatches(motif_ix) # Extract matches matrix from result
colnames(motif_counts) = names(PWMlist)
rownames(motif_counts) = names(enh_candidate_fa)
motif_se = SummarizedExperiment(assays = list(counts = motif_counts),
                                rowRanges = enh_candidate_gr)

save_motif_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.noncoding_peak.HOCOMOCO11_motifMatches.se.rds')
saveRDS(motif_se, save_motif_fn)


##########################################################
## count motifs for cell types in candidate enhancers

## get motifs that only match in 1 cell type
genesymbol_unique = motif_df %>% 
  group_by(label) %>% filter(!duplicated(MOTIF_ID)) %>% filter(!is.na(HGNC)) %>% 
  group_by(HGNC) %>% filter(n() == 1) %>% arrange(HGNC) %>% 
  split(x = .$HGNC, f =.$label)

motifid_unique = motif_df %>% 
  group_by(label) %>% filter(!duplicated(MOTIF_ID)) %>% filter(!is.na(HGNC)) %>% 
  group_by(HGNC) %>% filter(n() == 1) %>% arrange(HGNC) %>% 
  split(x = .$MOTIF_ID, f =.$label) 

### count how many celltype-specific motifs 
names(PWMlist) %in% unlist(motifid_unique) %>% table()
unlist(motifid_unique) %in% names(PWMlist) %>% table()

motif_unique_sums = sapply(motifid_unique, function(motif_ids){
  if(length(motif_ids) > 1)
    ret = rowSums(motif_counts[,motif_ids])
  else 
    ret = motif_counts[,motif_ids]
  return(ret)
})

motif_unique_sums_se = SummarizedExperiment(assays = list(counts = motif_unique_sums),
                                rowRanges = enh_candidate_gr)

save_motif_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.noncoding_peak.HOCOMOCO11_celltypeUniqueMotifs.se.rds')
saveRDS(motif_unique_sums_se, save_motif_fn)


##########################################################
## count motifs for cell types in candidate enhancers

## get motifs found in any
genesymbol_celltype = motif_df %>% 
  group_by(label) %>% filter(!duplicated(MOTIF_ID)) %>% filter(!is.na(HGNC)) %>% 
  arrange(HGNC) %>% split(x = .$HGNC, f =.$label)

motifid_celltype = motif_df %>% 
  group_by(label) %>% filter(!duplicated(MOTIF_ID)) %>% filter(!is.na(HGNC)) %>% 
  split(x = .$MOTIF_ID, f =.$label) 

### count how many celltype-specific motifs 
names(PWMlist) %in% unlist(motifid_celltype) %>% table()
unlist(motifid_celltype) %>% unique() %in% names(PWMlist) %>% table()

motif_celltype_zscore = sapply(motifid_celltype, function(motif_ids){
  if(length(motif_ids) > 1)
    ret = rowSums(motif_counts[,motif_ids])
  else 
    ret = motif_counts[,motif_ids]
  ret = (ret - mean(ret))/sd(ret)
  return(ret)
})

motif_celltype_zscore_se = SummarizedExperiment(assays = list(zscore = motif_celltype_zscore),
                                     rowRanges = enh_candidate_gr)

save_motif_fn = here(DATADIR, 'rdas', 'rheMac10_DLPFC.noncoding_peak.HOCOMOCO11_celltypeMotifs.se.rds')
saveRDS(motif_celltype_zscore_se, save_motif_fn)


