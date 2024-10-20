library(rtracklayer)
library(Biostrings)
library(here)

# 
fa_501 = import(here('data/tidy_data/celltype_specific_enhancers/fasta/rheMac10_DLPFC.noncoding_peak.fa'), format = 'fasta')
all(width(fa_501)==501) # TRUE

# subset from the 2nd base to the 501th base
fa_500 = subseq(fa_501, 2, 501)
all(width(fa_500)==500) # TRUE

export(fa_500, con = here('data/tidy_data/celltype_specific_enhancers/fasta/rheMac10_DLPFC.noncoding_peak.500.fa'), format = 'fasta')
