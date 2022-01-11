ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
suppressMessages(library(ArchR))
library(tidyverse)
library(RcppML)
library(RColorBrewer)
DATADIR= '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers'

#####################################
### read in command line options ####
library(optparse)
option_list = list(
  make_option(c("-i", "--iteration"), type="numeric", default=NULL, 
              help="integer of which iteration to use", metavar="numeric"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## factor a lot of models and find optimal ranks
tol = 1e-7
param_search = seq(0, .99, .1)
l1_h = param_search[(opt$iteration - 1) / length(param_search) + 1]
l1_w = param_search[(opt$iteration - 1) %% length(param_search) + 1]

save_model_fn = file.path(DATADIR, 'rdas', 'NMF_models', 
                     paste0('DLFPC_PeakMatrix_NMF_k10-100_tol',as.character(tol),
                            '_l1h',l1_h, '_l1w', l1_w, '.rds'))
if(file.exists(save_model_fn)) break

#################################################
## grab the peak matrix of the data
save_peakmat_fn = file.path(DATADIR, 'rdas', 'DLFPC_PeakMatrix_allPeaks.rds')
print('Getting peak matrix.')
if(file.exists(save_peakmat_fn)){
  peakMat = readRDS(save_peakmat_fn)
} else{
  proj = loadArchRProject(path = file.path('data/tidy_data/ArchRProjects','ArchR_DLPFC_scATAC'))
  peakMat = getMatrixFromProject(proj, 'PeakMatrix')
  rowRanges(peakMat) = getPeakSet(proj)
  saveRDS(peakMat, save_peakmat_fn)
}

## filter to be intronic/intergenic peaks
indKeep = which(mcols(peakMat)$peakType %in% c('Distal', 'Intron'))
peakMat = peakMat[indKeep,]

########################################################
# calculate the TF-log(IDF), Adapted from Stuart et al.
print('Computing the TF-log(IDF) matrix.')
mat = assay(peakMat)
colSm <- Matrix::colSums(mat)
rowSm <- Matrix::rowSums(mat)
mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
idf   <- as(ncol(mat) / rowSm, "sparseVector") #IDF
mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat #TF-IDF
#Log transform TF-IDF
mat@x <- log(mat@x * 10000 + 1)  

##################################################
## do non-negative NMF on non-coding/non-promoter peaks
print('Running NMF across k = 10 to 100.')
dir.create(file.path(DATADIR, 'rdas','NMF_models'), recursive = T, showWarnings = F)
threads = min(parallel::detectCores()/2 %>% round(), 12)
setRcppMLthreads(threads)

print(paste0('Using L1 values for w=', l1_w, ' and h=', l1_h, '.'))
model_list = lapply(seq(10, 100, 5), function(k){
  print(paste0('Solving with rank k:', k))
  model = nmf(mat, k = k, tol = tol, maxit = 200,
              L1 = c(l1_w, l1_h), seed = 1:10,
              mask_zeros = FALSE, verbose = TRUE,
              diag = TRUE, nonneg = TRUE)
  return(model)
})
names(model_list) = paste0('K', seq(2, 100))
saveRDS(model_list, save_model_fn)
