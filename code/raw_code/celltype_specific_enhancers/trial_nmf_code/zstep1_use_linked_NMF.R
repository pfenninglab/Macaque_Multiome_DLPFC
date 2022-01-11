ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
suppressMessages(library(ArchR))
library(tidyverse)
library(RcppML)
library(RColorBrewer)
library(ggplot2)
DATADIR= '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers'


lnmf <- function(data, k_wh, k_uv, tol = 1e-4, maxit = 100, L1 = c(0, 0), L2 = c(0, 0), nonneg = TRUE, seed = NULL, mask = NULL, verbose = FALSE ,...) {
  if (length(k_uv) != length(data)) stop("number of ranks specified in 'k_uv' must equal the length of the list of datasets in 'data'")
  if (length(data) == 1) stop("only one dataset was provided, linked NMF is only useful for multiple datasets")
  
  # define initial "h", the linking matrix
  set_pointers <- c(sapply(data, function(x) ncol(x)))
  for (i in 2:length(set_pointers))
    set_pointers[i] <- set_pointers[i] + set_pointers[i - 1]
  set_pointers <- c(0, set_pointers)
  k_pointers <- c(k_wh, k_uv)
  for (i in 2:length(k_pointers))
    k_pointers[i] <- k_pointers[i] + k_pointers[i - 1]
  n_samples <- sum(sapply(data, function(x) ncol(x)))
  link_matrix <- matrix(0, sum(k_wh, k_uv), n_samples)
  link_matrix[1:k_wh, 1:n_samples] <- 1
  
  for (i in 1:length(k_uv))
    link_matrix[(k_pointers[i] + 1):k_pointers[i + 1], (set_pointers[i] + 1):set_pointers[i + 1]] <- 1
  
  # combine data into a single matrix
  if (!all(sapply(data, function(x) class(x)) == class(data[[1]]))) stop("'data' contains items of different classes")
  if (!all(sapply(data, function(x) nrow(x)) == nrow(data[[1]]))) stop("'data' contains items with different numbers of rows")
  data <- do.call(cbind, data)
  
  link_matrix <- as(link_matrix, "ngCMatrix")
  
  model <- nmf(data, nrow(link_matrix), tol = tol, maxit = maxit, L1, L2, nonneg = nonneg, seed = seed, mask = mask, link_h = TRUE, link_matrix_h = link_matrix, sort_model = FALSE, verbose = verbose)
  diag_order_wh <- order(model@d[1:k_wh], decreasing = TRUE)
  w <- model@w[, diag_order_wh]
  rownames(model@h) <- paste0("h", 1:nrow(model@h))
  colnames(w) <- paste0("w", 1:ncol(w))
  u <- v <- h <- d_wh <- d_uv <- list()
  for (i in 1:length(k_uv)) {
    u[[i]] <- as.matrix(model@w[, (k_pointers[i] + 1):k_pointers[i + 1]])
    d_uv[[i]] <- model@d[(k_pointers[i] + 1):k_pointers[i + 1]]
    diag_order_uv <- order(d_uv[[i]], decreasing = TRUE)
    d_uv[[i]] <- d_uv[[i]][diag_order_uv]
    v[[i]] <- as.matrix(model@h[(k_pointers[i] + 1):k_pointers[i + 1], (set_pointers[i] + 1):set_pointers[i + 1]])
    if(ncol(v[[i]]) == 1) v[[i]] <- t(v[[i]])
    if(nrow(v[[i]]) > 1){
      v[[i]] <- v[[i]][diag_order_uv,]
      u[[i]] <- u[[i]][, diag_order_uv]
    }
    h[[i]] <- as.matrix(model@h[diag_order_wh, (set_pointers[i] + 1):set_pointers[i + 1]])
    if(ncol(h[[i]]) == 1) h[[i]] <- t(h[[i]])
    d_wh[[i]] <- model@d[diag_order_wh]
    scale_h <- rowSums(h[[i]])
    h[[i]] <- apply(h[[i]], 2, function(x) x / scale_h)
    d_wh[[i]] <- d_wh[[i]] * scale_h
    names(d_wh[[i]]) <- NULL
    colnames(u[[i]]) <- paste0("u", i, ".", 1:ncol(u[[i]]))
    rownames(v[[i]]) <- paste0("v", i, ".", 1:nrow(v[[i]]))
    rownames(u[[i]]) <- rownames(w)
    colnames(v[[i]]) <- colnames(h[[i]])
  }
  return(new("lnmf", w = w, u = u, v = v, h = h, d_wh = d_wh, d_uv = d_uv, misc = model@misc))
}


#####################################
### read in command line options ####
library(optparse)
option_list = list(
  make_option(c("-i", "--iteration"), type="numeric", default=25, 
              help="integer of which iteration to use", metavar="numeric"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


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


## factor a lot of models and find optimal ranks
tol = 1e-3
param_search = seq(0, .99, .1)
l1_h = param_search[(opt$iteration - 1) / length(param_search) + 1]
l1_w = param_search[(opt$iteration - 1) %% length(param_search) + 1]

save_model_fn = file.path(DATADIR, 'rdas', 'NMF_models', 
                          paste0('DLFPC_PeakMatrix_NMF_k10-100_tol',as.character(tol),
                                 '_l1h',l1_h, '_l1w', l1_w, '.rds'))
if(file.exists(save_model_fn)) break

groups = case_when(grepl('L2|L3|L4', as.character(colData(peakMat)$Celltype2)) ~ 'EXC.Upper', 
                   grepl('L5|L6', as.character(colData(peakMat)$Celltype2)) ~ 'EXC.Lower', 
                   grepl('LAMP|NDNF|VIP', as.character(colData(peakMat)$Celltype2)) ~ 'INH.CGE', 
                   grepl('SST|PV', as.character(colData(peakMat)$Celltype2)) ~ 'INH.MGE', 
                   TRUE ~ 'GLIA')
indList = split(seq(ncol(mat)), groups )
matList = lapply(indList, function(ind) mat[, ind])

lnmf_model <- lnmf(matList, k_wh = 15, k_uv = rep(5, length(indList)),  verbose = TRUE,
                   tol = tol, maxit = 200, L1 = c(.5, .75),
                   L2 = c(.1, .1), seed = 1:10, nonneg = TRUE)
nmf_model <- as(lnmf_model, "nmf")



## get 
celltypes1 = c("EXC","INH_LAMP5","INH_PVALB","INH_SST" ,
               "INH_VIP", "Astro", "Endo", "Microglia","Oligo", "OPC")
celltypes2 = c("L2.CUX2.MEIS2", "L3.CUX2.RORB" , "L4.ALPL" , "L4.TYR" , "L4.5.TBX15" ,
               "L5.PCP4" ,  "L5.POU3F1", "L5.6.NR4A2" , "L6.ITGA8", "L6.NKD1",  "L6.SYT6", 
               "LAMP5",  "NDNF",  "VIP",  "PV.BC" , "PV.ChC" ,"SST", 
               "Astro", "Endo" , "Microglia", "Mural" , "Oligo" , "OPC" )
celltypes2_cols = setNames( c(brewer.pal(5, 'PuBu'), brewer.pal(6, 'Purples'), 
                              brewer.pal(6, 'YlOrRd'), brewer.pal(6, 'Greys')), celltypes2)
pdf('tmp.pdf', height = 8, width = 12)
plot(summary(nmf_model, 
             group_by = factor(peakMat$Celltype2, celltypes2), 
             stat = "mean")) + 
  scale_fill_manual(values = celltypes2_cols)
dev.off()


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
