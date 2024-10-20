library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(rtracklayer)
library(ArchR)
library(future)
library(harmony)

plan("multiprocess", workers = 8)


ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#loading in annotations
annotations = readRDS(file = "/projects/pfenninggroup/allen_human_brain/allen_human_mca/21_08_11/neuron_meta.rds")

#loading in reference seurat object and taking a look at it
file = '/projects/pfenninggroup/singleCell/allen_human_brain/allen_human_mca/data/allen_human_mca.seurat.rds'
obj = readRDS(file)
meta2 = obj[[]]
head(meta2,1)

#adding the labels to the reference Seurat object
rename_vec = annotations %>% dplyr::select(sample_name, combined_subclass_label) %>% deframe()
obj$combined_subclass_label = rename_vec[obj$sample_name]
obj$combined_subclass_label = make.names(obj$combined_subclass_label)

#need to relabel the NAs in combined_subclass_label, since a lot of them could be glial cells 
#that didn't get labeled correctly. So, relabel them with the celltype in subclass_label, which is accurate
#this works because we've verified that every NA in combined_subclass_label belongs to a non-neuronal cell type!
obj$subclass.id <- with(obj, ifelse(obj$combined_subclass_label == 'NA.', obj$subclass_label,obj$combined_subclass_label))
#obj$subclass.id = obj[[]][obj$combined_subclass_label=='NA.' & obj$class_label == 'Non-neuronal']
table(obj$subclass.id)

# #read in annotations (if necessary)
# obj = readRDS('/projects/pfenninggroup/allen_human_brain/data/human_annotations.rds')

# add in mitochondrial percentage in meta data
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
# run SCTransform and other dim reduction procedures
obj <- SCTransform(obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

#saving as a rds object to be used later
save_seurat_fn = '/projects/pfenninggroup/allen_human_brain/data/6_28_human_annotations.rds'
saveRDS(obj, save_seurat_fn)

# #appending annotations to mouse object
# #load in mouse object
# mouse = readRDS('/projects/pfenninggroup/singleCell/allen_mouse_brain/smartseq_data/aibs_mouse_ctx-hpf_smart-seq_Seurat.ss.rda')
# annotations = readRDS('/projects/pfenninggroup/allen_human_brain/data/human_annotations.rds')
# 
# ###########################################################
# ## 0) # Basic function to convert human to mouse gene names
# convertHumanGeneList <- function(x){
#   
#   require("biomaRt")
#   human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#   mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#   
#   genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
#                    values = x , mart = human, attributesL = c("mgi_symbol"), 
#                    martL = mouse, uniqueRows=T)
#   
#   # get 1-1 gene orthologs by ensembl
#   oneToOne = split(genesV2$HGNC.symbol, genesV2$MGI.symbol) %>% 
#     lengths()
#   oneToOne = names(oneToOne)[which(oneToOne ==1)]
#   genesV2 = genesV2[genesV2$MGI.symbol %in% oneToOne,]
#   
#   # Print the first 6 genes found to the screen
#   print(head(genesV2))
#   return(genesV2)
# }
# 
# #######################################################
# ## 1) load in the Seurat objects  
# ## create new count object changing human gene names to 1-1 mouse ortholog
# rm_to_mm_genes <- convertHumanGeneList(rownames(mouse))
# 
# #make sparse matrix easier to work with
# # dcg_format = mouse@assays$RNA@counts
# # as_matrix <- function(mat){
# #   
# #   tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
# #   
# #   row_pos <- mat@i+1
# #   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
# #   val <- mat@x
# #   
# #   for (i in seq_along(val)){
# #     tmp[row_pos[i],col_pos[i]] <- val[i]
# #   }
# #   
# #   row.names(tmp) <- mat@Dimnames[[1]]
# #   colnames(tmp) <- mat@Dimnames[[2]]
# #   return(tmp)
# # }
# # 
# # mouse = as_matrix(dcg_format)
# 
# counts_rna = round(2^mouse@assays$RNA@counts)
# counts_rna = counts_rna[rownames(counts_rna) %in% rm_to_mm_genes$HGNC.symbol,]
# rownames(counts_rna) = rm_to_mm_genes$MGI.symbol[match(rownames(counts_rna), 
#                                                        rm_to_mm_genes$HGNC.symbol)]
# mouse[["RNA"]] <- CreateAssayObject(counts = counts_rna)
# 
# ## save this version of the dataset
# SaveH5Seurat(mouse, misc = F, tools = F, overwrite = T,
#              filename = '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/rdas/allen_human_labels_integrated.h5seurat')
# 
# sce_rm = as.SingleCellExperiment(mouse)
# saveRDS(sce_rm, '/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/rdas/allen_human_labels_integrated.sce.rds')
# 






















