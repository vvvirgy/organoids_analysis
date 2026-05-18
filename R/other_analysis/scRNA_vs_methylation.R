setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(SingleCellExperiment)

# dictionary with the sample ids for RNA and genomic
dict = readRDS('data/full_dict_dna_rna_prot.rds')

# load data as a list and name it 
data_path = 'data/scRNA' # data path

# load seurat objects
data = lapply(list.files(data_path, full.names = T)[1], function(x) {
  readRDS(x)
})

samples = list.files(data_path) %>% 
  gsub('_filtered.rds', '', .) %>% 
  gsub('Sample_', '', .)

names(data) = samples[1]

# set thresholds for filtering
min_lib_size <- 900
min_n_genes <- 250
max_pct_mt <- 22.5
min_cells <- 3

# compute the mitocondrial expression per cell 
data = lapply(data, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  
  return(x)
})

seurat_list = lapply(data %>% names, function(x) {
  
  print(x)
  
  x = data[[x]]
  
  metadata_before_qc = x@meta.data
  
  # detect low quality cells
  is_low_quality <- 
    x$nCount_RNA < min_lib_size |
    x$nFeature_RNA < min_n_genes |
    x$percent.mt > max_pct_mt
  
  # filter low quality cells 
  x <- subset(x, cells = colnames(x)[!is_low_quality])
  
  n_cells <- Matrix::rowSums(x[["RNA"]]$counts > 0)
  kept_genes <- rownames(x)[n_cells > min_cells]
  x = subset(x, features = kept_genes)
  
  return(x)
})
names(seurat_list) = names(data)

# run normalisation
seurat_norm = lapply(seurat_list %>% names, function(s) {
  
  print(s)
  seurat = seurat_list[[s]]
  seurat$sample_id = rep(s, seurat$orig.ident %>% length)
  
  seurat <- seurat %>%
    NormalizeData(
      normalization.method = "LogNormalize",
      scale.factor = 10000
    ) %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% 
    ScaleData(features = rownames(seurat)) %>%
    RunPCA() %>%
    # RunHarmony(group.by.vars = "sample_id", # reduction = "pca", 
    #            # dims = 1:30
    #            ) #%>%
    RunUMAP(dims = 1:30, reduction = "pca")
  return(seurat)
})
names(seurat_norm) = names(seurat_list)

# extract count matrix per sample
count_matrix = lapply(seurat_norm, function(x) {
  x@assays$RNA$data
})


