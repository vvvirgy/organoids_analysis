.libPaths(new = '/orfeo/cephfs/home/area/vgazziero/R/rstudio_4_4')
setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix.utils)
library(DESeq2)
library(harmony)
# library(RUVSeq)
source('organoids_analysis/R/functions_utils/scRNA_utils.R')

# load data as a list and name it 
data_path = 'data/scRNA'
data = lapply(list.files(data_path, full.names = T, pattern = 'PDO45'), function(x) {
  readRDS(x)
})

samples = list.files(data_path) %>% 
  gsub('_filtered.rds', '', .) %>% 
  gsub('Sample_', '', .)

names(data) = samples

# merge all seurat objects in a single one
alldata = merge(x = data[[1]], y = data[-1], add.cell.ids = samples, project = 'COAD_PDO', merge.data = F, na.rm = F)

# perform data integration -- add batch column to differentiate the two cohorts
alldata$batch = alldata@meta.data %>% 
  dplyr::mutate(batch = case_when(
    str_detect(sample, 'PDO') ~ 'HSR', 
    .default = 'ICR')
    ) %>% 
  dplyr::pull(batch)

# add the mitocondrial expression (even if they are already filtered)
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

# normalize data and run pca + batch correction
alldata = alldata %>% 
  NormalizeData(
    normalization.method = "LogNormalize",
    scale.factor = 10000
  ) %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% 
  ScaleData(features = rownames(alldata)) %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "batch",
             dims = 1:30
             ) %>%
  RunUMAP(dims = 1:30, reduction = "harmony") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6)



