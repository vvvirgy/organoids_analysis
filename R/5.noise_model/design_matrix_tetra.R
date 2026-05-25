library(tidyverse)
library(dplyr)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')

args <- commandArgs(trailingOnly = TRUE)

# Default values if not provided
# sf_method  <- 'psinorm'
# use_stable <- FALSE
# 
# cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))
# 
# rm(list = setdiff(ls(), c("sf_method", "use_stable")))

source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/constants.R")
source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/make_groups.R")

IMG_PATH = file.path(res_path, 'noise_model')
RES_PATH = file.path(data_path, 'noise_model')

dir.create(IMG_PATH, recursive = T)
dir.create(RES_PATH, recursive = T)

# using metadata with removed genes that are highly fragmented

genes_cna_status = readRDS(file.path(data_path, 'processed_data/dna/genes_filtered_karyo_mut_status_filt_ccf_08.rds'))

# selecting only some genes to run the dge --> removing genes without or with low number of diploids

diploid_genes = genes_cna_status %>% 
  filter(karyotype %in% '1:1') %>% 
  group_by(hgnc_symbol) %>% 
  dplyr::count() %>% 
  filter(n >= 5) %>% 
  pull(hgnc_symbol) %>% 
  unique 

design_matrix_samples_groups = make_groups(genes_cna_status, diploid_genes, karyo = '1:1', n = 100)
saveRDS(design_matrix_samples_groups, file.path(RES_PATH, 'design_matrix_diploid.rds'))

design_matrix_samples_groups = readRDS(file.path(RES_PATH, 'design_matrix_diploid.rds'))

filtered_design_matrix = lapply(design_matrix_samples_groups, function(df) {
  good_iters = df %>%
    dplyr::mutate(sample_group = paste0(sample, ":", group)) %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(id = paste(sample_group, collapse = ","), .groups = "drop") %>%
    dplyr::group_by(id) %>%
    dplyr::filter(iteration == min(iteration)) %>%
    dplyr::pull(iteration)
  
  df %>% dplyr::filter(iteration %in% good_iters)
})

# names(filtered_design_matrix) = names(design_matrix)
saveRDS(filtered_design_matrix, file.path(RES_PATH, 'filtered_design_matrix_diploids.rds'))
