rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')

library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(mosaic)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/rna_vs_proteomics/utils.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/classification/utils.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/test_belonging.R')

set.seed(6792820)

Sys.time()

print('loading data')
# consider synonimous mutations as wild-type
rna = readRDS('data/transcriptomics_data_all_genes_v4.rds') %>% 
  dplyr::mutate(is_mutated = ifelse(mut_consequence == 'synonymous_variant', FALSE, is_mutated)) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(rna_expression))

# consider synonimous mutations as wild-type
prot = readRDS('data/proteogenomics_data_all_genes_new_norm_v3.rds') %>%
  dplyr::mutate(is_mutated = ifelse(mut_consequence == 'synonymous_variant', FALSE, is_mutated)) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>%
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>%
  dplyr::rename(protein_expression = mean_intensity) %>%
  dplyr::filter(!is.na(protein_expression))

rna_expr_by_alleles = join_allele_expr(rna, column = 'rna_expression', which = 'Wild-type', use_2n = TRUE) 

prot_expr_by_alleles = join_allele_expr(prot, column = 'protein_expression', which = 'Wild-type', use_2n = TRUE) 

Sys.time()
print('launching bootraps')

# test if any element from the distribution of non altered belong to the distribution of 2nwt samples
test_expr_rna = test_belonging_alterations(rna_expr_by_alleles, pth = 0.1, n_bootstrap = 10000)
Sys.time()
print('Finished RNA bootstrap and test')

test_expr_prot = test_belonging_alterations(prot_expr_by_alleles, pth = 0.1, n_bootstrap = 10000)
Sys.time()
print('Finished protein bootstrap and test')
# # classify each point --> LATER
# test_expr_rna = classify_elements(test_expr_rna) %>% 
#   full_join(., rna_expr_by_alleles) %>% 
#   mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
#   mutate(cls_dosage = gsub('Dose', 'Dosage', cls_dosage)) %>% 
#   rename(expression = rna_expression)
# 
# test_expr_prot = classify_elements(test_expr_prot) %>% 
#   full_join(., prot_expr_by_alleles) %>% 
#   mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
#   mutate(cls_dosage = gsub('Dose', 'Dosage', cls_dosage)) %>% 
#   rename(expression = protein_expression)

saveRDS(test_expr_rna, 'data/test_expr_rna_bootstrap_v3.rds')
saveRDS(test_expr_prot, 'data/test_expr_prot_bootstrap_v3.rds')