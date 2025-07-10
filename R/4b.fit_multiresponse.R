setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('organoids_analysis/R/functions_utils/glm_function.R')

proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')

all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))

rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
                        by = join_by('sample' == 'sample', 
                                     'hgnc_symbol' == 'hgnc_symbol', 
                                     'karyotype' == 'karyotype', 
                                     'Major' == 'Major', 
                                     'minor' == 'minor', 
                                     'tot_cna' == 'tot_cna', 
                                     'chr' == 'chr', 
                                     'mut_consequence' == 'mut_consequence', 
                                     'driver_label' == 'driver_label', 
                                     'is_mutated' == 'is_mutated', 
                                     'is_driver_intogen' == 'is_driver_intogen', 
                                     'CGC_role_COAD' == 'CGC_role_COAD', 
                                     'CGC_role_PANCANCER' == 'CGC_role_PANCANCER')) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(!is.na(proteomics_code)) %>% 
  dplyr::filter(hgnc_symbol %in% all_genes) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(protein_expression = mean_intensity) %>% 
  dplyr::rename(rna_expression = value)

set.seed(176613)

fit_multiresponse_model = glm_fit_multiresponse(rna_vs_prot,
                           response = c('protein_expression', 'rna_expression'),
                           model = as.formula('~ tot_cna + mutation_status'), 
                           alphas = 1, 
                           lambda = 'lambda.min')
saveRDS(fit_multiresponse_model, 'data/multiresponse_glm_lasso_v2.rds')
