setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('organoids_analysis/R/functions_utils/glm_function.R')

proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')

all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))

# rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
#                         by = join_by('sample' == 'sample', 
#                                      'hgnc_symbol' == 'hgnc_symbol', 
#                                      'karyotype' == 'karyotype', 
#                                      'Major' == 'Major', 
#                                      'minor' == 'minor', 
#                                      'tot_cna' == 'tot_cna', 
#                                      'chr' == 'chr', 
#                                      'mut_consequence' == 'mut_consequence', 
#                                      'driver_label' == 'driver_label', 
#                                      'is_mutated' == 'is_mutated')) %>% 
#   dplyr::filter(!is.na(value)) %>% 
#   dplyr::filter(!is.na(proteomics_code)) 


# genes = readRDS('data/genes_to_check.rds')

# trying to fit a glm to account for mutational status in the model - mutated, not mutated gene

# preparing the data
proteogenomics_data = proteogenomics_data %>% 
  dplyr::filter(hgnc_symbol %in% all_genes) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) 

# prepare the function - per gene!

# perform cross validation
set.seed(176613)
# fit_model_full = glm_fit(proteogenomics_data,
#               response = 'mean_intensity',
#               model = as.formula('~ tot_cna + mutation_status'), 
#               alphas = 0)
# 
# saveRDS(fit_model_full, 'data/glm_fit_ridge_proteomics.rds')

fit_model_single = glm_fit(proteogenomics_data,
              response = 'mean_intensity',
              model = as.formula('~ tot_cna'), 
              alphas = 0)

saveRDS(fit_model_single, 'data/glm_fit_ridge_proteomics_cna_only.rds')
