setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('organoids_analysis/R/functions_utils/glm_function.R')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')

all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))
# genes = readRDS('data/genes_to_check.rds')

# trying to fit a glm to account for mutational status in the model - mutated, not mutated gene

# preparing the data
transcriptomics_data = transcriptomics_data %>% 
  dplyr::filter(hgnc_symbol %in% all_genes) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::mutate(tot_cna = factor(tot_cna))

# prepare the function - per gene!

# perform cross validation
set.seed(176613)
fit_model_full = glm_fit(transcriptomics_data,
              response = 'value',
              model = as.formula('~ tot_cna + mutation_status'),
              alphas = 0)
saveRDS(fit_model_full, 'data/glm_fit_ridge_transcriptomics.rds')


fit_model_single = glm_fit(transcriptomics_data,
              response = 'value', 
              model = as.formula('~ tot_cna'), 
              alphas = 0)
saveRDS(fit_model_single, 'data/glm_fit_ridge_transcriptomics_cna_only.rds')
