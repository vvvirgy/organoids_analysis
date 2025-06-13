setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('organoids_analysis/R/functions_utils/glm_function.R')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')

# trying to fit a glm to account for mutational status in the model - mutated, not mutated gene

# preparing the data
transcriptomics_data = transcriptomics_data %>% 
  mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated')))

# prepare the function - per gene!

# perform cross validation
fit = glm_fit(transcriptomics_data,
              response = 'value', 
              model = as.formula('~ tot_cna + mutation_status'))

saveRDS(fit, 'data/glm_fit_transcriptomics.rds')
