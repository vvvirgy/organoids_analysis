setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('organoids_analysis/R/functions_utils/glm_function.R')

proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')

# trying to fit a glm to account for mutational status in the model - mutated, not mutated gene

# preparing the data
proteogenomics_data = proteogenomics_data %>% 
  mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated')))

# prepare the function - per gene!

# perform cross validation
fit = glm_fit(proteogenomics_data, 
              response = 'mean_intensity', 
              model = as.formula('~ tot_cna + mutation_status'))

saveRDS(fit, 'data/glm_fit_proteomics.rds')
