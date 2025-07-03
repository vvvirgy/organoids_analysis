# prepare proteomics data for fitting 
library(tidyverse)
source('organoids_analysis/R/functions_utils/proteomics_utils.R')
# source('organoids_analysis/R/functions_utils/fit_plots.R')
# source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

# get correct sample names
new_mapping_experiment = readRDS("data/mapping_samples.rds")

# get the pair proteomics-genomics correct sample name
samples_check = new_mapping_experiment %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
  dplyr::distinct()

data = readxl::read_excel('data/Results_Organoids_NoIsoforms.xlsx', sheet = 1) %>% 
  as.data.frame()

data[which(is.na(data[,'PG.Genes'])), 'PG.Genes'] = 'Unknown'
raw_data = data %>% 
  dplyr::select(everything(), -PG.ProteinGroups) %>% 
  tibble::column_to_rownames('PG.Genes') %>% 
  dplyr::mutate(across(where(is.character), as.numeric)) %>% 
  dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

raw_data_melted = reshape2::melt(raw_data)

proteomics_data_norm = compute_tic(raw_data = raw_data, 
                                   gene_column = 'Genes', 
                                   log2 = TRUE, 
                                   filter_NA = TRUE) 

saveRDS(proteomics_data_norm, 'data/proteomics_tic_normalised.rds')

