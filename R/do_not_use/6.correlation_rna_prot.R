rm(list=ls())
.libPaths()
setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(tidyverse)
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/glm_multiresponse.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/classify_muts.R')
source('/orfeo/scratch/cdslab/vgazziero/TLS/CLL_Cu/copper/CuRes/R/utils/expression/data_manipulation/corr_functions.R')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds') %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(rna_expression))

transcriptomics_data = classify_mutations(transcriptomics_data) 

proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds') %>%
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>%
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>%
  dplyr::rename(protein_expression = mean_intensity) %>%
  dplyr::filter(!is.na(protein_expression))
proteogenomics_data = classify_mutations(proteogenomics_data)

all_data = full_join(proteogenomics_data, transcriptomics_data) %>% 
  dplyr::filter(!is.na(rna_expression)) %>% 
  dplyr::filter(!is.na(protein_expression))
# all_data = all_data %>% 
#    %>% 
#   mutate(sample_type = case_when(str_detect(sample, 'ATN') ~ 'ATN', .default = 'Cu')) %>% 
#   mutate(method = case_when(str_detect(sample, 'PROT') ~ 'PROT', .default = 'RNA')) 

all_data = all_data %>% 
  group_by(hgnc_symbol) %>% 
  group_split()
names(all_data) = lapply(all_data, function(x) {
  x$hgnc_symbol %>% unique
}) %>% unlist

all_data = lapply(all_data, function(x) {
  x %>% 
    mutate(across(matches('expression'), as.numeric)) %>% 
    dplyr::select(hgnc_symbol, sample, rna_expression, protein_expression) %>% 
    # pivot_longer(names_to = 'sample',
    #              values_to = 'expression', 
    #              cols = c('rna_expression', 'protein_expression'))
    reshape2::melt() %>% 
    dplyr::rename(expression = value) %>% 
    dplyr::rename(method = variable)
})

correlations = get_corr_stats(corr = all_data, method = 'spearman', cond1 = 'protein_expression', cond2 = 'rna_expression')
correlations = correlations %>% 
  dplyr::relocate(gene_condition, .before = pValue)

correlations %>% 
  mutate(cor_type = ifelse(correlation > 0, 'positive', 'negative')) %>% 
  mutate(cor_type = ifelse(qValue > 0.05, NA, cor_type)) %>% 
  ggplot(aes(y = correlation, x = -log(qValue, 10), color = cor_type)) + 
  geom_point() + 
  theme_bw() +
  scale_color_manual(values = c('positive' = 'goldenrod2', 'negative' = 'steelblue4')) + 
  theme(legend.position = 'bottom') 
ggsave('res/corr_rna_prot_spearman.png', width = 5, height = 4, bg = 'white')
 
