.libPaths()
rm(list=ls())
library(tidyverse)

source('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_Cu/copper/CuRes/R/utils/expression/data_manipulation/corr_functions.R')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes_new_norm_v4.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes_v5.rds')
all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))

# proteogenomics_data = proteogenomics_data %>% 
#   group_by(proteomics_code, hgnc_symbol) %>% 
#   mutate(mean_intensity = mean(norm_intensity)) %>% 
#   dplyr::select(everything(), -c(replicate, norm_intensity)) %>% 
#   distinct()

rna_vs_prot = inner_join(transcriptomics_data, proteogenomics_data, 
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
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(protein_expression)) %>% 
  dplyr::filter(!is.na(rna_expression))
  

sd_genes = rna_vs_prot %>% 
  dplyr::group_by(hgnc_symbol) %>% 
  dplyr::summarise(sd_prot = sd(protein_expression), sd_rna = sd(rna_expression)) %>% 
  dplyr::filter(!is.na(sd_prot)) %>% 
  dplyr::filter(!is.na(sd_rna)) %>% 
  dplyr::filter(sd_prot != 0) %>% 
  dplyr::filter(sd_rna != 0)

rna_vs_prot_corr = rna_vs_prot %>% 
  filter(hgnc_symbol %in% sd_genes$hgnc_symbol)

data = rna_vs_prot_corr %>% 
  dplyr::select(hgnc_symbol, rna_expression, protein_expression) %>% 
  reshape2::melt() %>% 
  dplyr::rename(method = variable) %>% 
  dplyr::rename(expression = value) %>% 
  group_by(hgnc_symbol) %>% 
  group_split()
names(data) = lapply(data, function(s) {s$hgnc_symbol %>% unique}) %>% unlist

corr = get_corr_stats(data,
                            cond1 = 'rna_expression', 
                            cond2 = 'protein_expression', 
                            method = 'spearman')
saveRDS(corr, 'data/correlations_protein_rna.rds')

corr = corr %>% 
  mutate(corr_type = case_when(
    (qValue <= 0.05 & correlation > 0) ~ 'positive', 
    (qValue <= 0.05 & correlation < 0) ~ 'negative', 
    .default = 'ns'
  ))

corr %>% 
  ggplot(aes(x = -log(qValue, 10), y = correlation, color = corr_type)) +
  geom_point(alpha = 0.7) + 
  theme_bw() + 
  labs(y = 'Correlation (Spearmann)',
       x = '-log10(pValue)') + 
  scale_color_manual(values = c('positive' = '#3291B6', 'negative' = '#DE1A58', 'ns' = 'grey45')) + 
  guides(color = guide_legend(title = 'Correlation')) + 
  theme(legend.position = 'bottom')
ggsave('res/correlation_spearmann.png', width = 6, height = 5, dpi = 300)

corr %>% 
  group_by(corr_type) %>% 
  count()
