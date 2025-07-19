library(tidyverse)

source('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_Cu/copper/CuRes/R/utils/expression/data_manipulation/corr_functions.R')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
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

corr %>% 
  ggplot(aes(y = -log(qValue, 10), x = correlation)) +
  geom_point()
