library(tidyverse)

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds') 

# proteogenomics_data = proteogenomics_data %>% 
#   group_by(proteomics_code, hgnc_symbol) %>% 
#   mutate(mean_intensity = mean(norm_intensity)) %>% 
#   dplyr::select(everything(), -c(replicate, norm_intensity)) %>% 
#   distinct()

rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
                        by = join_by('sample' == 'sample', 
                                     'hgnc_symbol' == 'hgnc_symbol', 
                                     'karyotype' == 'karyotype', 
                                     'Major' == 'Major', 
                                     'minor' == 'minor', 
                                     'tot_cna' == 'tot_cna', 
                                     'chr' == 'chr', 
                                     'is_mutated' == 'is_mutated', 
                                     'mut_consequence' == 'mut_consequence', 
                                     'driver_label' == 'driver_label')) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(!is.na(proteomics_code)) 

sd_genes = rna_vs_prot %>% 
  dplyr::group_by(hgnc_symbol) %>% 
  dplyr::summarise(sd_prot = sd(mean_intensity), sd_rna = sd(value)) %>% 
  dplyr::filter(!is.na(sd_prot)) %>% 
  dplyr::filter(!is.na(sd_rna)) %>% 
  dplyr::filter(sd_prot != 0) %>% 
  dplyr::filter(sd_rna != 0)

rna_vs_prot_corr = rna_vs_prot %>% 
  filter(hgnc_symbol %in% sd_genes$hgnc_symbol)

data = rna_vs_prot_corr %>% 
  dplyr::select(hgnc_symbol, value, mean_intensity) %>% 
  reshape2::melt() %>% 
  dplyr::rename(method = variable) %>% 
  dplyr::rename(expression = value) %>% 
  group_by(hgnc_symbol) %>% 
  group_split()
names(data) = lapply(data, function(s) {s$hgnc_symbol %>% unique}) %>% unlist

test = get_corr_stats(data,
                            cond1 = 'value', 
                            cond2 = 'mean_intensity', 
                            method = 'spearman')
