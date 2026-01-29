library(tidyverse)
source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

dict = readRDS('data/full_dict_dna_rna_prot.rds')
genes_cna_status = readRDS('data/processed_data/karyotypes_mutations_all_genes_qc_ccf.rds')

samples_check = dict %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
  dplyr::distinct()

genes_to_check = genes_cna_status$hgnc_symbol %>% unique

coad_genes = readRDS('data/all_genes_positions_info.rds')
genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::select(hgnc_symbol, chr, from, to)

genes_cna_status = genes_cna_status %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence))

# removing mutations falling on the same segment + multihit mutations + selecting the segments with the maximum overlap over the gene + segments with a minimun length
genes_filtered_test = filter_fragmented_cnas(genes_cna_status, 
                                        samples_list = (dict$fixed_name %>% unique), 
                                        genes_to_check = genes_to_check, 
                                        min_length = 10^6, 
                                        ccf_thr = .8, 
                                        filter_multihit = FALSE, 
                                        genes_position = genes, 
                                        strategy = 'MOv')

genes_filtered_v2 = genes_filtered_test_v2 %>% 
  select(sample, hgnc_symbol, karyotype, mut_consequence) %>% 
  ungroup() %>% 
  distinct()

genes_filtered_v2 %>% 
  group_by(hgnc_symbol, sample) %>% 
  count() %>% 
  arrange(desc(n))

saveRDS(genes_filtered_v2, 'data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')
saveRDS(genes_filtered_test_v2, 'data/processed_data/genes_filtered_karyo_mut_status_full_filt_ccf_08.rds')
