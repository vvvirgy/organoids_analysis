# prepare data
library(tidyverse)
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')
source('organoids_analysis/R/functions_utils/fit_plots.R')

matching_samples = readxl::read_excel('data/scRNA_samples.xlsx', sheet = 1) %>% 
  as.data.frame() %>% 
  dplyr::mutate(genomics_code = ifelse(genomics_code == 'NA', NA, genomics_code)) %>% 
  dplyr::mutate(fixed_name = ifelse(fixed_name == 'NA', NA, fixed_name))

samples_check = matching_samples %>% 
  dplyr::select(fixed_name, scRNA_sample) %>%
  dplyr::distinct()

new_mapping_experiment = readRDS("data/mapping_samples.rds")

all_samples_dict = full_join(matching_samples, new_mapping_experiment, 
                             by = join_by('genomics_code' == 'genomics_code', 
                                          'fixed_name' == 'fixed_name'), na_matches = 'na')
# saveRDS(all_samples_dict, 'data/samples_names_dictionary.rds')

normalized_res = readRDS('data/normalized_res_pseudobulk_v2.rds')

# genes_cna_status = readRDS('data/karyotypes_full_cohort_cgs.rds')
genes_cna_status = readRDS('data/karyotypes_mutations_all_genes_qc_ccf.rds')
genes_to_check = genes_cna_status$hgnc_symbol %>% unique

# filter genes that are highly fragmented and remove multi-ploidy genes
# filter genes that are highly fragmented
coad_genes = readRDS('data/all_genes_positions.rds')
coad_genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) 

genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::select(hgnc_symbol, chr, from, to)

genes_cna_status = genes_cna_status %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence))

drivers_to_check_correct_samples = filter_fragmented_cnas(genes_cna_status, 
                                                          samples_list = (matching_samples$fixed_name %>% unique), 
                                                          genes_to_check = genes_to_check, 
                                                          genes_position = genes, 
                                                          min_length = 10^6, 
                                                          strategy = 'MOv')
drivers_to_check_correct_samples = drivers_to_check_correct_samples %>% 
  mutate(mutation_multiplicity = ifelse(mut_consequence == 'wild-type', 0, mutation_multiplicity))

# Filter the data to keep only selected genes and prepare the names to be correct 
expression_genes = normalized_res %>%
  as.data.frame() %>% 
  tibble::rownames_to_column('Genes') %>% 
  dplyr::filter(Genes %in% genes_to_check) %>% 
  reshape2::melt() %>% 
  dplyr::rename(PDO = variable) %>% 
  dplyr::group_by(PDO, Genes)
# dplyr::mutate(PDO = gsub('_HSR', 'HSR', PDO)) %>% 
# dplyr::rename(replicate = variable)

transcriptomics_data = drivers_to_check_correct_samples %>%
  dplyr::select(chr, hgnc_symbol, karyotype, sample, is_mutated, mut_consequence, driver_label, CGC_role_COAD, CGC_role_PANCANCER, is_driver_intogen, mutation_multiplicity, VEP.IMPACT) %>% 
  distinct(.keep_all = F) %>% 
  # add correct sample names to map genomics and proteomics
  left_join(., samples_check, by = join_by("sample" == "fixed_name")) %>% 
  # add proteomics data
  full_join(., expression_genes, by = join_by("scRNA_sample" == "PDO", "hgnc_symbol" == "Genes")) %>% 
  # remove NAs and do some filtering
  dplyr::filter(!is.na(scRNA_sample)) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(karyotype)) %>% 
  tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = TRUE) %>% 
  dplyr::mutate(tot_cna = Major+minor) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(karyotype = paste(Major, minor, sep = ':'))

saveRDS(transcriptomics_data, 'data/transcriptomics_data_all_genes_v2.rds')

old_data = readRDS('data/transcriptomics_data_all_genes.rds')

old_p = old_data %>% 
  dplyr::filter(hgnc_symbol == 'KRAS') %>% 
  ggplot(aes(x = tot_cna, 
             y = value, colour = karyotype)) + 
  geom_point() +
  # facet_wrap(hgnc_symbol ~ is_mutated, scales = 'free_x') +
  theme_bw() + 
  scale_color_brewer(palette = 'Set1') + 
  ggtitle('old data')

new = transcriptomics_data %>% 
  dplyr::filter(hgnc_symbol == 'KRAS') %>% 
  ggplot(aes(x = tot_cna, 
             y = value, colour = karyotype)) + 
  geom_point() +
  # facet_wrap(hgnc_symbol ~ is_mutated, scales = 'free_x') +
  theme_bw() + 
  scale_color_brewer(palette = 'Set1') + 
  ggtitle('new data')

old_p / new
