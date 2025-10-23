# prepare proteomics data for fitting 
library(tidyverse)
source('organoids_analysis/R/functions_utils/proteomics_utils.R')
source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

# load proteomics data after normalization
proteomics_data_norm = readRDS('data/proteomics_normalized.rds')

# get correct sample names
new_mapping_experiment = readRDS("data/mapping_samples.rds")

# get the pair proteomics-genomics correct sample name
samples_check = new_mapping_experiment %>% 
  dplyr::select(proteomics_code, fixed_name) %>%
  dplyr::distinct()

# genes_cna_status = readRDS('data/karyotypes_full_cohort_cgs.rds')

# load genes with cna associated (only genes in common with transcriptomics)
genes_cna_status = readRDS('data/karyotypes_mutations_all_genes_qc_ccf.rds') #%>% 
  # dplyr::bind_rows() 
genes_to_check = genes_cna_status$hgnc_symbol %>% unique

# proteomics_data_norm %>% 
#   dplyr::mutate(sample = gsub('_a|_b', '', variable)) %>% 
#   # dplyr::filter(sample == '2L') %>% 
#   ggplot(aes(norm_intensity, colour = sample, y = after_stat(count))) + 
#   # geom_histogram(binwidth = 0.1, position = 'identity') +
#   geom_density(alpha = 0.1) +
#   theme_bw()

# raw_data_melted %>% 
#   dplyr::mutate(sample = gsub('_a|_b', '', variable)) %>% 
#   dplyr::mutate(value = log(value, base = 2)) %>% 
#   # dplyr::filter(sample == '2L') %>% 
#   ggplot(aes(value, colour = sample, y = after_stat(count))) + 
#   # geom_histogram(binwidth = 0.1, position = 'identity') +
#   geom_density(alpha = 0.1) +
#   theme_bw()

# filter genes that are highly fragmented
# coad_genes = readRDS('data/all_genes_positions.rds')
coad_genes = readRDS('data/all_genes_positions_info.rds')
genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::select(hgnc_symbol, chr, from, to)

genes_cna_status = genes_cna_status %>% 
  dplyr::mutate(mut_consequence = ifelse(is.na(mut_consequence), 'wild-type', mut_consequence))
# removing mutations falling on the same segment + multihit mutations + selecting the segments with the maximum overlap over the gene 
drivers_to_check_correct_samples = filter_fragmented_cnas(genes_cna_status, 
                                                          samples_list = (samples_check$fixed_name %>% unique), 
                                                          genes_to_check = genes_to_check, 
                                                          min_length = 10^6, 
                                                          genes_position = genes, 
                                                          strategy = 'MOv')
drivers_to_check_correct_samples = drivers_to_check_correct_samples %>% 
  mutate(mutation_multiplicity = ifelse(mut_consequence == 'wild-type', 0, mutation_multiplicity))

# Filter the data to keep only selected genes and prepare the names to be correct 
proteomic_genes = proteomics_data_norm %>%
  dplyr::filter(PG.Genes %in% genes_to_check) %>% 
  # dplyr::mutate(PDO = gsub("_a|_b", "", variable)) %>% 
  dplyr::group_by(PDO, PG.Genes) %>% 
  # dplyr::mutate(PDO = gsub('_HSR', 'HSR', PDO)) %>% 
  dplyr::rename(replicate = sample) %>% 
  group_by(PDO, PG.Genes) %>%
  mutate(mean_intensity = mean(transformedIntensity)) %>%
  dplyr::select(everything(), -c(replicate, transformedIntensity, log2_Intensity, sampleName, Intensity)) %>%
  distinct()

proteogenomics_data = drivers_to_check_correct_samples %>%
  dplyr::select(chr, hgnc_symbol, karyotype, sample, is_mutated, mut_consequence, driver_label, CGC_role_COAD, CGC_role_PANCANCER, is_driver_intogen, mutation_multiplicity, VEP.IMPACT) %>% 
  distinct(.keep_all = F) %>% 
  # add correct sample names to map genomics and proteomics
  left_join(., samples_check, by = join_by("sample" == "fixed_name")) %>% 
  # add proteomics data
  full_join(., proteomic_genes, by = join_by("proteomics_code" == "PDO", "hgnc_symbol" == "PG.Genes")) %>% 
  # remove NAs and do some filtering
  dplyr::filter(!is.na(proteomics_code)) %>% 
  distinct() %>% 
  dplyr::filter(!is.na(karyotype)) %>% 
  tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':', convert = TRUE) %>% 
  dplyr::mutate(tot_cna = Major+minor) %>% 
  dplyr::mutate(karyotype = paste(Major, minor, sep = ':')) 

# proteogenomics_data_low_CNAs = proteogenomics_data %>% 
#   dplyr::filter(tot_cna <= 6)

# proteogenomics_data_v2 = proteogenomics_data %>%
#   group_by(proteomics_code, hgnc_symbol) %>%
#   mutate(mean_intensity = mean(norm_intensity)) %>%
#   dplyr::select(everything(), -c(replicate, norm_intensity)) %>%
#   distinct()

saveRDS(proteogenomics_data, 'data/proteogenomics_data_all_genes_new_norm.rds')

p_old = proteogenomics_data_v2 %>% 
  filter(hgnc_symbol == 'PIK3CA') %>% 
  ggplot(aes(x = tot_cna, 
             y = mean_intensity, colour = karyotype)) + 
  geom_point(size = 8) +
  # facet_wrap(hgnc_symbol ~ is_mutated, scales = 'free_x') +
  theme_bw() + 
  scale_color_brewer(palette = 'Set1') + 
  ggtitle('old normalization')

p_new = proteogenomics_data %>% 
  filter(hgnc_symbol == 'PIK3CA') %>% 
  ggplot(aes(x = tot_cna, 
             y = mean_intensity, colour = karyotype)) + 
  geom_point(size = 8) +
  # facet_wrap(hgnc_symbol ~ is_mutated, scales = 'free_x') +
  theme_bw() + 
  scale_color_brewer(palette = 'Set1') + 
  ggtitle('new normalization')
  