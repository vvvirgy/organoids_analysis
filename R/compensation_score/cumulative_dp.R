rm(list=ls())
library(CNAqc)
library(tidyverse)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/preprocessing/get_cnas_muts.R')

cnas = readRDS('data/cnaqc_v2/cnas_list_v2.rds')
genes_cna_status = readRDS('data/processed_data/genes_filtered_karyo_mut_status_full_filt_ccf_08.rds')
karyos = c('1:1', '2:1', '1:0', '2:0', '2:2')

multi_omics = readRDS('data/lfc_prot_and_rna_bind.rds')

ggenes = multi_omics %>% 
  group_by(name, karyotype) %>% 
  filter(n() > 1) %>% 
  pull(name) %>% 
  unique
  
genes_cna_status_filt = genes_cna_status %>% 
  filter(hgnc_symbol %in% ggenes) %>% 
  select(hgnc_symbol, sample, karyotype, starts_with('segment'), chr) %>% 
  distinct() %>% 
  split(.$sample)

cnas_filt = cnas[names(genes_cna_status_filt)]

cum_dp = lapply(names(cnas_filt), function(x) {
  
  # print(x)
  
  mut = cnas[[x]]$mutations
  
  cumulative_dp = mut %>% 
    select(VEP.SYMBOL, DP, segment_id) %>% 
    group_by(segment_id) %>% 
    summarise(segment_mean_DP = mean(DP, na.rm = T),
              segment_median_DP = median(DP, na.rm = T), 
              segment_genes = list(VEP.SYMBOL), 
              segment_DP = list(DP)
              # segment_genes = paste(VEP.SYMBOL, collapse = '/'), 
              # segment_DP = paste(DP, collapse = '/'), 
              ) %>% 
    tidyr::separate(segment_id, into = c('chr', 'from', 'to', 'Major', 'minor'), sep = ':', remove = F, convert = T) %>% 
    relative_to_abs_coords(., ref = 'GRCh38') 
  
  gg = genes_cna_status_filt[[x]]
    
  cumulative_dp %>% 
    full_join(., gg, by = join_by(
      'chr' == 'chr', 
      'from' == 'segment_from', 
      'to' == 'segment_to'
    )) %>% 
    filter(!is.na(hgnc_symbol)) %>% 
    filter(!is.na(segment_id))
} )

cum_dp_red = cum_dp %>% 
  bind_rows() %>%
  select(chr, from, to, karyotype, segment_mean_DP, segment_median_DP, sample) %>% 
  distinct()

cum_dp_summarised = cum_dp %>% 
  bind_rows()

saveRDS(cum_dp_summarised, 'data/cumulative_dp.rds')

cum_dp_summarised = lapply(cum_dp, function(x) {
  x %>%
    # bind_rows() %>% 
    separate_rows(c(segment_genes, segment_DP), sep = '/', convert = T) %>% 
    group_by(hgnc_symbol) %>% 
    mutate(
      segment_genes = list(segment_genes),
      segment_DP    = list(segment_DP)
    ) %>% 
    distinct()
}) %>% bind_rows()

# cum_dp_red %>% 
#   ggplot(aes(
#     x = karyotype, 
#     y = segment_mean_DP
#   )) + 
#   geom_boxplot() + 
#   theme_bw()



