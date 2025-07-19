rm(list=ls())
.libPaths()
setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
library(glmnet)
library(tidyverse)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/glm_multiresponse.R')

proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')

all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))

rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
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
                                     'CGC_role_PANCANCER' == 'CGC_role_PANCANCER', 
                                     'mutation_multiplicity' == 'mutation_multiplicity', 
                                     'VEP.IMPACT' == 'VEP.IMPACT')) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(!is.na(proteomics_code)) %>% 
  dplyr::filter(hgnc_symbol %in% all_genes) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(protein_expression = mean_intensity) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(protein_expression)) %>% 
  dplyr::filter(!is.na(rna_expression))

# rna_vs_prot %>% 
#   filter(sample == '62_VIII', hgnc_symbol == 'TTN') %>% 
#   count(mut_consequence, mutation_multiplicity, karyotype)
# 
# rna_vs_prot %>% 
#   group_by(sample, hgnc_symbol, tot_cna) %>%
#   summarise(n = n()) %>% 
#   # filter(n  > 1 ) %>% 
#   dim
rna_vs_prot = rna_vs_prot %>% 
  dplyr::mutate(
    category = case_when(
      mut_consequence %in% c('upstream_gene_variant', 'downstream_gene_variant', 
                             'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 
                             'synonymous_variant') ~ 'low_effect', 
      mut_consequence %in% c("stop_gained", "start_lost", "stop_lost") ~ 'truncating',
      mut_consequence %in% c("missense_variant", "splice_region_variant", "splice_acceptor_variant", 
                             "non_coding_transcript_exon_variant", "frameshift_variant", "splice_donor_variant",
                             'inframe_insertion', 'inframe_deletion') ~ 'alterating', 
      mut_consequence == 'wild-type' ~ 'wild-type'
    )
  )

rna_vs_prot = rna_vs_prot %>% 
  dplyr::filter(!is.na(mutation_multiplicity)) %>% 
  dplyr::mutate(
    n_low = case_when(
      category == 'low_effect' ~ mutation_multiplicity, 
      .default = 0
    ), 
    n_alt = case_when(
      category == 'alterating' ~ mutation_multiplicity, 
      .default = 0
    ), 
    n_trunc = case_when(
      category == 'truncating' ~ mutation_multiplicity, 
      .default = 0
    ), 
    n_wt = case_when(
      category == 'wild-type' ~ tot_cna, 
      .default = tot_cna - mutation_multiplicity
    )
  ) 

saveRDS(rna_vs_prot, 'data/rna_vs_prot_data.rds')

# rna_vs_prot %>% 
#   filter(category == 'low_effect') %>% 
#   view
# 

# APC has two mutations that induce a stop codon and a third one that is an intronic mutation 
# view(rna_vs_prot %>% filter(sample == 'A4164_1004', hgnc_symbol == 'APC'))



# fit_multiresponse_model = glm_fit_multiresponse(rna_vs_prot,
#                            response = c('protein_expression', 'rna_expression'),
#                            model = as.formula('~ tot_cna + mutation_status'), 
#                            alphas = 1, 
#                            lambda = 'lambda.min')
# saveRDS(fit_multiresponse_model, 'data/multiresponse_glm_lasso_v2.rds')

# new model fitted!
set.seed(176613)
fit_multiresponse_model = glm_fit_multiresponse(rna_vs_prot,
                                                nobs = 4,
                                                response = c('protein_expression', 'rna_expression'),
                                                model = as.formula('~ n_low + n_alt + n_trunc + n_wt'), 
                                                alphas = 1, 
                                                lambda = 'lambda.min')

saveRDS(fit_multiresponse_model, 'data/multiresponse_glm_lasso_multiplicity_v1.rds')

# fit a second time
# set.seed(9016869)
# 
# # new model fitted!
# fit_multiresponse_model = glm_fit_multiresponse(rna_vs_prot,
#                                                 nobs = 4,
#                                                 response = c('protein_expression', 'rna_expression'),
#                                                 model = as.formula('~ n_low + n_alt + n_trunc + n_wt'), 
#                                                 alphas = 1, 
#                                                 lambda = 'lambda.min')
# 
# saveRDS(fit_multiresponse_model, 'data/multiresponse_glm_lasso_multiplicity_v2.rds')


