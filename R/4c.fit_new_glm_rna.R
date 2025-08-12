rm(list=ls())
.libPaths()
setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")
# library(glmnet)
library(tidyverse)
library(parallel)
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/glm_multiresponse.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/classify_muts.R')

# proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds') %>% 
#   dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
#   dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
#   dplyr::rename(protein_expression = mean_intensity) %>% 
#   dplyr::filter(!is.na(protein_expression)) 

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds') %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(rna_expression))

transcriptomics_data = classify_mutations(transcriptomics_data) 
# proteogenomics_data = classify_mutations(proteogenomics_data) 

# all_data = list(
#   'RNA' = transcriptomics_data, 
#   'Protein' = proteogenomics_data
# )

response = grep('expression', colnames(transcriptomics_data), value = T)
model = as.formula(paste0(response, '~ n_low + n_alt + n_trunc + n_wt'))

glm_fit_v2 = mclapply(transcriptomics_data$hgnc_symbol %>% unique, function(g) {
  print(g)
  df = transcriptomics_data %>% dplyr::filter(hgnc_symbol == g)
  res = glm(formula = model, data =  df)
  return(res)
}, mc.cores = 4, mc.preschedule = T)
names(glm_fit_v2) = transcriptomics_data$hgnc_symbol %>% unique

saveRDS(glm_fit_v2, 'data/glm_fit_v2_rna.rds')

# glm_fit_v2 = lapply(all_data, function(x) {
#   response = grep('expression', colnames(x), value = T)
#   model = as.formula(paste0(response, '~ n_low + n_alt + n_trunc + n_wt'))
#   
#   gg = lapply(x$hgnc_symbol %>% unique, function(g) {
#     print(g)
#     df = x %>% dplyr::filter(hgnc_symbol == g)
#     res = glm(formula = model, data =  x)
#     return(res)
#   })
#   names(gg) = x$hgnc_symbol %>% unique
#   return(gg)
# })

# all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))
# 
# rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
#                         by = join_by('sample' == 'sample', 
#                                      'hgnc_symbol' == 'hgnc_symbol', 
#                                      'karyotype' == 'karyotype', 
#                                      'Major' == 'Major', 
#                                      'minor' == 'minor', 
#                                      'tot_cna' == 'tot_cna', 
#                                      'chr' == 'chr', 
#                                      'mut_consequence' == 'mut_consequence', 
#                                      'driver_label' == 'driver_label', 
#                                      'is_mutated' == 'is_mutated', 
#                                      'is_driver_intogen' == 'is_driver_intogen', 
#                                      'CGC_role_COAD' == 'CGC_role_COAD', 
#                                      'CGC_role_PANCANCER' == 'CGC_role_PANCANCER', 
#                                      'mutation_multiplicity' == 'mutation_multiplicity', 
#                                      'VEP.IMPACT' == 'VEP.IMPACT')) %>% 
#   dplyr::filter(!is.na(value)) %>% 
#   dplyr::filter(!is.na(proteomics_code)) %>% 
#   dplyr::filter(hgnc_symbol %in% all_genes) %>% 
#   dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
#   dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
#   dplyr::rename(protein_expression = mean_intensity) %>% 
#   dplyr::rename(rna_expression = value) %>% 
#   dplyr::filter(!is.na(protein_expression)) %>% 
#   dplyr::filter(!is.na(rna_expression))

# rna_vs_prot %>% 
#   filter(sample == '62_VIII', hgnc_symbol == 'TTN') %>% 
#   count(mut_consequence, mutation_multiplicity, karyotype)
# 
# rna_vs_prot %>% 
#   group_by(sample, hgnc_symbol, tot_cna) %>%
#   summarise(n = n()) %>% 
#   # filter(n  > 1 ) %>% 
#   dim
 

# saveRDS(rna_vs_prot, 'data/rna_vs_prot_data.rds')

# new model fitted!
# set.seed(176613)
# 
# fit_multiresponse_model = glm_fit_multiresponse(rna_vs_prot,
#                                                 nobs = 4,
#                                                 response = c('protein_expression', 'rna_expression'),
#                                                 model = as.formula('~ n_low + n_alt + n_trunc + n_wt'),
#                                                 alphas = 1,
#                                                 lambda = 'lambda.min')
# 
# saveRDS(fit_multiresponse_model, 'data/multiresponse_glm_lasso_multiplicity_v1.rds')

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


