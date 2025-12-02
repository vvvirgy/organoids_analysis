.libPaths()
rm(list=ls())
library(tidyverse)
library(CNAqc)
# source('~/Desktop/lade_scratch/CDSlab/brca_cro_prj/BRCA_CRO/new_analysis/functions/CNA/get_cna_status_region_specific.R')
source('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/BRCA_CRO/new_analysis/functions/CNA/get_cna_status_region_specific.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/preparation/get_cna_function.R')

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')
cnas = readRDS('data/cnaqc/cnas_list.rds')
# cnas = lapply(cnas, function(x) {CNAqc::wt_mutant_alleles(x)})
# cnas = lapply(cnas, function(x) {CNAqc::compute_CCF(x, method = 'ROUGH')})
# saveRDS(cnas, 'data/cnaqc/cnas_list_rough_ccf.rds')
# 
# cnas = readRDS('data/cnaqc/cnas_list_rough_ccf.rds')

# CNAqc::plot_multisample_CNA(cnas)
# ggsave('res/all_organoids_cna_cohort.png')

# coad_genes = readRDS('~/Google Drive/Il mio Drive/PhD/organoids/cancer_gene_census_positions.rds')
coad_genes = readRDS('data/all_genes_positions_info.rds')
coad_genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol))

genes_karyo = lapply(names(cnas), function(x) {
  
  lapply(names(cnas[[x]]$per_chromosome), function(chr) {
    
    dat = cnas[[x]]$per_chromosome[[chr]]
    genes = coad_genes %>% 
      dplyr::filter(chr ==!! chr)
    
    get_cna_status(dat, drivers = genes) %>% 
      dplyr::mutate(sample = x)
    
  }) %>% 
    dplyr::bind_rows()
  # get_cna_status(cnas[[x]], drivers = coad_genes) %>% s
  #   dplyr::mutate(sample = x)
}) %>% 
  dplyr::bind_rows()

# saveRDS(genes_karyo, '~/Desktop/cdslab_scratch/organoids_prj/data/karyotypes_all_genes.rds')
saveRDS(genes_karyo, '/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_all_genes_qc_v3.rds')
# usa solo informazione del karyotype

# include mutation types

# phased_drivers = lapply(names(cnas), function(x) {
#   # add eventually a column sample + ogni gene in cgs per coad è driver
#   print(x)
#   cnas[[x]]$mutations = cnas[[x]]$mutations %>% 
#     dplyr::mutate(sample = x) %>% 
#     # dplyr::filter(VEP.SYMBOL %in% (coad_genes$hgnc_symbol %>% unique)) %>% 
#     dplyr::mutate(is_driver = ifelse(VEP.SYMBOL %in% (coad_genes$hgnc_symbol %>% unique), TRUE, FALSE)) %>% 
#     dplyr::mutate(driver_label = ifelse(is_driver == TRUE, VEP.SYMBOL, ''))
#   
#   CNAqc::wt_mutant_alleles(cnas[[x]], gene_column = 'VEP.SYMBOL')
# })

# muts = lapply(names(cnas), function(x) {
#     # add eventually a column sample + ogni gene in cgs per coad è driver
#     print(x)
#     cnas[[x]]$mutations %>%
#       dplyr::mutate(sample = x) %>% 
#       dplyr::filter(VEP.SYMBOL %in% (coad_genes$hgnc_symbol %>% unique)) %>% 
#       dplyr::select(sample, chr, VAF, VEP.SYMBOL, karyotype)
#     # CNAqc::wt_mutant_alleles(cnas[[x]], gene_column = 'VEP.SYMBOL')
# }) %>% 
#   dplyr::bind_rows()

# genes_karyo_v2 = right_join(genes_karyo, muts, 
#                            by = join_by('sample' == 'sample', 
#                                         'karyotype' == 'karyotype',
#                                         'region' == 'VEP.SYMBOL', 
#                                         'chr' == 'chr'
#                            ))
# 
# genes_karyo_v2 = genes_karyo_v2 %>% 
#   dplyr::group_by(region, sample, chr, karyotype) %>% 
#   dplyr::distinct(.keep_all = FALSE)
