source('/orfeo/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/preparation/get_cnas_muts.R')
setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')

cnas = readRDS('data/cnaqc/cnas_list.rds')
genes_karyo = readRDS('data/karyotypes_all_genes_qc.rds')

coad_genes = readRDS('data/all_genes_positions_info.rds')
coad_genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::pull(hgnc_symbol)

genes_karyo_muts = lapply(cnas, function(x) {
  extract_mutational_status(x, karyotypes = genes_karyo, genes_pos = coad_genes)
}) %>% 
  bind_rows()

saveRDS(genes_karyo_muts, 'data/karyotypes_mutations_all_genes_qc.rds')

# 
# genes_karyo_muts %>% 
#   filter(hgnc_symbol == 'KRAS') %>% 
#   filter(QC_PASS == TRUE)
