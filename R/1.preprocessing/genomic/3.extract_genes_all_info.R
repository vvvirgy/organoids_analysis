source('/orfeo/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/preprocessing/get_cnas_muts.R')
setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')

# cnas = readRDS('data/cnaqc/cnas_list_rough_ccf.rds')
cnas = readRDS('data/cnaqc_v2/cnas_list_v2.rds')
genes_karyo = readRDS('data/processed_data/karyotypes_all_genes_qc.rds')

coad_genes = readRDS('data/all_genes_positions_info.rds')
coad_genes = coad_genes %>% 
  dplyr::relocate(hgnc_symbol, .after = to) %>% 
  # dplyr::filter(chr %in% c('chr5', 'chr17')) %>% 
  dplyr::filter(!grepl('LINC', hgnc_symbol)) %>% 
  dplyr::pull(hgnc_symbol)

# x = cnas[[1]]

genes_karyo_muts = lapply(cnas, function(x) {
  print(x$sample)
  tryCatch({extract_mutational_status(x, 
                                      karyotypes = genes_karyo, 
                                      genes_pos = coad_genes, 
                                      which = 'CCF')}, 
           error = function(e) {
             NA
           })
}) 
genes_karyo_muts = Filter(function(x) {!is.null(nrow(x))}, genes_karyo_muts)
genes_karyo_muts = genes_karyo_muts %>% 
  bind_rows()

saveRDS(genes_karyo_muts, 'data/processed_data/karyotypes_mutations_all_genes_qc_ccf.rds')


# extract using the mutations and not caring of multiplicity
genes_karyo_muts = lapply(cnas, function(x) {
  print(x$sample)
  tryCatch({extract_mutational_status(x, 
                                      karyotypes = genes_karyo, 
                                      genes_pos = coad_genes, 
                                      which = 'muts')}, 
           error = function(e) {
             NA
           })
}) 
genes_karyo_muts = Filter(function(x) {!is.null(nrow(x))}, genes_karyo_muts)
genes_karyo_muts = genes_karyo_muts %>% 
  bind_rows()

saveRDS(genes_karyo_muts, 'data/karyotypes_mutations_all_genes_qc_only_muts.rds')

# 
# genes_karyo_muts %>% 
#   filter(hgnc_symbol == 'KRAS') %>% 
#   filter(QC_PASS == TRUE)
