library(tidyverse)
library(biomaRt)

get_grch38_genomics_positions = function(genes, cnaqc_obj_path) {
  
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")   
  attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position")
  gene_coords <- getBM(attributes = attributes,
                       filters = "hgnc_symbol",
                       values = genes,
                       mart = ensembl)
  
  genes_pos = gene_coords %>% 
    dplyr::arrange(chromosome_name) %>%
    dplyr::rename(chr = chromosome_name, from = start_position, to = end_position) %>% 
    dplyr::mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::filter(chr %in% paste0('chr', c(seq(1:22), 'X', 'Y')))
  
  cnaqc = readRDS(cnaqc_obj_path)
  
  genes_pos = CNAqc:::relative_to_absolute_coordinates(cna = genes_pos, x = cnaqc)  
  
  return(genes_pos)
}