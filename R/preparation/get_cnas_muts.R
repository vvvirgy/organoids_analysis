library(tidyverse)
library(CNAqc)

extract_mutational_status = function(x, #cnaqc obj
                                     karyotypes, #result of the previous iteration
                                     genes_pos # dataframe or list of genes
) {
  
  muts = relative_to_abs_coords(Mutations(x), x$reference_genome)
  cna_sample = karyotypes %>% 
    dplyr::filter(sample == x$sample) %>% 
    tidyr::separate(segment_id, sep = ':', into = c('.chr', 'segment_from', 'segment_to'), convert = TRUE) %>% 
    dplyr::mutate(.chr = NULL) %>% 
    tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':',  convert = TRUE) %>% 
    dplyr::rename(from_gene = from) %>%
    dplyr::rename(to_gene = to)
    # dplyr::mutate(from = NULL, to = NULL)
  
  muts = muts %>%  
    dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, VEP.SYMBOL, VEP.Consequence, driver_label, is_driver, segment_id, QC_PASS, karyotype) %>% 
    tidyr::separate(karyotype, into = c('Major', 'minor'), sep = ':',  convert = TRUE) %>% 
    dplyr::filter(VEP.SYMBOL != '.') %>% 
    tidyr::separate(segment_id, into = c('.chr', 'segment_from', 'segment_to'), sep = ':', convert = TRUE) %>% 
    dplyr::mutate(.chr = NULL) %>% 
    tidyr::separate(VEP.Consequence, into = 'mut_consequence', sep = '&') %>% 
    dplyr::rename(mut_from = from) %>% 
    dplyr::rename(mut_to = to) %>% 
    dplyr::filter(VEP.SYMBOL %in% unique(genes_pos)) %>% 
    dplyr::rename(from = segment_from) %>% 
    dplyr::rename(to = segment_to) %>% 
    relative_to_abs_coords(., x$reference_genome) %>% 
    dplyr::rename(segment_from = from) %>% 
    dplyr::rename(segment_to = to) %>% 
    dplyr::mutate(sample = x$sample)
  
  muts_with_cna = cna_sample %>% 
    dplyr::full_join(., muts, by = join_by(
      'chr' == 'chr', 
      'hgnc_symbol' == 'VEP.SYMBOL', 
      'Major' == 'Major',
      'minor' == 'minor',
      'segment_from' == 'segment_from', 
      'segment_to' == 'segment_to', 
      'QC_PASS' == 'QC_PASS', 
      'sample' == 'sample'
    ))
  
  muts_with_cna = muts_with_cna %>% 
    dplyr::mutate(segment_id = paste(chr, segment_from, segment_to, sep = ':'), 
                  karyotype = paste(Major, minor, sep = ':'), 
                  is_mutated = ifelse(is.na(mut_consequence), FALSE, TRUE), 
                  Major = NULL, 
                  minor = NULL, 
                  segment_from = NULL, 
                  segment_to = NULL)
    
  return(muts_with_cna)

}


# mapping mutations on segments
relative_to_abs_coords = function(x, ref) {
  # get reference genome coordinates
  reference_genome = CNAqc:::get_reference(ref)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  x %>% 
    mutate(from = from + vfrom[chr], to = to + vfrom[chr])
}


  
