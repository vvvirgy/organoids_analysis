classify_mutations = function(data) {
  data = data %>% 
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
  
  data = data %>% 
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
  return(data)
}