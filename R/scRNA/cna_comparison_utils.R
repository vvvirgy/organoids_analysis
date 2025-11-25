# utilities when comparing expression and copy numbers

filter_fragmented_cnas = function(x, 
                                  samples_list, 
                                  genes_to_check,
                                  genes_position = genes_pos,
                                  min_length = 10^6, 
                                  strategy = c('MMSl', # max segment length
                                               'FS', # first segment that covers the gene 
                                               'MOv' # maximise the overlap
                                               )) {
  # treat multihits
  x = x %>% 
    dplyr::filter(sample %in% samples_list) %>%
    dplyr::filter(hgnc_symbol %in% genes_to_check) %>%
    tidyr::separate(segment_id, into = c('chr', 'segment_from', 'segment_to'), sep = ':', convert = TRUE) %>% 
    dplyr::group_by(sample, hgnc_symbol, is_mutated) %>%
    group_by(sample, hgnc_symbol) %>% 
    mutate(
      has_multihit = n() > 1 & any(is_mutated == TRUE),
      mut_consequence = case_when(
        has_multihit ~ "multihit",
        TRUE ~ mut_consequence
      )
    ) %>% 
    dplyr::mutate(length = segment_to - segment_from) %>% 
    filter(mut_consequence != 'multihit')
    # dplyr::mutate(mut_consequence = case_when(
    #   (n()>1 & is_mutated == TRUE) ~ 'multihit', 
    #   .default = mut_consequence
    # )) %>% 
    # group_by(sample, hgnc_symbol) %>% 
    # mutate(mut_consequence = ifelse(any(mut_consequence) == 'multihit', 'multihit', mut_consequence))
    
    # 
    # dplyr::filter(n() == 1)
  
  # handle multisegments situations
  if (strategy == 'MMSl') {
    x = x %>% 
      group_by(sample, hgnc_symbol) %>% 
      slice_max(length)
  } 
  
  if(strategy == 'FS') {
    x = x %>% 
      group_by(sample, hgnc_symbol) %>% 
      slice_max(segment_from)
  } 
  
  if (strategy == 'MOv') {
    
    # find genes with different karyotypes 
    
    # removing genes with multiple hits on the same segments (can't handle them in this model!)
    multihit_genes = x %>% 
      distinct(sample, hgnc_symbol, QC_PASS, mut_consequence, segment_from, segment_to, chr, karyotype) %>% 
      dplyr::group_by(sample, hgnc_symbol, segment_from, segment_to, karyotype) %>%
      dplyr::count() %>% # detect genes that appear more than once
      dplyr::filter(n > 1) %>% 
      dplyr::select(-n)
    
    # filter out the combination of genes, sample and segments with multihit 
    x = x %>% 
      dplyr::anti_join(., multihit_genes, by = join_by(
        'sample' == 'sample', 
        'hgnc_symbol' == 'hgnc_symbol', 
        'segment_from' == 'segment_from',  
        'segment_to' == 'segment_to', 
        'karyotype' == 'karyotype'
      )) %>% 
      dplyr::filter(length > min_length) %>% 
      dplyr::filter(QC_PASS == TRUE)
    
    multiple_segments = x %>% 
      dplyr::select(-c(from_gene, to_gene)) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(sample, hgnc_symbol, segment_from, segment_to, karyotype) %>%
      dplyr::count() %>% # detect genes that appear more than once
      dplyr::filter(n == 1) %>% # remove multihit genes (same segment)
      dplyr::group_by(sample, hgnc_symbol, karyotype) %>%
      dplyr::count() %>%  # count how many genes have a different cna state among segments
      dplyr::filter(n == 1) %>%  # select only genes that have multiple segments associated
      dplyr::group_by(sample, hgnc_symbol) %>%
      dplyr::count() %>%
      dplyr::filter(n > 1) %>%
      dplyr::group_by(sample) %>% 
      dplyr::group_split()
    names(multiple_segments) = lapply(multiple_segments, function(x) {x$sample %>% unique}) %>% unlist
    
    multiple_segments = lapply(multiple_segments, function(s) {
      pos = genes_position %>% 
        dplyr::filter(hgnc_symbol %in% unique(s$hgnc_symbol)) %>% 
        dplyr::full_join(., s, by = join_by('hgnc_symbol' == 'hgnc_symbol'))
      return(pos)
    })
    
    x = lapply(names(multiple_segments), function(tt) {
      
      df = x %>% 
        dplyr::filter(sample == tt) %>% 
        dplyr::filter(hgnc_symbol %in% unique(multiple_segments[[tt]]$hgnc_symbol))
      
      df = df %>% 
        dplyr::full_join(., multiple_segments[[tt]], 
                         by = join_by('hgnc_symbol' == 'hgnc_symbol', 
                                      'chr' == 'chr', 
                                      'sample' == 'sample'))
      
      df = df %>% 
        dplyr::mutate(OV_gene_segment = case_when(
          segment_to > to ~ (to - segment_from),
          to > segment_to ~ (segment_to - from), 
          .default = NA
        )) %>% 
        dplyr::filter(OV_gene_segment > 0) %>% 
        dplyr::group_by(hgnc_symbol) %>% 
        dplyr::slice_max(OV_gene_segment)
      
      x_v2 = x %>% 
        dplyr::filter(sample == tt) %>% 
        dplyr::filter(!hgnc_symbol %in% unique(multiple_segments[[tt]]$hgnc_symbol))
      
      df = df %>% 
        dplyr::select(all_of(colnames(x_v2)))
      
      bind_rows(x_v2, df)
      
    }) %>% 
      bind_rows()
    
  }
  
  return(x)
  
}



extract_lm_per_gene = function(genes_cna, formula) {
  
  genes_cna = genes_cna %>% 
    group_by(hgnc_symbol) %>% 
    group_split()
  
  names(genes_cna) = lapply(genes_cna, function(x) {
    x$hgnc_symbol %>% unique
  }) %>% unlist
  
  tt = lapply(genes_cna %>% names, function(x) {
    mm = lm(as.formula(formula), genes_cna[[x]])
    tt = broom::tidy(mm) %>% 
      dplyr::mutate(gene = x) %>% 
      dplyr::mutate(r.squared = summary(mm)$r.squared)
    return(tt)
  }) %>% 
    dplyr::bind_rows()
  
  return(tt)
  
}

get_good_candidates_rsquared = function(x, n) {
  x %>% 
    dplyr::arrange(desc(r.squared)) %>% 
    head(n)
}

get_good_candidates_estimate = function(x, n, which) {
  if(which == 'up') {
    x %>% 
      dplyr::arrange(desc(estimate)) %>% 
      head(n)
    } else {
        x %>% 
          dplyr::arrange(desc(estimate)) %>% 
          tail(n)
      }
}
