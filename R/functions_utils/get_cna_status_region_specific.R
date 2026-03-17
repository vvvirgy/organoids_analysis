# define internal function to check which are the real segments of the gene
checking = function(pos, maxFrom, minTo, seg, label) {
  driver_seg = seg %>%
    dplyr::filter(from >= maxFrom$from & to <= minTo$to) %>% 
    dplyr::mutate(region = label)
  return(driver_seg)
}

# get the cn state for each desidered gene across the genome

get_cna_status_reduced = function(x, drivers) {
  
  segments = CNAqc:::relative_to_absolute_coordinates(x, CNA(x))
  
  # segments = CNAqc::CNA(x)
  
  
  print(x$sample)
  
  res = apply(drivers, 1, function(y) {
    
    # y[3]
    seg = segments %>%
      dplyr::filter(chr == as.character(y[1]))
    # 
    seg_wt_driver_from = seg %>%
      filter(from <= as.numeric(y[2])) %>%
      slice_max(from)
    
    if(nrow(seg_wt_driver_from)==0) { 
      seg_wt_driver_from = seg %>% 
        slice_min(from)
      }
    
    # print(seg_wt_driver_from)
    # print(y[2])
    seg_wt_driver_to = seg %>%
      filter(to >= as.numeric(y[3])) %>%
      slice_min(to)
    
    if(nrow(seg_wt_driver_to) == 0) {
      seg_wt_driver_to = seg %>% 
        slice_max(to)
    }
    # 
    # seg_wt_driver_from
    # seg_wt_driver_to
    
    print(as.character(y[4]))
    result = checking(pos = y, 
                      maxFrom = seg_wt_driver_from, 
                      minTo = seg_wt_driver_to, 
                      seg, 
                      label = as.character(y[4]))
    
    result = result %>% 
      dplyr::mutate(karyotype = paste(Major, minor, sep = ':')) %>% 
      # dplyr::mutate(segment_id = paste(chr, from, to, sep = ':')) %>% 
      dplyr::rename(segment_from = from) %>% 
      dplyr::rename(segment_to = to) %>% 
      dplyr::rename(hgnc_symbol = region) %>% 
      dplyr::mutate(from = y[2], to = y[3]) %>% 
      # dplyr::relocate(from, .after = chr) %>% 
      # dplyr::relocate(to, .after = from) %>% 
      # dplyr::relocate(karyotype, .after = hgnc_symbol) %>% 
      # dplyr::relocate(segment_id, .after = karyotype) %>% 
      dplyr::select(c(chr, from, to, hgnc_symbol, karyotype, segment_from, segment_to, QC_PASS))
    
    return(result)
  }) %>% bind_rows()
  
  return(res)
}
