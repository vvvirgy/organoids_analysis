# get the cn state for each desidered gene across the genome

get_cna_status = function(x, drivers, filter = TRUE) {
  
  segments = CNAqc:::relative_to_absolute_coordinates(x, CNA(x))
  segments = segments %>% 
    dplyr::relocate(QC_PASS, .after = minor)

  # segments = CNAqc::CNA(x)
  
  print(x$sample)
  
  res = apply(segments, 1, function(y) {
     
    # y[3]
    dr = drivers %>%
      dplyr::filter(chr == as.character(y[1])) %>% 
      dplyr::arrange(from) 
    # 
    drivers_on_segment = dr %>%
      dplyr::filter(from >= as.numeric(y[2])) %>%
      dplyr::filter(to <= as.numeric(y[3])) %>% 
      dplyr::mutate(karyotype = paste(y[4], y[5], sep = ':')) %>% 
      # dplyr::mutate(segment_id = paste(y[1], y[2], y[3], sep = ':')) %>% 
      # dplyr::mutate(segment_id = gsub(' ', '', segment_id)) %>% 
      dplyr::mutate(segment_from = as.numeric(y[2]),segment_to = as.numeric(y[3])) %>% 
      dplyr::mutate(QC_PASS = as.logical(y[6]))
    
    return(drivers_on_segment)
      
    # if(nrow(drivers_on_segment)==0) { 
    #   drivers_on_segment = dr %>%
    #     dplyr::filter(from <= as.numeric(y[2])) %>%
    #     dplyr::filter(to <= as.numeric(y[3]))
    # }
    
    # print(seg_wt_driver_from)
    # print(y[2])
    # seg_wt_driver_to = seg %>%
    #   filter(to >= as.numeric(y[3])) %>%
    #   slice_min(to)
    # 
    # if(nrow(seg_wt_driver_to) == 0) {
    #   seg_wt_driver_to = seg %>% 
    #     slice_max(to)
    # }
    # 
    # seg_wt_driver_from
    # seg_wt_driver_to
    
    # print(as.character(y[4]))
    # result = checking(pos = y, 
    #                   maxFrom = seg_wt_driver_from, 
    #                   minTo = seg_wt_driver_to, 
    #                   seg, 
    #                   label = as.character(y[4]))
    # return(result)
  }) %>% bind_rows()
  
  missing_drivers = setdiff(unique(drivers$hgnc_symbol), unique(res$hgnc_symbol))
  
  if(length(missing_drivers) > 0) {
    
    missing_drivers = drivers %>% 
      dplyr::filter(hgnc_symbol %in% missing_drivers)
    
    missing_dr = get_cna_status_reduced(x, missing_drivers) %>% 
      dplyr::mutate(from = as.numeric(from), to = as.numeric(to))
    
    res = bind_rows(res, missing_dr)
  }

  return(res)
}
