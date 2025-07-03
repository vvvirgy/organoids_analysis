library(tidyverse)
library(CNAqc)

samples = readRDS('~/Dropbox/Organoids_Accelerator/data/samples_sheet.rds')

cnaqc_obj = lapply(samples$PDO %>% unique, function(x) {
  readRDS(samples$CNAqc[x])
})

names(cnaqc_obj) = lapply(cnaqc_obj, function(x) {
  x$sample
}) %>% unlist

genes = readRDS('~/Google Drive/Il mio Drive/PhD/organoids/cancer_gene_census_positions.rds')

cnaqc_obj_v2 = lapply(cnaqc_obj, function(x) {
  x$mutations = x$mutations %>% 
    dplyr::mutate(is_driver = case_when(VEP.SYMBOL %in% genes$hgnc_symbol ~ TRUE, .default = is_driver)) %>% 
    dplyr::mutate(driver_label = case_when(VEP.SYMBOL %in% genes$hgnc_symbol ~ VEP.SYMBOL, .default = NA))
  return(x)
})

genes_phasing = lapply(cnaqc_obj, function(x) {
  x$phasing %>% 
    dplyr::filter(VEP.SYMBOL %in% genes$hgnc_symbol)
})

saveRDS(genes_phasing, 'genes_phasing.rds')

# define internal function to check which are the real segments of the gene

checking = function(pos, maxFrom, minTo, seg) {
  driver_seg = seg %>%
    filter(from >= maxFrom$from & to <= minTo$to) %>% 
    mutate(driver = as.character(pos["hgnc_symbol"]))
  return(driver_seg)
}

# get the cn state for each desidered gene across the genome

get_driver_cnstate = function(x, drivers) {
  segments = CNAqc::CNA(x)
  
  res = apply(drivers, 1, function(y) {
    seg = segments %>% 
      filter(chr == as.character(y["chr"])) 
    
    seg_wt_driver_from = seg %>% 
      filter(from <= as.numeric(y["from"])) %>% 
      slice_max(from)
    seg_wt_driver_to = seg %>% 
      filter(to >= as.numeric(y["to"])) %>%
      slice_min(to)       
    
    result = checking(pos = y, maxFrom = seg_wt_driver_from, minTo = seg_wt_driver_to, seg)
    return(result)
  }) %>% bind_rows()
  
  return(res)
}

driver_cnstate_cohort = lapply(cnaqc_obj, function(x) {
  x$cna = CNAqc:::relative_to_absolute_coordinates(x = x, cna = x$cna)
  get_driver_cnstate(x, genes)
})

driver_cnstate_cohort = lapply(names(driver_cnstate_cohort), function(x) {
  driver_cnstate_cohort[[x]] %>% 
    mutate(sample = x)
}) %>% bind_rows()

saveRDS(driver_cnstate_cohort, '~/Google Drive/Il mio Drive/PhD/organoids/cgs_pos_cna_state.rds')
