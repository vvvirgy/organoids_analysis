library(CNAqc)
library(tidyverse) 

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')
cnas_path = list.files('data/cnaqc', full.names = T)
cnas = lapply(cnas_path, readRDS)
names(cnas) = lapply(cnas, function(x) {x$sample}) %>% unlist

cnas = lapply(cnas, function(x) {
  CNAqc::analyze_peaks(x)
})

saveRDS(cnas, 'data/cnaqc/cnas_list.rds')

cnas = lapply(cnas, function(x) {
  x$per_chromosome = x %>% 
    CNAqc::split_by_chromosome() %>% 
    lapply(., function(y) {
      CNAqc::analyze_peaks(y)
    })
  
  return(x)
})

saveRDS(cnas, 'data/cnaqc/cnas_list.rds')
