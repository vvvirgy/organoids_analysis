rm(list=ls())

# .libPaths()
library(tidyverse)
library(parallel)
library(DEP)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
source('organoids_analysis/R/functions_utils/constants.R')


print('reading data')
dep_res = readRDS(file.path(data_path, 'noise_model/protein/dep_tetraploid_genes_results.rds'))

# getting the results of dge for the specific gene 
get_dge_res = function(dep, gene) {
  
  dep %>% 
    get_df_long() %>% 
    filter(PG.Genes == gene) 
  
}

print('data read')

dep_res_fc = parallel::mclapply(dep_res, function(x) {
  
  print(x$gene)
  
  tryCatch(
    {lapply(x$iterations, function(s) {
      get_dge_res(dep = s$dep, gene = x$gene) %>% 
        mutate(iteration = s$iteration)
    }) %>% 
        bind_rows()}, 
    error = function(e) {
      NULL
    })
  
}, mc.cores = 2)
print('all fc extracted')

# bind and filter out null element
dep_res_fc <- Filter(Negate(is.null), dep_res_fc) 
dep_res_fc = bind_rows(dep_res_fc)
saveRDS(dep_res_fc, file.path(data_path, 'noise_model/protein/dep_res_fc.rds'))

# merge everything together
dep_res_fc_red = dep_res_fc %>%
  dplyr::select(PG.Genes, condition, iteration, ends_with('CI.L'), ends_with('CI.R'), ends_with('diff'),
         ends_with('p.adj'), ends_with('p.val'), ends_with('significant'))
dep_res_fc_red = dep_res_fc_red %>%
  distinct()

print('starting to merge all results')

fc_tb_clean = dep_res_fc_red %>%
  dplyr::filter(condition == 'B') %>%
  group_by(PG.Genes) %>%
  mutate(mean_fc = mean(B_vs_A_diff),
         sd_fc = sd(B_vs_A_diff))

print('results merged')
saveRDS(fc_tb_clean, file.path(data_path, 'noise_model/protein/fc_tb_clean.rds'))
