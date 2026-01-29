rm(list=ls())

# .libPaths()
library(tidyverse)
library(parallel)
library(DEP)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
# source('organoids_analysis/R/functions_utils/fit_plots.R')
# source('organoids_analysis/R/scRNA/cna_comparison_utils.R')
# source('organoids_analysis/R/proteomics_de_utils.R')

dep_res = readRDS('data/dep_by_karyotype_results_v2.rds')

# getting the results of dge for the specific gene 
get_dge_res = function(dep, gene) {
  
  dep %>% 
    get_df_long() %>% 
    filter(PG.Genes == gene) 
  
}

dep_res_fc = parallel::mclapply(dep_res, function(x) {
  get_dge_res(dep = x$dep, gene = x$gene)
}, mc.cores = 2)

print('all fc extracted')

dep_res_fc = bind_rows(dep_res_fc)
saveRDS(dep_res_fc, 'data/dep_res_fc_v2.rds')

dep_res_fc_red = dep_res_fc %>% 
  select(PG.Genes, condition, ends_with('CI.L'), ends_with('CI.R'), ends_with('diff'), 
         ends_with('p.adj'), ends_with('p.val'), ends_with('significant'))
dep_res_fc_red = dep_res_fc_red %>% 
  distinct()

print('starting to merge all results')
get_fc_tb_clean = function(tb) {
  
  lapply(tb$condition %>% unique, function(x) {
    df = tb %>% 
      select(PG.Genes, starts_with(x)) %>%
      distinct()
    colnames(df) = gsub(paste0(x, '_vs_X1.1_'), '', colnames(df))
    
    df = df %>% 
      mutate(condition = x) %>% 
      relocate(condition, .after = PG.Genes) %>% 
      filter(condition != 'X1.1')
    
    return(df)
    
   }) %>% 
    bind_rows()
  
}

fc_tb_clean = get_fc_tb_clean(dep_res_fc_red)
print('results merged')
saveRDS(fc_tb_clean, 'data/fc_tb_clean_v2.rds')
  