setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
rm(list = ls())
.libPaths()
library(tidyverse)
library(dplyr)
# library(boot)
# library(clusterProfiler)
# library(ReactomePA)
# library(org.Hs.eg.db)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/dge_utils_and_plots.R')

min_samples = 1

# dividing the genes in different groups ---- 
## dividing by karyotypes -----

df = readRDS(paste0('data/compensation/cs_results_min_', min_samples,'.rds'))
karyos = c('1:0', '2:0', '2:1', '2:2')

df %>% 
  ggplot(
    aes(
      x = karyotype, 
      y = CS, 
      fill = omic
    )
  ) + 
  geom_violin() + 
  theme_bw() +
  ggpubr::stat_compare_means(comparisons = list(c('1:0', '2:1')), ) + 
  facet_wrap(vars(omic))

df = df %>% 
  mutate(omic = gsub('protein', 'Protein', omic)) %>% 
  mutate(omic = factor(omic, levels = c('RNA', 'Protein')))

# setting the thr karyotype specific 
df_groups = lapply(karyos, function(x) {
  define_th(df, q = .5, karyo = x)  
})

df_groups = df_groups %>% 
  bind_rows()

saveRDS(df_groups, 'data/compensation/cs_classification_by_karyo_min_1.rds')

## karyotypes aggregated -----

df_mean = readRDS(paste0('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/compensation/cs_mean_across_karyos_min_', min_samples,'.rds')) %>% 
  mutate(omic = gsub('protein', 'Protein', omic)) %>% 
  mutate(omic = factor(omic, levels = c('RNA', 'Protein')))
  # dplyr::select(-c(lfc, n, DNA_lfc, CS)) %>% 
  # tidyr::pivot_wider(values_from = mean_CS, names_from = omic)

q = .5

rna_cuts_v2 = df_mean %>%
  filter(omic == 'RNA') %>% 
  dplyr::mutate(s = sign(mean_CS)) %>%
  dplyr::group_by(s, omic) %>%
  dplyr::summarise(q = quantile(mean_CS, q)) %>%
  dplyr::pull(q) %>%
  sort()

prot_cuts_v2 = df_mean %>%
  filter(omic == 'Protein') %>% 
  dplyr::mutate(s = sign(mean_CS)) %>%
  dplyr::group_by(s) %>%
  dplyr::summarise(q = quantile(mean_CS, q)) %>%
  dplyr::pull(q) %>%
  sort()

# assign genes to cs class

df_groups <- df_mean %>%
  dplyr::select(-n) %>% 
  pivot_wider(names_from = omic, values_from = c(lfc, DNA_lfc, CS, mean_CS)) %>% 
  mutate(reg_group = case_when(
    # Group 1: High RNA CS (>65th) and Low Protein CS (<35th)
    CS_RNA > rna_cuts[2] & CS_Protein < prot_cuts[1] ~ "RNA-comp.",
    CS_RNA > rna_cuts[2] & CS_Protein > prot_cuts[2] ~ "Fully compensated",
    
    # Group 2: Low RNA CS (<35th) and High Protein CS (>65th)
    CS_RNA < rna_cuts[1] & CS_Protein > prot_cuts[2] ~ "Prot-comp.",
    CS_RNA < rna_cuts[1] & CS_Protein < prot_cuts[1] ~ "Hyper-responders",
    
    TRUE ~ "Intermediate/Other"
  ))
saveRDS(df_groups, paste0('data/compensation/cs_classification_min_', min_samples, '.rds'))
