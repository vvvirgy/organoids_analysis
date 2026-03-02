rm(list = ls())
.libPaths()
library(tidyverse)
library(dplyr)
library(boot)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(googlesheets4)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/dge_utils_and_plots.R')

# computing the CS -----
DATA_PATH = '/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results'
SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_karyo_all_organoids_filt.rds"
META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

MIN_SAMPLES = 2

# Get gene/karyotypes with at least two samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>% 
  dplyr::group_by(hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= MIN_SAMPLES) %>% 
  dplyr::rename(name = hgnc_symbol)

df_dna = readRDS(paste(DATA_PATH, "DNA_lfc.rds", sep = '/')) %>% dplyr::rename(DNA_lfc = lfc)

df = readRDS(paste(DATA_PATH, "sf_psinorm_stable_FALSE/lfc_prot_and_rna_bind.rds", sep = '/'))  

df = df %>% 
  dplyr::group_by(karyotype, name) %>% 
  dplyr::filter(n() == 2)

df = df %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

df_dna = df_dna %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

df %>% 
  dplyr::select(name, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise("n_distinct_karyotypes" = n()) %>% 
  dplyr::count(n_distinct_karyotypes) %>% 
  dplyr::rename(n_genes = n) 

df %>% 
  dplyr::select(name, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() %>% 
  dplyr::select(karyotype) %>% 
  dplyr::count(karyotype) 

df = df %>% 
  dplyr::filter(omic != "DNA")

df = df %>% 
  ungroup() %>% 
  # dplyr::left_join(df_dna %>% dplyr::select(!omic), by = join_by("name" == "name", "karyotype" == "karyotype")) %>% 
  dplyr::left_join(df_dna) 

# computing CS per karyotype
df = df %>%
  # mutate(CS = DNA_lfc - lfc) #%>%
  dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc)) 

# weighting the CS per number of samples for each karyotype
df_tot = df %>% 
  group_by(name, omic) %>% 
  mutate(
    CS_tot = sum(CS*n)/sum(n)
  )
  

df_mean_CS = df %>%
  dplyr::group_by(name, omic) %>%
  dplyr::mutate(mean_CS = mean(CS)) #%>%
  # tidyr::pivot_wider(values_from = mean_CS, names_from = omic)

saveRDS(df_mean_CS, 'data/compensation/cs_mean_across_karyos.rds')


p1 = df_tot %>% 
  ggplot(aes(
    y = CS_tot, 
    x = omic
  )) + 
  geom_violin()

p2 = df_mean_CS %>% 
  ggplot(aes(
    y = mean_CS, 
    x = omic
  )) + 
  geom_violin()


df %>% 
  dplyr::select(karyotype, name, CS, omic, DNA_lfc) %>% 
  tidyr::pivot_wider(values_from = CS, names_from = omic) %>% 
  ggplot(mapping = aes(x = RNA, y = protein)) +
  geom_point() + 
  facet_wrap(~karyotype) + 
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red')

saveRDS(df, 'data/compensation/cs_results.rds')

# trying to compute a cs per karyotype

