rm(list=ls())
library(tidyverse)
library(ggupset)

META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

MIN_SAMPLES = 3

sf_method = 'psinorm'
use_stable = FALSE
RES_PATH = paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results/sf_", sf_method, "_stable_", use_stable)

df = readRDS(RES_PATH, "lfc_prot_and_rna_bind.rds")
df = df %>% 
  dplyr::group_by(karyotype, name) %>% 
  dplyr::filter(n() == 2)

ggenes = df$name %>% unique
# Get gene/karyotypes with at least two samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>% 
  dplyr::group_by(hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  # dplyr::summarise(n = n()) %>% 
  dplyr::filter(n() >= MIN_SAMPLES) %>% 
  dplyr::rename(name = hgnc_symbol)

karyotypes_df_good %>% 
  filter(name %in% ggenes) %>% 
  dplyr::select(name, karyotype) %>%
  distinct() %>% 
  mutate(tt = 1) %>% 
  pivot_wider(names_from = karyotype, values_from = tt, values_fill = 0) %>% 
  tibble::column_to_rownames('name') %>% 
  UpSetR::upset(., order.by = 'freq', nintersects = 25, main.bar.color = '#088395', matrix.color = '#088395', sets.bar.color = '#088395') 


png('res/upsets/karyotypes_dist_200_genes.png', width = 45, height = 20, units = 'cm', res = 300)
karyotypes_df_good %>% 
  dplyr::select(name, karyotype) %>% 
  group_by(name) %>% 
  summarise(karyotype = list(karyotype), .groups = "drop") %>%  
  sample_n(size = 200) %>%
  ggplot(aes(karyotype)) + 
  geom_bar(fill = '#088395', width = .8) + 
  scale_x_upset() + 
  theme_light() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "#088395",
    combmatrix.label.extra_spacing = 12,
    combmatrix.panel.line.color = '#088395',
    combmatrix.label.make_space = TRUE
  )
dev.off()


karyotypes_df_good %>% 
  ggplot(aes(
    karyotype
  )) + 
  geom_bar(stat = 'count') + 
  theme_bw()


karyotypes_df_good %>% 
  dplyr::select(name, karyotype) %>% 
  group_by(name) %>% 
  distinct() %>% 
  mutate(karyotype = factor(karyotype)) %>% 
  arrange(karyotype) %>% 
  summarise(karyo_dist = paste(karyotype, collapse = ', ')) %>% 
  group_by(name) %>% 
  # summarise(karyo_dist = list(karyo_dist), .groups = "drop") %>%  
  # sample_n(size = 200) %>% 
  ggplot(aes(karyo_dist)) + 
  geom_bar(fill = '#088395', width = .8) + 
  coord_flip()
  

scale_x_upset() + 
  theme_light() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "#088395",
    combmatrix.label.extra_spacing = 12,
    combmatrix.panel.line.color = '#088395',
    combmatrix.label.make_space = TRUE
  )
  
