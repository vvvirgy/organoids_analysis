.libPaths()
rm(list=ls())
library(tidyverse)

fc_tb_clean = readRDS('data/fc_tb_clean.rds')
head(fc_tb_clean)

fc_tb_clean = fc_tb_clean %>% 
  mutate(condition = gsub('X', '', condition)) %>% 
  mutate(condition = gsub('\\.', ':', condition)) %>% 
  rename(genes = PG.Genes) %>% 
  as_tibble()

fc_tb_clean %>% 
  mutate(sign = ifelse(p.adj <= 0.05, TRUE, FALSE)) %>% 
  filter(!is.na(diff)) %>% 
  ggplot(aes(
    x = diff, 
    y = -log(p.adj, base = 10), 
    color = sign
  )) + 
  geom_point() + 
  facet_wrap(~condition, nrow = 1) + 
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(y = '-log10(padj)',
       x = 'log2(FC)')
ggsave('res/proteomics_de_fc.png', width = 10, height = 6)

all_lfc_res <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results/all_lfc_res.rds")

all_lfc_res %>% 
  mutate(sign = ifelse(adj_pval <= 0.05, TRUE, FALSE)) %>% 
  filter(!is.na(lfc)) %>% 
  ggplot(aes(
    x = lfc, 
    y = -log(adj_pval, base = 10), 
    color = sign
  )) + 
  geom_point() + 
  facet_wrap(~coef, nrow = 1) + 
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(y = '-log10(padj)',
       x = 'log2(FC)')
  

corr_omics_fc = all_lfc_res %>% 
  mutate(condition = gsub('karyotype', '', coef)) %>% 
  full_join(., fc_tb_clean, by = join_by('name' == 'genes', 'condition' == 'condition'))

corr_omics_fc %>% 
  mutate(diff = ifelse(is.na(diff), 0, diff)) %>% 
  mutate(lfc = ifelse(is.na(lfc), 0, lfc)) %>% 
  mutate(significance = case_when(
    (adj_pval <= .05 & p.adj <= .05) ~ 'both', 
    (adj_pval <= .05 & p.adj > .05 | adj_pval <= .05 & is.na(p.adj)) ~ 'RNA', 
    (adj_pval > .05 & p.adj <= .05 | is.na(adj_pval) & p.adj <= .05) ~ 'Protein', 
    .default = 'None'
  )) %>% 
  ggplot(aes(lfc, diff, color = significance)) + 
  geom_point() + 
  facet_wrap(significance~condition) + 
  theme_bw()



