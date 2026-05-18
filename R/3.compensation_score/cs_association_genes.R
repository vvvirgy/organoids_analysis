library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
source('organoids_analysis/jovoniR/constants.R')
df_groups = readRDS('data/compensation/cs_classification_gene_labels.rds')
karyos = readRDS('data/karyotypes_genes_filtered_scrna.rds')

df_groups = df_groups %>% 
  full_join(., CNAqc::gene_coordinates_GRCh38, by = join_by('name' == 'gene')) %>% 
  filter(!is.na(reg_group))

df_groups %>%
  dplyr::select(name, chr, reg_group) %>% 
  dplyr::filter(chr %in% paste0('chr', c(seq(1:22), 'X', 'Y'))) %>% 
  mutate(chr = factor(chr, levels = paste0('chr', c(seq(1:22), 'X', 'Y')) %>% rev)) %>% 
  group_by(chr) %>% 
  mutate(n_genes = n()) %>% 
  group_by(chr, reg_group) %>% 
  mutate(prop = n()/n_genes) %>% 
  ggplot(aes(y = chr, 
             x = prop, 
             fill = reg_group)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  theme_bw() + 
  facet_wrap(~reg_group, ncol = 1, scales = 'free_x') + 
  theme(#axis.text.x = element_text(angle = 66, hjust = 0.5, vjust = 0.5),
        legend.position = 'none') + 
  scale_fill_manual(values = category_colors) 
ggsave('res/compensation_score/prop_categories_by_chr.png', width = 6, height = 15, units = 'in', dpi = 300)
ggsave('res/compensation_score/prop_categories_by_chr.pdf', width = 6, height = 15)


hyper_res = df_groups %>% 
  dplyr::filter(reg_group == 'Hyper-responders') %>% 
  pull(name) %>% 
  unique

karyos %>% 
  filter(hgnc_symbol %in% hyper_res) %>% 
  mutate(karyotype = factor(karyotype, levels = rev(c('1:0', '1:1', '2:0', '2:1', '2:2')))) %>% 
  group_by(hgnc_symbol) %>% 
  filter(hgnc_symbol %in% sample(hyper_res, size = 30)) %>% 
  ggplot(aes(hgnc_symbol, fill = karyotype)) + 
  geom_bar(position = 'dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Pastel2') + 
  coord_flip()


karyos %>% 
  filter(hgnc_symbol %in% hyper_res) %>% 
  mutate(karyotype = factor(karyotype, levels = rev(c('1:0', '1:1', '2:0', '2:1', '2:2')))) %>% 
  filter(karyotype != '1:1') %>% 
  group_by(hgnc_symbol) %>% 
  mutate(samples_per_gene = n()) %>% 
  group_by(hgnc_symbol, karyotype) %>% 
  summarise(prop = n()/samples_per_gene) %>% 
  distinct() %>% 
  dplyr::filter(karyotype == '1:0') %>% 
  ggplot(aes(
    y = hgnc_symbol, 
    x = prop
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw()
  

