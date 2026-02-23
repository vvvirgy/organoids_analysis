library(tidyverse)
library(clusterProfiler)

essential_genes = read.table('data/homo_sapiens_essential_genes.csv', sep = ';', header = F)
essential_genes = essential_genes %>% 
  rename(gene_symbol = V3) %>%
  filter(V8 == 'Homo sapiens')
# essential_genes = readxl::read_excel('data/homo_sapiens_essential_genes.xlsx')
# essential_genes = essential_genes %>% 
#   rename(gene_function = `function`)

essential_genes_list = essential_genes$gene_symbol

multi_omics = readRDS('data/lfc_prot_and_rna_bind.rds')

# volcano plot of expression
pth = .05
fth = .75

multi_omics = classify_genes(multi_omics)
multi_omics = multi_omics %>% 
  mutate(omic = factor(omic, levels = c('RNA', 'protein')))

multi_omics = multi_omics %>% 
  dplyr::mutate(essential_gene = ifelse(name %in% essential_genes_list, TRUE, FALSE))

multi_omics %>% 
  filter(essential_gene) %>% 
  filter(significance == 'significant') %>% 
  plot_sankey(.,  strat = 'karyotype', facet = 'omic', cols = sankey_cols)

multi_omics %>% 
  filter(essential_gene) %>% 
  plot_barplots(cols = expr_cols) + 
  ylim(c(0,13))

multi_omics %>% 
  filter(essential_gene) %>% 
  pull(name) %>% 
  unique %>% length()


# housekeeping genes 

gmt = clusterProfiler::read.gmt('~/Downloads/msigdb.v2026.1.Hs.symbols.gmt')

housekeeping_genes = gmt %>% 
  filter(str_detect(pattern = 'HOUSEKEEPING', term)) %>% 
  pull(gene) %>% 
  unique

multi_omics = multi_omics %>% 
  mutate(housekeeping = ifelse(name %in% housekeeping_genes, TRUE, FALSE))

hk = multi_omics %>%
  filter(housekeeping) %>% 
  # filter(significance == 'significant') %>%
  plot_sankey(.,  strat = 'karyotype', facet = 'omic', cols = sankey_cols) +
  ggtitle('housekeeping genes')

ess = multi_omics %>%
  filter(essential_gene) %>% 
  # filter(significance == 'significant') %>%
  plot_sankey(.,  strat = 'karyotype', facet = 'omic', cols = sankey_cols) +
  ggtitle('Essential genes')

hk + ess

hk_barplot = multi_omics %>%
  filter(housekeeping) %>% 
  # filter(significance == 'significant') %>%
  plot_barplots(filter_significance = F, cols = expr_cols) + 
  ggtitle('housekeeping genes')


