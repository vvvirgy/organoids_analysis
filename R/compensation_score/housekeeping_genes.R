rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
library(tidyverse)
library(clusterProfiler)

essential_genes = read.table('data/utilities/homo_sapiens_essential_genes.csv', sep = ';', header = F)
essential_genes = essential_genes %>% 
  rename(gene_symbol = V3) %>%
  filter(V8 == 'Homo sapiens')
# essential_genes = readxl::read_excel('data/homo_sapiens_essential_genes.xlsx')
# essential_genes = essential_genes %>% 
#   rename(gene_function = `function`)

essential_genes_list = essential_genes$gene_symbol

# check the cs of these genes across karyotypes
df = readRDS('data/compensation/cs_mean_across_karyos.rds') %>% 
  select(-c(DNA_lfc, n, mean_CS))

karyos = c('1:0', '2:0', '2:1', '2:2')

karyos_groups = lapply(karyos, function(x) {
  define_th(df, q = .5, karyo = x)  
})

karyos_groups = karyos_groups %>% 
  bind_rows()

# plot proportions
plot_props_by_karyo(karyos_groups, colors = category_colors, gene_list = essential_genes_list)


df_groups = readRDS('data/compensation_score/sf_psinorm_stable_FALSE/CS_tables/cs_mean_groups.rds')

plot_props(df_groups, colors = category_colors, gene_list = essential_genes_list)

df_groups %>% 
  filter(name %in% essential_genes_list) %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") +
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2, 
                              title.position = 'top', 
                              title.hjust =0.5)) + 
  xlim(-3,4) + 
  ylim(-3,4)

# checking homeostasis genes 
homeostatis_genes = clusterProfiler::read.gmt('data/utilities/GOBP_HOMEOSTATIC_PROCESS.v2026.1.Hs.gmt')
homeostatis_genes = homeostatis_genes$gene %>% unique

plot_props_by_karyo(karyos_groups, colors = category_colors, gene_list = homeostatis_genes)

df_groups %>% 
  filter(name %in% homeostatis_genes) %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") +
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2, 
                              title.position = 'top', 
                              title.hjust =0.5)) + 
  xlim(-3,4) + 
  ylim(-3,4)


dm_essential = read.table('data/utilities/deepmap_colon_essential_genes', header = F, sep = '\t')
dm_essential = dm_essential$V1 %>% unique

plot_props_by_karyo(karyos_groups, colors = category_colors, gene_list = dm_essential)


df_groups %>% 
  filter(name %in% dm_essential) %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") +
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2, 
                              title.position = 'top', 
                              title.hjust =0.5)) + 
  xlim(-3,4) + 
  ylim(-3,4)

dm = read.table('data/utilities/dm_v2.txt', header = F, sep = '\t')
dm = dm$V1 %>% unique
df_groups %>% 
  filter(name %in% dm) %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") +
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2, 
                              title.position = 'top', 
                              title.hjust =0.5)) + 
  xlim(-3,4) + 
  ylim(-3,4)

################################################################################


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



