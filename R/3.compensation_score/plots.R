setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj')
rm(list=ls())
.libPaths()
library(tidyverse)
library(patchwork)
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/colors.R')
source('organoids_analysis/jovoniR/constants.R')

df_groups = readRDS('data/compensation/cs_classification_min_1.rds')

df = readRDS('data/compensation/cs_classification_by_karyo_min_1.rds')

df_v2 = df_groups %>% 
  rename(reg_group_mean = reg_group) %>%
  dplyr::select(-starts_with('DNA')) %>% 
  left_join(., df) %>% 
  dplyr::select(-c(RNA, protein)) %>% 
  # dplyr::select(-lfc) %>% 
  rename(DNA = DNA_lfc) %>%
  rename(RNA = lfc_RNA) %>% 
  rename(Protein = lfc_protein) %>% 
  pivot_longer(cols = c(RNA, DNA, Protein), names_to = 'omic', values_to = 'lfc')

groups_plt = df_v2 %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
  ggplot(aes(
    y = lfc, 
    x = karyotype, 
    fill = omic
  )) + 
  geom_boxplot(position = 'dodge', outliers = T) + 
  geom_smooth(
    aes(
      x = karyo_num, 
      group = omic, 
      color = omic, 
      y = lfc, 
      fill = omic
    ),
    method = "lm",  
    # se = FALSE,
    inherit.aes = FALSE, show.legend = F, 
    alpha = .3, linetype = 1
  ) +
  theme_bw() +
  facet_grid(~reg_group) + 
  scale_fill_manual(values = omic_colors) + 
  scale_color_manual(values = omic_colors) + 
  # scale_fill_manual(values = c('DNA' = '#237227', 'RNA' = '#088395', 'protein' = '#FFAA00')) + 
  # scale_color_manual(values = c('DNA' = '#237227', 'RNA' = '#088395', 'protein' = '#FFAA00')) + 
  guides(color = guide_legend(title = 'Omic'), 
         fill = guide_legend(title = 'Omic')) + 
  labs(x = 'Karyotype', y = 'Log2FC')
ggsave('res/compensation_score/omics_fc_comparison_reg_groups_v2.png', plot = groups_plt, width = 12, height = 7)
ggsave('res/compensation_score/omics_fc_comparison_reg_groups_v2.pdf', plot = groups_plt, width = 12, height = 7)

# checking specific gene sets

# checking on TSG and OG

# cancer_genes = read.table('data/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
#   dplyr::as_tibble()
# 
# # genes associated with COAD
# cancer_genes_somatic_colon = cancer_genes %>% 
#   dplyr::filter(Somatic == 'yes') %>% 
#   dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
#   dplyr::select(Gene.Symbol, Role.in.Cancer) %>% 
#   dplyr::distinct() %>% 
#   mutate(role = case_when(
#     role == 'oncogene, fusion' ~ 'oncogene', 
#     role == 'TSG, fusion' ~ 'TSG',
#     
#   ))

cols = c('Hugo_Symbol', 	
        'Entrez_Gene_ID', 
        'GRCh37_Isoform', 	
        'GRCh37_RefSeq',	
        'GRCh38_Isoform', 	
        'GRCh38_RefSeq', 	
        'Gene_Type', 
        'number_occurrences', 
        'OncoKB_Annotated',
        'MSK_IMPACT', 
        'MSK_HEME',
        'FOUNDATION_ONE',
        'FOUNDATION_ONE_HEME',
        'Vogelstein',
        'COSMIC_CGC',
        'Gene_Aliases')

cancer_genes = read.table('data/utilities/cancerGeneList.tsv', sep = "\t", header = F, skip = 1, col.names = cols) %>%
  dplyr::as_tibble() %>% 
  filter(Gene_Type %in% c('ONCOGENE', 'TSG', 'ONCOGENE_AND_TSG')) %>% 
  dplyr::select(Hugo_Symbol, Gene_Type) %>% 
  distinct() %>% 
  mutate(Gene_Type = case_when(
    Gene_Type == 'ONCOGENE' ~ 'Oncogene', 
    Gene_Type == 'ONCOGENE_AND_TSG' ~ 'Oncogene/TSG', 
    .default = Gene_Type
  )) %>% 
  mutate(Gene_Type = factor(Gene_Type, levels = c('TSG', 'Oncogene', 'Oncogene/TSG')))

# str_to_title(cancer_genes$Gene_Type)

gene_roles_groups = df_groups %>% 
  dplyr::right_join(., cancer_genes, by = join_by('name' == 'Hugo_Symbol')) %>% 
  rename(role = Gene_Type) %>% 
  filter(!is.na(reg_group)) 


gene_roles_groups %>% 
  filter(role == 'Oncogene') %>% 
  filter(reg_group == '(RNA-Prot light)') %>% 
  pull(name) %>% unique


gene_roles_groups %>% 
  group_by(role,reg_group) %>% 
  summarise(n_genes = n()) %>% 
  group_by(role) %>% 
  mutate(tot = sum(n_genes), 
         frac = n_genes/tot) %>% 
  filter(reg_group != 'Intermediate/Other') %>% 
  ggplot(aes(
    y = role, 
    x = frac,
    fill = reg_group
  )) + 
  geom_bar(stat = 'identity', 
           # position = 'dodge',
           position = position_dodge2(preserve = "single", width = 1),
           alpha = 1, 
           # width = 
           ) + 
  theme_bw()  +
  scale_fill_manual(values = category_colors) + 
  coord_flip() +
  labs(
    y = 'Gene role', 
    x = '% genes'
  ) + 
  guides(fill = guide_legend(title = 'Regulatory group')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_x_continuous(labels = scales::percent)



df %>% 
  filter(name %in% c('MERTK', 'KRAS', 'PTK7', 'FYN', 'ABCB1', 'CDX2', 'APOBEC3C')) %>% 
  ggplot(aes(
    x = karyotype, 
    y = lfc, 
    fill = omic
  )) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  theme_bw() + 
  scale_fill_manual(values = omic_colors) + 
  facet_wrap(~name, nrow = 1) 

# checking drivers
intogen_drivers = read.table('../2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep = '\t', header = T)
coad_drivers = intogen_drivers %>% 
  dplyr::filter(CANCER_TYPE == 'COAD') %>% 
  dplyr::select(SYMBOL, IS_DRIVER) %>% 
  dplyr::distinct() %>% 
  pull(SYMBOL)

df_groups












# toy example
toy_plt = readRDS('data/figures_data/cs/toy_plot.rds')

omics_scatter = df_groups %>% 
  mutate(gene_id = ifelse(reg_group != "Intermediate/Other", name, '')) %>% 
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

pt_layout = '
AABB
CCCC
'
wrap_plots(list(toy_plt, omics_scatter, groups_plt), design = pt_layout) + 
  plot_annotation(tag_levels = 'A')
ggsave('figures/cs/fig1.pdf', width = 10, height = 10)
ggsave('figures/cs/fig1.png', width = 10, height = 10)


# checking coad driver genes 

coad_genes = readRDS('data/utilities/COADREAD_drivers.rds') %>% 
  pull(SYMBOL) %>% unique


df_groups %>% 
  filter(name %in% coad_genes) %>% 
  mutate(gene_id = ifelse(reg_group != "Intermediate/Other", name, '')) %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group, label = gene_id)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") +
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2)) + 
  ggrepel::geom_label_repel(show.legend = F)


df_groups %>% 
  left_join(., df) %>%  
  dplyr::select(-DNA_lfc) %>% 
  bind_rows(., df_dna) %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
  filter(name %in% coad_genes) %>% 
  ggplot(aes(
    y = lfc, 
    x = karyotype, 
    fill = omic
  )) + 
  geom_boxplot(position = 'dodge') + 
  # theme_bw() + 
  # facet_grid(~reg_group, scales = 'free')
  geom_smooth(
    aes(
      x = karyo_num, 
      group = omic, 
      color = omic, 
      y = lfc
    ),
    method = "lm",        # or "loess"
    # se = FALSE,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  facet_grid(~reg_group, scales = "free") + 
  scale_fill_manual(values = c('DNA' = '#537D96', 'RNA' = '#44A194', 'protein' = '#F26076')) + 
  scale_color_manual(values = c('DNA' = '#537D96', 'RNA' = '#44A194', 'protein' = '#F26076'))


df_groups %>% 
  left_join(., df) %>%  
  dplyr::select(-DNA_lfc) %>% 
  bind_rows(., df_dna) %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
  filter(name == 'PTEN') %>% 
  ggplot(aes(
    x = karyotype, 
    y = lfc, 
    fill = omic
  )) + 
  geom_bar(stat = 'identity', position = 'dodge')


df_groups %>% 
  left_join(., df) %>%  
  dplyr::select(-DNA_lfc) %>% 
  bind_rows(., df_dna) %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
  filter(name %in% coad_genes) %>%
  # dplyr::select(name, karyotype, n) %>% 
  ggplot(aes(
    y = karyotype, 
    x = lfc, 
    fill = name
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_wrap(~omic, nrow= 1)




