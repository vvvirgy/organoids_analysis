library(tidyverse)
setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/constants.R')

# looking at intogen drivers

intogen_drivers = read.table('/orfeo/scratch/cdslab/vgazziero/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep = '\t', header = T)
coad_drivers = intogen_drivers %>% 
  dplyr::filter(CANCER_TYPE == 'COAD') %>% 
  dplyr::select(SYMBOL, IS_DRIVER) %>% 
  dplyr::distinct() %>% 
  pull(SYMBOL)

genes_to_plot = c(genes_to_plot, coad_drivers)

df_groups = readRDS('data/compensation/cs_classification.rds')

metadata = readRDS('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds')

# include also og and tsg
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

cancer_genes = read.table('/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/utilities/cancerGeneList.tsv', sep = "\t", header = F, skip = 1, col.names = cols) %>%
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

# add information on the reg groups 

df_groups = df_groups %>% 
  mutate(is_driver = ifelse(name %in% coad_drivers, 'Driver', 'Not a driver')) %>% 
  dplyr::full_join(., cancer_genes, by = join_by('name' == 'Hugo_Symbol')) %>% 
  dplyr::rename(role = 'Gene_Type') %>% 
  filter(!is.na(reg_group)) %>%
  mutate(status = paste(is_driver, role, sep = ', '))


tsg_og_reg_group_dist = df_groups %>% 
  filter(!is.na(role) | is_driver == 'Driver') %>%
  group_by(status,reg_group) %>% 
  summarise(n_genes = n()) %>% 
  group_by(status) %>% 
  mutate(tot = sum(n_genes), 
         frac = n_genes/tot) %>% 
  filter(reg_group != 'Intermediate/Other') %>% 
  ggplot(aes(
    y = status, 
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
  guides(fill = guide_legend(title = '')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        # legend.position = 'bottom',
        # legend.key.size = unit(.2, 'cm'), 
        legend.key.height = unit(.2, 'cm'),
        legend.key.width = unit(.2, 'cm'),
        legend.text = element_text(size=8)) + 
  scale_x_continuous(labels = scales::percent)

gene_prevalent_status = metadata %>% 
  group_by(hgnc_symbol) %>% 
  mutate(tot = n()) %>% 
  group_by(hgnc_symbol, mut_status) %>% 
  mutate(frac = n()/tot) %>% 
  dplyr::select(hgnc_symbol, mut_status, frac, tot) %>% 
  distinct() %>% 
  group_by(hgnc_symbol) %>%
  slice_max(frac) %>% 
  rename(prevalent_status = mut_status)

df_groups = df_groups %>% 
  right_join(., gene_prevalent_status, by = join_by('name' == 'hgnc_symbol')) %>% 
  filter(!is.na(reg_group))

tab <- df_groups %>% 
  ungroup() %>% 
  filter(!is.na(role) | is_driver == "Driver") %>%
  mutate(reg_group = ifelse(reg_group != 'Hyper-responders', 'other',reg_group )) %>% 
  # filter(reg_group != "Intermediate/Other") %>% 
  count(prevalent_status, reg_group) %>% 
  tidyr::pivot_wider(
    names_from = reg_group,
    values_from = n,
    values_fill = 0
  ) 

mat <- as.matrix(tab[,-1])
rownames(mat) <- tab$prevalent_status
mat

test = fisher.test(mat)
test$p.value

df_groups %>% 
  filter(!is.na(role) | is_driver == 'Driver') %>% 
  group_by(prevalent_status,reg_group) %>% 
  summarise(n_genes = n()) %>% 
  group_by(prevalent_status) %>% 
  mutate(tot = sum(n_genes), 
         frac = n_genes/tot) %>% 
  # filter(reg_group != 'Intermediate/Other') %>% 
  ggplot(aes(
    fill = prevalent_status, 
    x = frac,
    y = reg_group
  )) + 
  geom_bar(stat = 'identity', 
           # position = 'dodge',
           position = position_dodge2(preserve = "single", width = 1),
           alpha = 1, 
           # width = 
  ) +
  # ggplot(aes(prevalent_status, fill = reg_group)) + 
  # geom_bar(stat = 'count', position = 'dodge') + 
  # scale_fill_manual(values= category_colors) + 
  theme_bw() + 
  coord_flip() + 
  scale_x_continuous(labels = scales::percent) + 
  ggsignif::geom_signif(y_position = 0.1, 
                        xmin = 1.8, 
                        xmax = 2.2, 
                        annotation = signif(test$p.value, digits = 1), 
                        tip_length = 0.01)

  # facet_wrap(~reg_group, scales = 'free_y') + 
  


