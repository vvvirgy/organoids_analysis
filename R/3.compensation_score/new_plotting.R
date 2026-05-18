rm(list=ls())
library(tidyverse)
# library(decoupleR)

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/constants.R')
source('organoids_analysis/R/functions_utils/dge_utils_and_plots.R')


df_groups = readRDS('data/compensation/cs_classification.rds') 

# get drivers 
intogen_drivers = read.table('../2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep = '\t', header = T)
coad_drivers = intogen_drivers %>% 
  dplyr::filter(CANCER_TYPE == 'COAD') %>% 
  dplyr::select(SYMBOL, IS_DRIVER) %>% 
  dplyr::distinct() %>% 
  pull(SYMBOL)

hcc_drivers = intogen_drivers %>% 
  dplyr::filter(CANCER_TYPE == 'HCC') %>% 
  dplyr::select(SYMBOL, IS_DRIVER) %>% 
  dplyr::distinct() %>% 
  pull(SYMBOL)

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
  mutate(Gene_Type = factor(Gene_Type, levels = c('TSG', 'Oncogene', 'Oncogene/TSG'))) %>% 
  split(.$Gene_Type)
cancer_genes = lapply(cancer_genes, function(s) s$Hugo_Symbol)

# add information on the reg groups 

# df_groups = df_groups %>% 
#   mutate(is_driver = ifelse(name %in% coad_drivers, 'Driver', 'Not a driver')) %>% 
#   dplyr::full_join(., cancer_genes, by = join_by('name' == 'Hugo_Symbol')) %>% 
#   dplyr::rename(role = 'Gene_Type') %>% 
#   filter(!is.na(reg_group)) %>%
#   mutate(status = paste(is_driver, role, sep = ', '))

df_groups = df_groups %>% 
  mutate(coad_drivers = ifelse(name %in% coad_drivers, name, '')) %>% 
  mutate(hcc_drivers = ifelse(name %in% hcc_drivers, name, '')) %>% 
  mutate(TSG = ifelse(name %in% cancer_genes$TSG, name, '')) %>% 
  mutate(Oncogene = ifelse(name %in% cancer_genes$Oncogene, name, '')) %>% 
  mutate(Oncogene_TSG = ifelse(name %in% cancer_genes$`Oncogene/TSG`, name, '')) %>% 
  dplyr::select(-label) %>% 
  pivot_longer(cols = c('coad_drivers', 'hcc_drivers', 'TSG', 'Oncogene', 'Oncogene_TSG'), names_to = 'class', values_to = 'label')

df_groups %>% 
  ggplot(mapping = aes(x = RNA, y = Protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS", col = "") + 
  ggrepel::geom_label_repel(mapping = aes(label = label), show.legend = F, max.overlaps = Inf) + 
  theme(legend.position = 'bottom') + 
  guides(color = guide_legend(ncol = 3)) + 
  facet_wrap(~class)
  
saveRDS(df_groups, 'data/compensation/cs_classification_gene_labels.rds')

tsg_og_reg_group_dist = df_groups %>% 
  mutate(role = ifelse(label != '', class, NA)) %>% 
  filter(!is.na(role)) %>% 
  dplyr::select(-label) %>% 
  group_by(role,reg_group) %>% 
  summarise(n_genes = n()) %>% 
  group_by(role) %>% 
  mutate(tot = sum(n_genes), 
         frac = n_genes/tot) %>% 
  # filter(reg_group != 'Intermediate/Other') %>% 
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
  guides(fill = guide_legend(title = '')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        # legend.position = 'bottom',
        # legend.key.size = unit(.2, 'cm'), 
        legend.key.height = unit(.2, 'cm'),
        legend.key.width = unit(.2, 'cm'),
        legend.text = element_text(size=8)) + 
  scale_x_continuous(labels = scales::percent) 
