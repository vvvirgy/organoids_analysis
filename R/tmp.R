library(tidyverse)
library(decoupleR)

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/constants.R')
source('organoids_analysis/R/functions_utils/dge_utils_and_plots.R')


df_groups = readRDS('data/compensation/cs_classification.rds') 

data = df_groups %>% 
  group_by(name) %>% 
  mutate(mean_cs = mean(c(RNA, Protein))) %>% 
  filter(reg_group %in% c('Hyper-responders', 'Fully compensated')) %>% 
  dplyr::select(name, mean_cs, reg_group) %>% 
  pivot_wider(names_from = reg_group, values_from = mean_cs, values_fill = 0) %>% 
  tibble::column_to_rownames('name')

net = decoupleR::get_progeny(organism = 'human', top = 500)

# Run mlm
contrast_acts <- decoupleR::run_mlm(mat = data, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)

contrast_acts %>% 
  ggplot(aes(
    x = source, 
    y = score, 
    fill = condition
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() + 
  scale_fill_manual(values = category_colors)
  


# checking the alterations on specific genes

metadata = readRDS('data/processed_data/genes_filtered_karyo_mut_status_filt_ccf_08.rds')

metadata %>% 
  filter(hgnc_symbol == 'CDX2') %>% 
  ggplot(
    aes(karyotype, 
    fill = mut_consequence)
  ) + 
  geom_bar(position = 'dodge') + 
  theme_bw()


# correlations
df <- readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/results/sf_psinorm_stable_FALSE/lfc_prot_and_rna_bind.rds")

p1 = df %>% 
  group_by(name, karyotype) %>%
  filter(n() == 2) %>%
  # pull(name) %>% unique %>% length
  pivot_wider(names_from = omic, values_from = lfc) %>% 
  ggplot(aes(
    x = RNA, 
    y = Protein,
    color = karyotype
  )) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  facet_wrap(~karyotype, nrow = 1) + 
  labs(x = 'Log2FC RNA', 
       y = 'Log2FC Protein') + 
  geom_smooth(method = 'lm', show.legend = F, se = T, color = 'grey30', linetype = 'dashed') + 
  scale_color_brewer(palette = 'Set2')

data = df %>% 
  group_by(name, karyotype) %>%
  filter(n() == 2) %>% 
  pivot_wider(names_from = omic, values_from = lfc)

corr_res <- data %>%
  group_by(karyotype) %>%
    summarise(
      spearman_rho = cor(RNA, Protein, method = "spearman"),
      p_val = cor.test(RNA, Protein, method = "spearman")$p.value,
      n_genes = n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = p.adjust(p_val, "BH"))

p2 = corr_res %>% 
    ggplot(mapping = aes(x = karyotype, y = spearman_rho, fill = karyotype)) +
    geom_col(position = "dodge", width = .5) +
    scale_fill_brewer(palette = 'Set2') + 
    theme_bw() +
    labs(x = "Karyotype", y = "Spearman correlation")

pp = (p1 / p2) + 
  plot_annotation(tag_level = c('A'))
ggsave('/orfeo/scratch/cdslab/vgazziero/organoids_prj/res/compensation_score/corr_rna_prot.pdf', width = 10, height = 6)
ggsave('/orfeo/scratch/cdslab/vgazziero/organoids_prj/res/compensation_score/corr_rna_prot.png', width = 10, height = 6)

df = df %>% 
  dplyr::select(name, omic, karyotype, lfc)

write.csv(df, file = '/orfeo/scratch/cdslab/vgazziero/organoids_prj/res/genes_fc/dge_res.csv', quote = F, row.names = F, col.names = T)

df_groups = df_groups %>% 
  dplyr::select(-label)

write.csv(df_groups, file = '/orfeo/scratch/cdslab/vgazziero/organoids_prj/res/compensation_score/cs_classes_genes.csv', quote = F, row.names = F, col.names = T)

ann_res = cluster_comparison@compareClusterResult


write.csv(ann_res, file = '/orfeo/scratch/cdslab/vgazziero/organoids_prj/res/compensation_score/annotation_res_cs_groups.csv', quote = F, row.names = F, col.names = T)

# Create a label string
  # mutate(label = paste0("rho == ", round(rho, 2), "\nslope == ", round(slope, 2))) %>%
  # group_by(karyotype) %>%
  # mutate(y_pos = seq(max(data$Protein, na.rm=T), by = -0.5, length.out = n()))


