library(tidyverse)

df = readRDS('data/compensation/cs_mean_across_karyos.rds')

df_groups = readRDS('data/compensation_score/sf_psinorm_stable_FALSE/CS_tables/cs_mean_groups.rds')

df_dna = df_groups %>% 
  left_join(., df) %>% 
  dplyr::select(-c(omic, lfc)) %>% 
  rename(DNA = DNA_lfc) %>% 
  pivot_longer(DNA, names_to = 'omic', values_to = 'lfc')

df_groups %>% 
  left_join(., df) %>%  
  dplyr::select(-DNA_lfc) %>% 
  bind_rows(., df_dna) %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) %>% 
  mutate(karyo_num = as.numeric(karyotype)) %>% 
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
ggsave('res/compensation_score/omics_fc_comparison_reg_groups.png', width = 12, height = 7)


# checking coad driver genes 

coad_genes = readRDS('data/utilities/COADREAD_drivers.rds') %>% 
  pull(SYMBOL) %>% unique

category_colors <- c(
  "(RNA-prot heavy)" = "#AD002AB2",
  "(Prot-heavy)" = "#E18727B2",
  "(RNA-heavy)" = "#20854Eb2",
  "(RNA-Prot light)" = "#00468BB2",
  "Intermediate/Other" = "gainsboro"
)


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




