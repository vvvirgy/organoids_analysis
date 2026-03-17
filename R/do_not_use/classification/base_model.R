rm(list=ls())
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(mosaic)

# library(ggalluvial)
.libPaths()
setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/rna_vs_proteomics/utils.R')
source('organoids_analysis/R/rna_vs_proteomics/utils.R')

rna = readRDS('data/transcriptomics_data_all_genes_v2.rds') %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(rna_expression = value) %>% 
  dplyr::filter(!is.na(rna_expression))

prot = readRDS('data/proteogenomics_data_all_genes_new_norm.rds') %>%
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>%
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>%
  dplyr::rename(protein_expression = mean_intensity) %>%
  dplyr::filter(!is.na(protein_expression))

rna_expr_by_alleles = join_allele_expr(rna, column = 'rna_expression', which = 'Wild-type', use_2n = TRUE) 

prot_expr_by_alleles = join_allele_expr(prot, column = 'protein_expression', which = 'Wild-type', use_2n = TRUE) 

# test if any element from the distribution of non altered belong to the distribution of 2nwt samples
test_expr_rna = test_belonging_alterations(rna_expr_by_alleles, pth = 0.1, n_bootstrap = 10000)
test_expr_prot = test_belonging_alterations(prot_expr_by_alleles, pth = 0.1, n_bootstrap = 10000)

# classify each point
test_expr_rna = classify_elements(test_expr_rna) %>% 
  full_join(., rna_expr_by_alleles) %>% 
  mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
  rename(expression = rna_expression)
  
test_expr_prot = classify_elements(test_expr_prot) %>% 
  full_join(., prot_expr_by_alleles) %>% 
  mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
  rename(expression = protein_expression)

saveRDS(test_expr_rna, 'data/test_expr_rna_bootstrap.rds')
saveRDS(test_expr_prot, 'data/test_expr_prot_bootstrap.rds')

# trying to do some summary statistics and further analyses 

# do gene have a certain trend? --> proportions are changing between the two omics?
propotion_rna_by_gene = get_props(test_expr_rna, group = 'hgnc_symbol')
propotion_prot_by_gene = get_props(test_expr_prot, group = 'hgnc_symbol')

props_by_gene = full_join(propotion_rna_by_gene, propotion_prot_by_gene, suffix = c('_prot', '_rna'), 
                  by = join_by('hgnc_symbol' == 'hgnc_symbol', 
                               'cls_dosage' == 'cls_dosage'))

props_by_gene_plot = props_by_gene %>% 
  filter(!is.na(cls_dosage)) %>% 
  pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'assay', values_to = 'prop') %>% 
  mutate(assay = gsub('prop_', '', assay)) %>% 
  ggplot(aes(x = cls_dosage, y = prop, fill = assay)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'bottom') + 
  xlab('Dosage class') + 
  ylab('Proportion of genes') +
  scale_fill_manual(values = c('rna' = 'steelblue', 'prot' = 'goldenrod'))
ggsave('res/props_by_gene.png', props_by_gene_plot)

# do gene have a certain trend? --> for the same gene-sample couple, is the assigned class the same?
class_merge = merging_omics(prot_cls = test_expr_prot, rna_cls = test_expr_rna)

# plotting
alt_colors = setNames(
  # object = c('#4A70A9', '#658C58', '#B95E82', '#8D5F8C'),
  c('#476EAE', '#48B3AF', '#A7E399', '#F6FF99'), 
  nm = c("Wild-type, CNA", 'Mutated, CNA', 'Mutated, no CNA', 'Wild-type, no CNA')
  )

merge_plot = class_merge %>% 
  filter(!is.na(cls_comparison)) %>% 
  filter(alteration != 'Wild-type, no CNA') %>% 
  group_by(cls_comparison) %>% 
  count(alteration) %>% 
  # mutate(prop = n/sum(n)) %>% 
  ggplot(aes(
    y = cls_comparison, x = n, fill = alteration
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw()+ 
  scale_fill_manual(values = alt_colors) + 
  labs(x = 'Proportions', 
       y = 'Classes comparison')
merge_plot
saveRDS(merge_plot, 'data/plots_rds/cls_merge_plots.rds')
ggsave('res/cls_merge_plot.png', merge_plot)

# take a look to which are the genes with discordant classes
class_merge %>% 
  filter(cls_comparison == 'Different classes') %>% 
  ggplot(aes(
    cls_dosage_RNA, cls_dosage_Prot
  )) + 
  geom_jitter()

# alluvial plot 

class_merge %>% 
  filter(cls_comparison == 'Different classes') %>% 
  group_by(cls_dosage_RNA, cls_dosage_Prot) %>% 
  count(name = 'Freq') %>% 
  mutate(
    alluvium = interaction(cls_dosage_RNA, cls_dosage_Prot)
  ) %>% 
  pivot_longer(cols = starts_with('cls_'), values_to = 'class', names_to = 'assay') %>%
  mutate(assay = gsub('cls_dosage_', '', assay)) %>% 
  mutate(assay = factor(assay, levels = c('RNA', 'Prot'))) %>% 
  ggplot(aes(
    x = assay, 
    y = Freq, 
    # fill = class, 
    stratum = class, 
    alluvium = alluvium
  )) + 
  geom_flow() + 
  geom_stratum(aes(fill = class)) +
  theme_bw()
  
props %>% 
  filter(!is.na(cls_dosage)) %>% 
  ggplot(aes(prop_prot, prop_rna, color = cls_dosage, fill = cls_dosage)) + 
  geom_point() + 
  theme_bw() + 
  # geom_smooth(method = 'lm') + 
  facet_wrap(~cls_dosage) 
  
props_sample_prot = test_expr_prot %>% 
  filter(!is.na(sample)) %>% 
  group_by(sample) %>% 
  # mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
  count(cls_dosage) %>% 
  mutate(tot_genes = sum(n)) %>% 
  mutate(prop_genes = n/tot_genes)
  # select(-n, -tot_genes) %>% 
  # pivot_wider(names_from = cls_dosage, values_from = prop_genes) %>% 
  # filter(!is.na(sample)) %>% 
  # tibble::column_to_rownames('sample')

props_sample_prot %>% 
  mutate(cls_dosage = ifelse(is.na(cls_dosage), 'Not assigned', cls_dosage)) %>% 
  ggplot(aes(
    y = sample, x = prop_genes, fill = cls_dosage
  )) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme_bw() + 
  ggsci::scale_fill_bmj()

props_sample_rna = test_expr_rna %>% 
  filter(!is.na(sample)) %>% 
  group_by(sample) %>% 
  mutate(cls_dosage = ifelse((mutation_status == 'Wild-type' & tot_cna == 2), 'Control', cls_dosage)) %>% 
  count(cls_dosage) %>% 
  mutate(tot_genes = sum(n)) %>% 
  mutate(prop_genes = n/tot_genes)
# select(-n, -tot_genes) %>% 
# pivot_wider(names_from = cls_dosage, values_from = prop_genes) %>% 
# filter(!is.na(sample)) %>% 
# tibble::column_to_rownames('sample')

props_sample_rna %>% 
  mutate(cls_dosage = ifelse(is.na(cls_dosage), 'Not assigned', cls_dosage)) %>% 
  ggplot(aes(
    y = sample, x = prop_genes, fill = cls_dosage
  )) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme_bw() + 
  ggsci::scale_fill_bmj()


# -------------------------------------------------------------------------------------------
####################################################################################################################



# cluster_rna = kmeans_estimates(test_expr_rna)

  
# test_expr_rna %>% 
#   filter(hgnc_symbol == 'KRAS') %>% 
#   ggplot(aes(observed_expr)) + 
#   stat_function(fun = dnorm, 
#                 args = with(test_expr_rna, c(mean = mu_not_alt, sd = sd_not_alt))) 



################################################################################################################################################################## 
rna_expr_by_alleles %>% 
  mutate(alteration_classes = case_when(
    # (tot_cna == 2 & mutation_status != 'Wild-type') ~ 'Mutated, no CNA',
    # (tot_cna != 2 & mutation_status == 'Wild-type') ~ 'Wild-type, CNA',
    # (tot_cna != 2 & mutation_status != 'Wild-type') ~ 'Mutated, CNA', 
    (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'No alteration', 
    .default = 'Alteration'
  )) %>% 
  dplyr::filter(hgnc_symbol == 'AKNA') %>%
  ggplot(aes(observed_expr, fill = alteration)) + 
  geom_histogram() + 
  theme_bw() + 
  facet_grid(rows = vars(alteration_classes))



  # dplyr::filter(hgnc_symbol == 'KRAS') %>% 
  # mutate(alteration = ifelse((tot_cna == 2 & mutation_status == 'Wild-type'), 'Not altered', 'Altered')) 
rna = plot_single_allele_expr(rna_expr_by_alleles, gene = 'KRAS') + 
  ggtitle('KRAS, RNA')
prot = plot_single_allele_expr(prot_expr_by_alleles, gene = 'KRAS') + 
  ggtitle('KRAS, Protein')

(rna + prot) + 
  plot_layout(guides = 'collect') 
ggsave('res/kras_expr_by_allele.png', width = 10, height = 7)

# test 
rna_expr_by_alleles %>% 
  group_by(hgnc_symbol) %>% 
  mutate(class = case_when(
    (observed_expr > mean_wt_expr_observed+sd_expr) ~ 'Higher', 
    (observed_expr < mean_wt_expr_observed-sd_expr) ~ 'Lower', 
  )) %>% 
  group_by(hgnc_symbol, class, tot_cna) %>% 
  count()

# LET'S DO SOME PLOTTING
x = full_join(rna_expr_by_alleles, prot_expr_by_alleles, 
          by = join_by(
            'hgnc_symbol' == 'hgnc_symbol', 
            'chr' == 'chr', 
            'Major' == 'Major', 
            'minor' == 'minor', 
            'sample' == 'sample', 
            'is_mutated' == 'is_mutated', 
            'mut_consequence' == 'mut_consequence', 
            'driver_label' == 'driver_label', 
            'CGC_role_COAD' == 'CGC_role_COAD', 
            'CGC_role_PANCANCER' == 'CGC_role_PANCANCER',
            'is_driver_intogen' == 'is_driver_intogen', 
            'IMPACT' == 'IMPACT', 
            'tot_cna' == 'tot_cna', 
            'karyotype' == 'karyotype', 
            'mutation_status' == 'mutation_status', 
            'alteration' == 'alteration'
          ), 
          suffix = c('_RNA', '_Prot')
          ) %>% 
  mutate(alteration_classes = case_when(
    # (tot_cna == 2 & mutation_status != 'Wild-type') ~ 'Mutated, no CNA',
    # (tot_cna != 2 & mutation_status == 'Wild-type') ~ 'Wild-type, CNA',
    # (tot_cna != 2 & mutation_status != 'Wild-type') ~ 'Mutated, CNA', 
    (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'No alteration', 
    .default = 'Alteration'
  ))

alt_cols = setNames(
  nm = x$alteration %>% unique, 
  c('#1D546C', '#A72703', 'goldenrod', '', 'gainsboro')
)  

alphas = setNames(
  nm = x$alteration_classes %>% unique, 
  c(1, 0.1)
)

x %>% 
  filter(hgnc_symbol == 'TP53') %>% 
  ggplot(aes(observed_expr_RNA, observed_expr_Prot, color = alteration, alpha = alteration_classes)) + 
  geom_point() + 
  scale_color_manual(values = alt_cols) + 
  scale_alpha_manual(values = alphas) + 
  theme_bw() + 
  xlab('Observed_RNA/ploidy') + 
  ylab('Observed_Prot/ploidy')

x %>% 
  filter(hgnc_symbol == 'TP53') 

plot_single_allele_expr(rna_expr_by_alleles, gene = 'APC')

# try to test if the expression when present one alteration is different wrt expression wt

# rna %>% 
#   mutate(alteration = case_when(
#     (tot_cna == 2 & mutation_status == 'Wild-type') ~ 'Not altered',
#     .default = 'Altered'
#   )) %>% 
#   group_by(hgnc_symbol, alteration) %>% 
  



# df_long = df %>% 
#   pivot_longer(values_to = 'expression', cols = c(estimate_expression, protein_expression), names_to = 'origin') 
# 
# df_long %>% 
#   filter(origin == 'estimate_expression') %>% 
#   ggplot() + 
#   geom_point(aes(
#     y = expression, 
#     x = tot_cna
#   )) +
#   geom_line(aes(
#     y = expression, 
#     x = tot_cna
#   )) +
#   geom_point(data = df_long %>% 
#                filter(origin == 'protein_expression'), 
#              aes(
#                y = expression, 
#                x = tot_cna
#              ))





  
rna_expr_by_alleles %>%
  dplyr::filter(hgnc_symbol == 'KRAS') %>%
  ggplot(aes(mut_single_allele_expr)) +
  geom_rect(aes(xmin = mean_wt_expr-sd_expr,
                xmax = mean_wt_expr+sd_expr,
                ymax = unique(length(sample)),
                ymin = 0),
            alpha = 0.1,
            fill = 'skyblue4') +
  geom_histogram(binwidth = 0.1) +
  geom_vline(aes(xintercept = mean_wt_expr), linetype = 'dashed', show.legend = T, colour = 'firebrick') +
  theme_bw()

#   
# compute_expectation_expression(rna_2n_wt)


prot_2n_wt = prot %>% 
  # filter(karyotype == '1:1') %>%
  filter(tot_cna == 2) %>% 
  filter(mutation_status == 'Wild-type') %>% 
  mutate(single_allele_expr = protein_expression/2) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(n1 = mean(single_allele_expr), sd_expr = sd(single_allele_expr))

prot_2n_wt %>% 
  bind_cols(., 
            map_dfc(v, ~ prot_2n_wt$n1 * .x) %>%
              set_names(paste0("n", v))
  ) %>% 
  pivot_longer(cols = starts_with('n'), names_to = 'ploidy', values_to = 'expression') %>% 
  mutate(ploidy = gsub('n', '', ploidy))

