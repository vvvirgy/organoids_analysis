library(tidyverse)
library(ggpmisc)
library(wesanderson)
library(ggExtra)
library(factoextra)
library(kernlab)
rm(list=ls())
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/models_plotting_utils.R')
########### loading data
pth = 0.05
rna_fit = readRDS('data/transcriptomics_lm_fit_ploidy_diff_multiplicity_tidy.rds')
prot_fit = readRDS('data/proteomics_lm_fit_ploidy_diff_multiplicity_tidy.rds')

rna_fit_res = readRDS('data/transcriptomics_lm_fit_ploidy_diff_multiplicity.rds')
prot_fit_res = readRDS('data/proteomics_lm_fit_ploidy_diff_multiplicity.rds')

all_genes_functions = readRDS('data/all_genes_positions_info.rds') %>% 
  mutate(is_driver_intogen = ifelse(is_driver_intogen == 'True', TRUE, is_driver_intogen))

# rna_fit_metrics = lapply(names(rna_fit_res), function(x) {
#   rna_fit_res[[x]] %>% 
#     broom::glance() %>% 
#     mutate(gene = x)
# }) %>% 
#   bind_rows()
# 
# prot_fit_metrics = lapply(names(prot_fit_res), function(x) {
#   prot_fit_res[[x]] %>% 
#     broom::glance() %>% 
#     mutate(gene = x)
# }) %>% 
#   bind_rows()

all_fit = full_join(rna_fit, prot_fit, by = join_by('term' == 'term', 
                                                    'gene' == 'gene'), 
                    suffix = c('_RNA', '_prot'))

########### compare the omics
# look at genes present in both omics
all_fit_filt = all_fit %>% 
  dplyr::filter(!is.na(estimate_RNA)) %>% 
  dplyr::filter(!is.na(estimate_prot))

all_fit_v2 = all_fit %>% 
  dplyr::filter(term != '(Intercept)') %>%
  mutate(Significance = case_when(
    (p.value_RNA <= 0.05 & p.value_prot <= 0.05) ~ 'both', #) %>% 
    (p.value_RNA <= 0.05 & p.value_prot > 0.05 | p.value_RNA <= 0.05 & is.na(p.value_prot) ) ~ 'RNA',
    (p.value_RNA > 0.05 & p.value_prot < 0.05 | is.na(p.value_RNA)  & p.value_prot <= 0.05) ~ 'Prot', 
    (p.value_RNA > 0.05 & p.value_prot > 0.05) ~ 'None', 
    .default = 'None'
  ))

coeff_comparison = all_fit_v2 %>% 
  #   .default = NA
  # )) %>% 
  ggplot(aes(x = estimate_RNA, y = estimate_prot, color = Significance)) + 
  geom_point() + 
  facet_grid(term~Significance) + 
  theme_bw() +
  geom_smooth(method = 'lm', aes(fill = Significance), inherit.aes = T, alpha = 0.1) + 
  ggpmisc::stat_poly_eq(use_label(c('P','eq'))) + 
  ggsci::scale_color_futurama() +
  ggsci::scale_fill_futurama() +
  xlab('RNA coefficients') + 
  ylab('Protein coefficients') + 
  theme(legend.position = 'bottom') + 
  guides(color = guide_legend(title = 'Significance (P <= 0.05)'), fill = guide_legend(title = 'Significance (P <= 0.05)'))

ggsave(plot = coeff_comparison, filename = 'res/coefficients_comparison.png', width = 10, height = 6)

# visualizing coefficients distribution
plot_estimates_pvals(rna_fit, what = 'RNA', cols = c('#0D1164', 'grey60'), pval_th = pth)
ggsave('res/RNA_coefficients.png')

plot_estimates_pvals(prot_fit, what = 'Protein', cols = c('#0D1164', 'grey60'), pval_th = pth) 
ggsave('res/Protein_coefficients.png')

all_fit_v2 %>% 
  ggplot(aes(estimate_RNA, fill = Significance)) + 
  geom_histogram(binwidth = 0.1) +
  # geom_density() +
  facet_grid(Significance~term, scales = 'free_y') + 
  theme_bw()

# classify the genes according to the betas direction 
comparing_rank_coefficients = all_fit_v2 %>% 
  group_by(term) %>% 
  group_split()
 
names(comparing_rank_coefficients) = lapply(comparing_rank_coefficients, function(x) {x$term %>% unique}) %>% unlist

comparing_rank_coefficients = lapply(comparing_rank_coefficients, function(x) {
  x %>% 
    mutate(class = 
             case_when((estimate_RNA < 0 & estimate_prot < 0) ~ 'both negative', 
                       (estimate_RNA > 0 & estimate_prot > 0) ~ 'both positive', 
                       (estimate_RNA > 0 & estimate_prot < 0) ~ 'RNA pos, protein neg', 
                       (estimate_RNA < 0 & estimate_prot > 0) ~ 'RNA neg, protein pos',
                       (estimate_RNA == 0 & estimate_prot == 0) ~ 'both zero', 
                       (is.na(estimate_RNA) & !is.na(estimate_prot)) ~ 'only prot', 
                       (!is.na(estimate_RNA) & is.na(estimate_prot)) ~ 'only rna', 
                       .default = 'no estimate')
    )
})

comparison = lapply(comparing_rank_coefficients, function(x) {
  x %>% 
    mutate(tot = n()) %>% 
    group_by(class, Significance) %>% 
    count() %>% 
    ungroup()
})
saveRDS(comparing_rank_coefficients, 'data/comparing_rank_coefficients.rds')
# TODO plot class coefficients

comparison = lapply(comparison %>% names, function(x) {
  comparison[[x]] %>% 
    mutate(tot = sum(n)) %>% 
    mutate(prop = n/tot) %>% 
    mutate(predictor = x)
}) %>% 
  bind_rows()

# cols = c('mutation_multiplicity' = '#9E1C60', 'ploidy_diff' = '#636CCB')
cols = c('mutation_multiplicity' = '#B3C8CF', 'ploidy_diff' = '#6C567B')

# '#99235C'
comparison %>% 
  ggplot(aes(y = class, x = prop, fill = Significance)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  # scale_fill_manual(values = cols) + 
  ggsci::scale_fill_futurama() + 
  theme_bw() +
  facet_wrap(~predictor)
  # facet_wrap(~Significance)
ggsave('res/comparison_type_betas_predictors.png')

rna_fit = rna_fit %>% 
  mutate(method = 'Transcriptomics')
prot_fit = prot_fit %>% 
  mutate(method = 'Proteomics')

method_colors = setNames(wes_palettes$Darjeeling1[2:3],
                         c('Transcriptomics', 'Proteomics'))

pp = bind_rows(rna_fit, prot_fit) %>% 
  plot_predictor_higher_impact(., palette = method_colors, pth = 0.05) + 
  theme(legend.position = 'bottom')
ggsave('res/higher_impact_pred_significant.png', width = 5, height = 6)

des = 'AAAB'
(coeff_comparison + pp ) + 
  plot_layout(design = des)
ggsave('res/patch_res.png', width = 12, height = 6)

rna_betas = betas_association(rna_fit, pth = 0.05, cols = cols) + 
  ggtitle('Betas estimation for transcriptomics \n(Significant only)') +
  # xlim(c(-5.5,5.5)) + 
  ylim(c(-4.5,4.5))
ggsave(plot = rna_betas, filename =  'res/betas_rna.png') 
prot_betas = betas_association(prot_fit, pth = 0.05, cols = cols) + 
  ggtitle('Betas estimation for proteomics \n(Significant only)')  +
  # xlim(c(-5.5,5.5)) + 
  ylim(c(-4.5,4.5))
ggsave(filename = 'res/betas_prot.png', prot_betas) 

(rna_betas + prot_betas) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
ggsave(filename = 'res/betas_prot_rna_sign.png', width = 10, height = 6) 

cls_betas_rna = classify_betas(rna_fit, pth = 0.05)
cls_betas_prot = classify_betas(prot_fit, pth = 0.05)

cls_genes = all_genes_functions %>% 
  dplyr::select(hgnc_symbol, is_driver_intogen, CGC_role_COAD, CGC_role_PANCANCER) %>%
  distinct() 

cls_betas_rna = cls_betas_rna %>% 
  left_join(., cls_genes, by = join_by('gene' == 'hgnc_symbol'))
cls_betas_prot = cls_betas_prot %>% 
  left_join(., cls_genes, by = join_by('gene' == 'hgnc_symbol'))

cls_betas = list('RNA' = cls_betas_rna, 'Prot' = cls_betas_prot)

genes_types = c('CGC_role_COAD', 'CGC_role_PANCANCER')
cols_gene_role = setNames(nm = 
                  c('fusion', 'mixed', 'TSG', 'oncogene'), 
                object = wes_palettes$GrandBudapest1
                )

cls_betas = bind_rows(cls_betas)
cls_betas %>% 
  select(gene, method, class, CGC_role_COAD, CGC_role_PANCANCER) %>% 
  # filter(!is.na(CGC_role_COAD)) %>% 
  # filter(!is.na(CGC_role_PANCANCER)) %>% 
  # filter(CGC_role_COAD != 'None') %>%
  # filter(CGC_role_PANCANCER != 'None') %>%
  # filter(CGC_role_COAD != FALSE) %>%
  # filter(CGC_role_PANCANCER != FALSE) %>%
  tidyr::pivot_longer(cols = c(CGC_role_COAD, CGC_role_PANCANCER), names_to = 'Database', values_to = 'Gene_role') %>% 
  filter(Gene_role != 'None') %>% 
  filter(!is.na(Gene_role)) %>% 
  ggplot(aes(class, fill = Gene_role)) + 
  geom_bar(stat = 'count', position = 'dodge') + 
  coord_flip() + 
  scale_fill_manual(values = cols) +
  theme_bw() + 
  # guides(fill = guide_legend(title = '')) + 
  theme(legend.direction = 'horizontal') + 
  facet_grid(method~Database) + 
  theme(legend.position = 'bottom')
ggsave('res/cgc_gene_roles_effect.png', width = 8, height = 6)

# cls_betas_rna_all = classify_betas(rna_fit, pth = 1)
# cls_betas_prot_all = classify_betas(prot_fit, pth = 1)

cls_genes = all_genes_functions %>% 
  dplyr::select(hgnc_symbol, is_driver_intogen, CGC_role_COAD, CGC_role_PANCANCER) %>%
  distinct() 

cls_betas_rna = rna_fit %>% 
  filter(term != '(Intercept)') %>% 
  left_join(., cls_genes, by = join_by('gene' == 'hgnc_symbol'))
cls_betas_prot = prot_fit %>% 
  filter(term != '(Intercept)') %>% 
  left_join(., cls_genes, by = join_by('gene' == 'hgnc_symbol'))

all_betas = bind_rows(cls_betas_rna, cls_betas_prot)
  
all_betas %>% 
  filter(CGC_role_PANCANCER != 'None') %>%
  filter(method == 'Transcriptomics') %>% 
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(sign = ifelse(p.value <= pth, 
                              paste0('Significant (P <= ', pth, ')'),
                              paste0('Not significant (P > ', pth, ')'))) %>% 
  ggplot(aes(y = -log10(p.value), x = estimate, alpha = sign, colour = CGC_role_PANCANCER)) + 
  geom_point() + 
  facet_grid(term~CGC_role_PANCANCER) + 
  theme_bw() + 
  scale_color_manual(values = cols_gene_role) +
  ggtitle('Effect on expression for transcriptomics') + 
  geom_hline(yintercept = -log10(pth), 
             linetype = 'dashed', 
             colour = 'Firebrick', 
             show.legend = T, 
             linewidth = 1
  )
ggsave('res/expression_effect_rna.png', width = 8, height = 5)

all_betas %>% 
  filter(CGC_role_PANCANCER != 'None') %>%
  filter(method == 'Proteomics') %>% 
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(sign = ifelse(p.value <= pth, 
                              paste0('Significant (P <= ', pth, ')'),
                              paste0('Not significant (P > ', pth, ')'))) %>% 
  ggplot(aes(y = -log10(p.value), x = estimate, alpha = sign, colour = CGC_role_PANCANCER)) + 
  geom_point() + 
  facet_grid(term~CGC_role_PANCANCER) + 
  theme_bw() + 
  scale_color_manual(values = cols_gene_role) +
  ggtitle('Effect on expression for proteomics') + 
  geom_hline(yintercept = -log10(pth), 
             linetype = 'dashed', 
             colour = 'Firebrick', 
             show.legend = T, 
             linewidth = 1
  )
ggsave('res/expression_effect_prot.png', width = 8, height = 5)


transcriptomics_data %>% 
  filter(hgnc_symbol == 'CDH1') %>% 
  mutate(ploidy_diff = tot_cna -2) %>%
  mutate(ploidy_diff = factor(ploidy_diff)) %>% 
  ggplot(aes(x = mutation_multiplicity, y = rna_expression)) + 
  geom_point() + 
  geom_smooth(method = 'lm')



proteogenomics_data %>% 
  filter(hgnc_symbol == 'MAP2K1') %>% 
  ggplot(aes(x = mutation_multiplicity, y = protein_expression)) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  theme_bw() + 
  ggtitle('MAP2K1 expression (Protein)')
ggsave('res/MAP2K1.png')

transcriptomics_data %>% 
  filter(hgnc_symbol == 'RXRA') %>% 
  mutate(ploidy_diff = tot_cna -2) %>%
  # mutate(ploidy_diff = factor(ploidy_diff)) %>% 
  ggplot(aes(x = ploidy_diff, y = rna_expression)) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  theme_bw() + 
  ggtitle('RXRA expression (RNA)') 
ggsave('res/RXRA.png')

cls_betas_plt = lapply(cls_betas %>% names, function(x) {
  lapply(genes_types, function(g) {
    cls_betas[[x]] %>%
      filter(!is.na(.data[[g]])) %>% 
      filter(.data[[g]] != 'None') %>%
      filter(.data[[g]] != FALSE) %>% 
      rename(Gene_role = .data[[g]]) %>% 
      ggplot(aes(class, fill = Gene_role)) + 
      geom_bar(stat = 'count') + 
      coord_flip() + 
      scale_fill_manual(values = cols) +
      theme_bw() + 
      ggtitle(g) + 
      guides(fill = guide_legend(title = '')) + 
      theme(legend.direction = 'horizontal')
  }) %>% 
    wrap_plots(guides = 'collect') + 
    plot_annotation(title = x) & 
    theme(legend.position = 'bottom')
})

ggsave(plot = cls_betas_plt[[1]], 'res/cgc_gene_roles_effect_rna.pdf', width = 10)
ggsave(plot = cls_betas_plt[[2]], 'res/cgc_gene_roles_effect_prot.pdf', width = 10)

#+ 
  # facet_wrap(~CGC_role_PANCANCER, scales = 'free_x')

# perform clustering on significant coefficients
# all_fit_res = list('RNA' = rna_fit, 'Protein' = prot_fit)
# 
# x = all_fit_res$RNA %>% 
#   filter(term != '(Intercept)') %>%
#   dplyr::filter(p.value <= pth) %>%
#   select(gene, term, estimate) %>%
#   tidyr::pivot_wider(names_from = term, values_from = estimate, values_fill = 0) %>%
#   tibble::column_to_rownames('gene')
# 
# # BIC <- mclustBIC(x)
# mod1 <- Mclust(x, modelNames = 'VVI', G = 1:5)
# plot(mod1)
# 
# 
# dist_mat = dist(x)
# cls_betas = kernlab::specc(x = dist_mat %>% as.matrix, centers = 5)
# 
# classes_betas = x %>% 
#   mutate(class = cls_betas@.Data)

  

# checking for the distribution conditioned on the pvalues
fit_clusters = lapply(all_fit_res, function(x){
  kmeans_estimates(x, seed = 67921)
})


clusters_2d = lapply(fit_clusters %>% names, function(x) {
  plot_clusters_2D(fit_clusters[[x]]) +
    ggtitle(x)
}) %>% 
  wrap_plots() &
  theme(legend.position = 'bottom')

x = lapply(names(fit_clusters), function(x) {
    fit_clusters[[x]][['clusters']] %>% 
      mutate(group = paste(x, sep = '_'))
}) %>% 
  bind_rows() %>% 
  mutate(cluster = paste0('C', cluster)) %>% 
  group_by(cluster, group) %>% 
  group_split()
names(x) = lapply(x, function(s) {paste(unique(s$group), unique(s$cluster), sep = '_')}) %>% unlist
x = lapply(x, function(s) {s %>% pull(gene)})

wrap_plots(
  list(coeff_comparison, clusters_2d), ncol = 1
)








# get the gene names for genes that have higher impact on the mutation multiplicity
genes_mutations = all_fit_v2 %>% 
  dplyr::filter(term != '(Intercept)') %>% 
  dplyr::select(term, gene, estimate, p.value, method) %>% 
  tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>% 
  mutate(higher_predictor = case_when(
    abs(estimate_ploidy_diff) > abs(estimate_mutation_multiplicity) ~ 'ploidy_diff', 
    abs(estimate_ploidy_diff) < abs(estimate_mutation_multiplicity) ~ 'mutation_multiplicity', 
    (is.na(estimate_ploidy_diff) & !is.na(estimate_mutation_multiplicity)) ~ 'mutation_multiplicity', 
    (!is.na(estimate_ploidy_diff) & is.na(estimate_mutation_multiplicity)) ~ 'ploidy_diff', 
    .default = NA
  )) %>% 
  filter(higher_predictor == 'mutation_multiplicity') %>% 
  # filter(p.value_mutation_multiplicity <= 0.05) %>% 
  dplyr::select(gene, method, estimate_mutation_multiplicity, p.value_mutation_multiplicity) %>%
  group_by(method) %>% 
  group_split()
names(genes_mutations) = lapply(genes_mutations, function(x) {x$method %>% unique}) %>% unlist

genes_muts_rna = transcriptomics_data %>% 
  dplyr::filter(hgnc_symbol %in% (genes_mutations$Transcriptomics %>% pull(gene))) %>% 
  dplyr::select(hgnc_symbol, sample, mutation_status, mutation_multiplicity)

genes_muts_rna %>% 
  group_by(hgnc_symbol) %>% 
  count(mutation_status) %>% 
  tidyr::pivot_wider(names_from = mutation_status, values_from = n) %>% 
  dplyr::full_join(., genes_mutations$Transcriptomics, by = join_by('hgnc_symbol' == 'gene')) %>% 
  ggplot(aes(estimate_mutation_multiplicity)) + 
  geom_histogram()

genes_muts_prot = proteogenomics_data %>% 
  dplyr::filter(hgnc_symbol %in% (genes_mutations$Proteomics %>% pull(gene))) %>% 
  dplyr::select(hgnc_symbol, sample, mutation_status, mutation_multiplicity)

rna_drivers = transcriptomics_data %>% 
  filter(is_driver_intogen == 'True') %>% 
  pull(hgnc_symbol) %>% 
  unique
prot_drivers = proteogenomics_data %>% 
  filter(is_driver_intogen == 'True') %>% 
  pull(hgnc_symbol) %>% 
  unique

rna_fit %>% 
  dplyr::filter(gene %in% rna_drivers) %>% 
  # filter(p.value <= 0.1) %>%
  dplyr::filter(term != '(Intercept)') %>% 
  mutate(class = ifelse(estimate > 0, 'positive', 'negative')) %>% 
  dplyr::select(term, gene, class, p.value, method) %>% 
  tidyr::pivot_wider(names_from = term, values_from = c(class, p.value)) %>% 
  select(gene, starts_with('class')) %>% 
  column_to_rownames('gene') %>% 
  as.matrix %>% 
  ComplexHeatmap::Heatmap()


all_fit_v2 %>% 
  dplyr::filter(term != '(Intercept)') %>% 
  ggplot(aes(x = estimate, fill = term)) + 
  geom_histogram() + 
  facet_wrap(~method) + 
  theme_bw()

# all_fit_v2 %>% 
#   filter(term != '(Intercept)') %>% 
#   select(gene, method, term, estimate) %>% 
#   pivot_wider(names_from = method, values_from = estimate) %>% 
#   mutate(concordant = sign(Transcriptomics) == sign(Proteomics)) %>%
#   count(term, concordant) %>%
#   ggplot(aes(x = term, y = n, fill = concordant)) +
#   geom_col(position = "fill") +
#   scale_y_continuous(labels = scales::percent) +
#   labs(y = "% of genes", title = "Sign concordance of RNA vs Protein β’s")

