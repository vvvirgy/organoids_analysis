library(tidyverse)
library(ggpmisc)
library(wesanderson)
library(ggExtra)
rm(list=ls())

rna_fit = readRDS('data/transcriptomics_lm_fit_ploidy_diff_multiplicity_tidy.rds')
prot_fit = readRDS('data/proteomics_lm_fit_ploidy_diff_multiplicity_tidy.rds')

rna_fit_res = readRDS('data/transcriptomics_lm_fit_ploidy_diff_multiplicity.rds')
prot_fit_res = readRDS('data/proteomics_lm_fit_ploidy_diff_multiplicity.rds')

rna_fit_metrics = lapply(names(rna_fit_res), function(x) {
  rna_fit_res[[x]] %>% 
    broom::glance() %>% 
    mutate(gene = x)
}) %>% 
  bind_rows()

prot_fit_metrics = lapply(names(prot_fit_res), function(x) {
  prot_fit_res[[x]] %>% 
    broom::glance() %>% 
    mutate(gene = x)
}) %>% 
  bind_rows()

all_fit = full_join(rna_fit, prot_fit, by = join_by('term' == 'term', 
                                                    'gene' == 'gene'), 
                    suffix = c('_RNA', '_prot'))

# look at genes present in both omics
all_fit_filt = all_fit %>% 
  dplyr::filter(!is.na(estimate_RNA)) %>% 
  dplyr::filter(!is.na(estimate_prot))

all_fit_filt %>% 
  dplyr::filter(term != '(Intercept)') %>%
  mutate(sign = case_when(
    (p.value_RNA <= 0.05 & p.value_prot <= 0.05) ~ 'both')) %>% 
  #   (p.value_RNA <= 0.05 & p.value_prot > 0.05) ~ 'RNA', 
  #   (p.value_RNA > 0.05 & p.value_prot <0 0.05) ~ 'Prot',
  #   .default = NA
  # )) %>% 
  ggplot(aes(x = estimate_RNA, y = estimate_prot)) + 
  geom_point() + 
  facet_wrap(~term) + 
  theme_bw() +
  geom_smooth(method = 'lm') + 
  ggpmisc::stat_poly_eq(use_label(c('P','eq'))) + 
  xlab('RNA coefficients') + 
  ylab('Protein coefficients')

ggsave('res/coefficients_comparison.png', width = 10, height = 6)

# visualizing coefficients distribution
plot_estimates_pvals(rna_fit, what = 'RNA', cols = c('#0D1164', 'grey60'))
ggsave('res/RNA_coefficients.png')

plot_estimates_pvals(prot_fit, what = 'Protein', cols = c('#0D1164', 'grey60')) 
ggsave('res/Protein_coefficients.png')


# classify the genes according to the betas direction 
comparing_rank_coefficients = all_fit %>% 
  dplyr::filter(term != '(Intercept)') %>%
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
    group_by(class) %>% 
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
  ggplot(aes(y = class, x = prop, fill = predictor)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(values = cols) + 
  theme_bw()
ggsave('res/comparison_type_betas_predictors.png')

rna_fit = rna_fit %>% 
  mutate(method = 'Transcriptomics')
prot_fit = prot_fit %>% 
  mutate(method = 'Proteomics')

method_colors = setNames(wes_palettes$Darjeeling1[2:3],
                         c('Transcriptomics', 'Proteomics'))

all_fit_v2 = bind_rows(rna_fit, prot_fit)
plot_predictor_higher_impact(all_fit_v2, palette = method_colors)

rna_betas = betas_association(rna_fit, color = method_colors[1], pth = 0.05) + 
  ggtitle('Betas estimation for transcriptomics')
ggsave(plot = rna_betas, filename =  'res/betas_rna.png') 
prot_betas = betas_association(prot_fit, color = method_colors[2], pth = 0.05) + 
  ggtitle('Betas estimation for proteomics')
ggsave(filename = 'res/betas_prot.png', prot_betas) 

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

