library(tidyverse)
library(ggpmisc)
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')
source('organoids_analysis/R/functions_utils/fit_plots.R')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
# transcriptomics_data_low_CNAs = transcriptomics_data %>% 
#   dplyr::filter(tot_cna <= 6)

# fit to extract the m of the genes wrt ploidy of the segment
# remove genes with few observations
transcriptomics_data = transcriptomics_data %>% 
  group_by(hgnc_symbol) %>% 
  add_count() %>% 
  filter(n > 5)


# med = transcriptomics_data %>% 
#   group_by(hgnc_symbol, tot_cna) %>% 
#   summarise(expr_median = median(value))
# 
# 
# lm_by_gene_med = extract_lm_per_gene(med, formula = 'expr_median ~ tot_cna')
# good_genes_median = lm_by_gene_med %>% filter(!is.nan(p.value)) %>% filter(term == 'tot_cna') %>%  get_good_candidates_rsquared(n = 100) %>% 
#   pull(gene)
# 
# med %>% 
#   # filter(hgnc_symbol %in% good_genes_median) %>% 
#   filter(hgnc_symbol == 'SLC12A7') %>% 
#   # mutate(tot_cna = as.factor(tot_cna)) %>% 
#   ggplot(aes(tot_cna, expr_median, fill = tot_cna)) + 
#   geom_point() + 
#   geom_smooth(method = 'lm') + 
#   theme_bw() + 
#   facet_wrap(vars(hgnc_symbol), scales = 'free') + 
#   stat_poly_eq(use_label(c('eq', 'R2'))) 
# 
# plot_fit_statistics(lm_by_gene_med, 'r.squared') 
  

lm_by_gene = extract_lm_per_gene(transcriptomics_data, formula = 'value ~ tot_cna')

plot_fit_statistics(lm_by_gene, 'r.squared') + 
  ggtitle('RNA vs CNA R2 distribution')
ggsave('res/RNA_vs_CNA/RNA_vs_cna_r_squared.png', width = 5, height = 4)

plot_fit_statistics(lm_by_gene, 'estimate') + 
  ggtitle('Transcriptomics vs CNA linear coefficient distribution')
ggsave('res/RNA_vs_CNA/RNA_vs_cna_lm_coefficient.png', width = 5, height = 4)

plot_fit_statistics(lm_by_gene, 'p.value') + 
  ggtitle('Proteomics vs CNA p-values distribution') + 
  geom_vline(xintercept = 0.05, colour = 'darkred', linetype = 'dashed')
ggsave('res/RNA_vs_CNA/RNA_vs_cna_pvalues.png', width = 5, height = 4)

## get interesting genes based on their association type with the cna - positive or negative 
good_candidates_up_best = lm_by_gene %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  # dplyr::filter(estimate > mean(estimate)) %>% 
  get_good_candidates_estimate(., n = 100, which = 'up')

good_genes_up_rna = lm_by_gene %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::filter(estimate > 0)

good_candidates_down_best = lm_by_gene %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  get_good_candidates_estimate(., n = 100, which = 'down') %>% 
  dplyr::filter(estimate < 0)

good_genes_down_rna = lm_by_gene %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::filter(estimate < 0)

saveRDS(lm_by_gene, 'res/protein_vs_CNA/rna_vs_cna_fit.rds')
saveRDS(good_genes_down_rna, 'data/good_candidates_transcriptomics_down.rds')
saveRDS(good_genes_up_rna, 'data/good_candidates_transcriptomics_up.rds')
saveRDS(good_candidates_down_best, 'data/best_good_candidates_transcriptomics_down.rds')
saveRDS(good_candidates_up_best, 'data/best_good_candidates_transcriptomics_up.rds')

# take a look to how known cancer genes behave 
cgc_colon = readRDS('data/cancer_genes_somatic_colon_positions.rds') %>% 
  dplyr::rename(gene = hgnc_symbol)
cgc_somatic = readRDS('data/cancer_genes_somatic_positions.rds') %>% 
  dplyr::rename(gene = hgnc_symbol)
drivers_coad = readRDS('data/coad_drivers_positions.rds') %>% 
  dplyr::rename(gene = hgnc_symbol)

interesting_genes = list(
  'cgc_colon' = cgc_colon, 
  'cgc_somatic' = cgc_somatic, 
  'drivers_coad' = drivers_coad
)

interesting_genes_fit_not_filtered = lapply(interesting_genes, function(x) {
  lm_by_gene %>% 
    dplyr::filter(term == 'tot_cna') %>%
    dplyr::filter(gene %in% (unique(x$gene))) 
})

saveRDS(interesting_genes_fit_not_filtered, 'res/RNA_vs_CNA/interesting_genes_fit_not_filtered.rds')

interesting_genes_fit = lapply(interesting_genes, function(x) {
  lm_by_gene %>% 
    dplyr::filter(term == 'tot_cna') %>%
    dplyr::filter(gene %in% (unique(x$gene))) %>% 
    dplyr::filter(p.value < 0.05) 
  })

interesting_genes_fit = c(
  list('good_down' = good_genes_down_rna,
       'good_up' = good_genes_up_rna,
       'good_best_down' = good_candidates_down_best,
       'good_best_up' = good_candidates_up_best
       ), 
  interesting_genes_fit
)
saveRDS(interesting_genes_fit, 'res/RNA_vs_CNA/interesting_genes_fit.rds')

lapply(names(interesting_genes_fit), function(x) { 
  p = interesting_genes_fit[[x]] %>% 
    plot_fit_statistics(., 'p.value') + 
    ggtitle(x)
  ggsave(paste0('res/RNA_vs_CNA/interesting_genes/', x, '_pvalues.png'), width = 6, height = 5, plot = p)
})

lapply(names(interesting_genes_fit), function(x) { 
  p = interesting_genes_fit[[x]] %>% 
    plot_fit_statistics(., 'estimate') + 
    ggtitle(x)
  ggsave(paste0('res/RNA_vs_CNA/interesting_genes/', x, '_linear_coefficients.png'), width = 6, height = 5, plot = p)
})

# fits plotting
# colors = setNames(Polychrome::kelly.colors(length(transcriptomics_data$tot_cna %>% unique)), nm = transcriptomics_data$tot_cna %>% unique)

# plot the fit of good genes - RNA vs CNA


pp_fit_rna = lapply(interesting_genes_fit, function(x) {
  p = plot_fit(transcriptomics_data %>% filter(hgnc_symbol %in% x$gene), 
           x = 'tot_cna', 
           y = 'value', 
           facet = 'hgnc_symbol', 
           add_fit_line = TRUE) +
    xlab('Gene ploidy') + 
    ylab('RNA expression') 
  return(p)
})

lapply(names(pp_fit), function(x) {
  ggsave(paste0('res/RNA_vs_CNA/interesting_genes/', x, '_fit.png'),pp_fit[[x]], width = 10, height = 5)
})

p = plot_fit(transcriptomics_data %>% filter(hgnc_symbol %in% c('CHD4', 'NRAS', 'SMAD2')), 
             x = 'tot_cna', 
             y = 'value', 
             facet = 'hgnc_symbol', 
             add_fit_line = TRUE) +
  xlab('Gene ploidy') + 
  ylab('RNA expression') 

ggsave('res/RNA_vs_CNA/rna_vs_cna_specific_genes.png', p, width = 10, height = 6)
# saveRDS(good_candidates, 'data/good_fit_genes.rds')

# x = transcriptomics_data_low_CNAs %>% 
#   filter(hgnc_symbol == 'A1BG')
# 
# broom::tidy(lm(value ~ tot_cna, x))

# fit RNA vs CNA
transcriptomics_data %>% 
  ggplot(aes(x = tot_cna, y = value)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  stat_poly_eq(use_label(c('eq', 'R2')))+
  theme_bw() + 
  xlab('Gene ploidy') + 
  ylab('RNA expression')

ggsave('res/RNA_vs_CNA/RNA_vs_CNA_all_genes.png', width = 15, height = 10)

# checking TSG and oncogenes 
cancer_genes = read.table('data/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
  dplyr::as_tibble()

cancer_genes_somatic_colon = cancer_genes %>% 
  dplyr::filter(Germline != 'yes') %>% 
  dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  dplyr::pull(Gene.Symbol) %>% 
  unique

transcriptomics_data %>% 
  mutate(tot_cna = factor(tot_cna)) %>% 
  filter(hgnc_symbol %in% c(cancer_genes_somatic_colon, 'TP53')) %>% 
  ggplot(aes(x = tot_cna, y = value, fill = tot_cna)) + 
  geom_boxplot() + 
  geom_jitter() + 
  scale_fill_manual(values = colors) +
  theme_bw() + 
  facet_wrap(vars(hgnc_symbol), scales = 'free')

ggsave('res/cgs_important_genes_rna_vs_cna.png', width = 19, height = 10)

lm_by_gene %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(gene %in% c(cancer_genes_somatic_colon, 'TP53')) %>% 
  ggplot(aes(r.squared)) + 
  geom_histogram(binwidth = 0.01) + 
  theme_bw()
ggsave('res/r_squared_cgs_genes_rna_vs_protein.png')

# colors = setNames(ggthemes::calc_pal()(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)

# x = x %>% 
#   dplyr::mutate(karyotype = factor(karyotype, 
#                                    levels = c("1:0", "1:1", "2:0", "2:1", "3:0", "2:2", "3:1", "4:0", "3:2", "5:0", "3:3")))


# geom_smooth(method=lm , se=TRUE)


transcriptomics_data %>% 
  filter(hgnc_symbol == 'KRAS') %>%
  ggplot(aes(x = tot_cna, y = value)) + 
  # geom_violin() +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_poly_eq(use_label(c('eq', 'R2'))) + 
  # geom_boxplot(outliers = F) +
  # geom_jitter(alpha = 0.7, size = 1) + 
  # facet_wrap(~hgnc_symbol, scales = 'free') +
  # scale_fill_manual(values = colors) + 
  theme_bw() + 
  ggtitle('KRAS, expression vs cna')

saveRDS(transcriptomics_data_low_CNAs, 'data/pseudobulk_vs_cna.rds')

ggsave('res/pseudobulk_expression_selected_genes.png', width = 15, height = 10)

transcriptomics_data_low_CNAs %>% 
  mutate(expr_by_cna = value/tot_cna) %>% 
  ggplot(aes(x = tot_cna, y = expr_by_cna, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~region, scales = 'free') +
  xlab('tot_cna') + 
  ylab('expr/tot_cna') +
  scale_color_manual(values = colors) + 
  ggtitle('Normalized data, NAs not imputed')

transcriptomics_data_low_CNAs %>% 
  filter(region=='TP53') %>%
  filter(!is.na(karyotype)) %>% 
  ggplot(aes(value, fill = karyotype)) +
  geom_histogram() +
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  # mutate(expr_by_cna = value/tot_cna) %>% 
  # ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
  # geom_point() +
  theme_bw() +
  facet_wrap(~region, scales = 'free') +
  # scale_color_brewer(palette = 'Paired')
  # ggsci::scale_color_futurama()
  scale_fill_manual(values = colors) 