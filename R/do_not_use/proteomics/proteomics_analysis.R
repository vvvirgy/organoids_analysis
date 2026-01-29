# new normalization method
library(tidyverse)
library(ggpmisc)
source('organoids_analysis/R/proteomics/proteomics_utils.R')
source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')


proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')

proteogenomics_data_v2 = proteogenomics_data %>%
  group_by(proteomics_code, hgnc_symbol) %>%
  mutate(mean_intensity = mean(norm_intensity)) %>%
  dplyr::select(everything(), -c(replicate, norm_intensity)) %>%
  distinct()

lm_by_gene_proteomics = extract_lm_per_gene(proteogenomics_data_v2, formula = 'mean_intensity ~ tot_cna')
saveRDS(lm_by_gene_proteomics, 'res/protein_vs_CNA/proteins_vs_cna_fit.rds')

plot_fit_statistics(lm_by_gene_proteomics, 'estimate') + 
  ggtitle('Proteomics vs CNA linear coefficient distribution')
ggsave('res/protein_vs_CNA/Protein_vs_cna_lm_coefficient.png', width = 5, height = 4)

plot_fit_statistics(lm_by_gene_proteomics, 'p.value') + 
  ggtitle('Proteomics vs CNA p-values distribution') + 
  geom_vline(xintercept = 0.05, colour = 'darkred', linetype = 'dashed')
ggsave('res/protein_vs_CNA/Protein_vs_cna_pvalues.png', width = 5, height = 4)

# get interesting genes based on their association type with the cna - positive or negative 
good_candidates_up_best = lm_by_gene_proteomics %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  # dplyr::filter(estimate > mean(estimate)) %>% 
  get_good_candidates_estimate(., n = 100, which = 'up')

good_genes_up_prot = lm_by_gene_proteomics %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::filter(estimate > 0)

good_candidates_down_best = lm_by_gene_proteomics %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  get_good_candidates_estimate(., n = 100, which = 'down') %>% 
  dplyr::filter(estimate < 0)

good_genes_down_prot = lm_by_gene_proteomics %>% 
  dplyr::filter(term == 'tot_cna') %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::filter(estimate < 0)

saveRDS(lm_by_gene_proteomics, 'res/protein_vs_CNA/protein_vs_cna_fit.rds')
saveRDS(good_genes_down_prot, 'data/good_candidates_proteomics_down.rds')
saveRDS(good_genes_up_prot, 'data/good_candidates_proteomics_up.rds')
saveRDS(good_candidates_down_best, 'data/best_good_candidates_proteomics_down.rds')
saveRDS(good_candidates_up_best, 'data/best_good_candidates_proteomics_up.rds')

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
  lm_by_gene_proteomics %>% 
    dplyr::filter(term == 'tot_cna') %>%
    dplyr::filter(gene %in% (unique(x$gene))) 
})

saveRDS(interesting_genes_fit_not_filtered, 'res/protein_vs_CNA/interesting_genes_fit_not_filtered.rds')

interesting_genes_fit = lapply(interesting_genes, function(x) {
  lm_by_gene_proteomics %>% 
    dplyr::filter(term == 'tot_cna') %>%
    dplyr::filter(gene %in% (unique(x$gene))) %>% 
    dplyr::filter(p.value < 0.05) 
})

interesting_genes_fit = c(
  list('good_up' = good_genes_up_prot, 
       'good_down' = good_genes_down_prot, 
       'good_best_up' = good_candidates_up_best, 
       'good_best_down' = good_candidates_down_best), 
  interesting_genes_fit
)
saveRDS(interesting_genes_fit, 'res/protein_vs_CNA/interesting_genes_fit.rds')

lapply(names(interesting_genes_fit), function(x) { 
  p = interesting_genes_fit[[x]] %>% 
    plot_fit_statistics(., 'p.value') + 
    ggtitle(x)
  ggsave(paste0('res/protein_vs_CNA/interesting_genes/', x, '_pvalues.png'), width = 6, height = 5, plot = p)
})

lapply(names(interesting_genes_fit), function(x) { 
  p = interesting_genes_fit[[x]] %>% 
    plot_fit_statistics(., 'estimate') + 
    ggtitle(x)
  ggsave(paste0('res/protein_vs_CNA/interesting_genes/', x, '_linear_coefficients.png'), width = 6, height = 5, plot = p)
})

# fits plotting
# colors = setNames(Polychrome::kelly.colors(length(transcriptomics_data$tot_cna %>% unique)), nm = transcriptomics_data$tot_cna %>% unique)

# plot the fit of good genes - RNA vs CNA
interesting_genes_fit_v2 = interesting_genes_fit[-c(1,2)]

pp_fit = lapply(interesting_genes_fit_v2, function(x) {
  p = plot_fit(proteogenomics_data_v2 %>% filter(hgnc_symbol %in% x$gene), 
               x = 'tot_cna', 
               y = 'mean_intensity', 
               facet = 'hgnc_symbol', 
               add_fit_line = TRUE) +
    xlab('Gene ploidy') + 
    ylab('Protein expression (mean of replicates)') 
  return(p)
})

lapply(names(pp_fit), function(x) {
  ggsave(paste0('res/protein_vs_CNA/interesting_genes/', x, '_fit.png'),pp_fit[[x]], width = 10, height = 6)
})



# later, move in a different script
# compare with transcriptomics

ggvenn::ggvenn(list(Proteomics = good_candidates_prot$gene, Transcritomics = good_candidates$gene), auto_scale = T)
ggsave('res/venn_best_fitted_genes_prot_transcritomics.png', bg = 'white', width = 6, height = 7)

# fits plotting
colors = setNames(Polychrome::kelly.colors(length(proteogenomics_data_v2$tot_cna %>% unique)), nm = proteogenomics_data_v2$tot_cna %>% unique)

# plot the fit of good genes - RNA vs CNA
best_genes = plot_fit(proteogenomics_data_v2 %>% 
                        filter(hgnc_symbol %in% good_candidates_prot$gene), 
                      x = 'tot_cna', 
                      y = 'mean_intensity', 
                      facet = 'hgnc_symbol', 
                      add_fit_line = TRUE) +
  xlab('Gene ploidy') + 
  ylab('Protein expression (mean of replicates)')

ggsave('res/protein_vs_CNA/best_genes_100_fit_cna_vs_protein.png', best_genes, width = 30, height = 20)

# plot the fit of the genes associated with coad
coad_genes = readRDS('data/cancer_genes_somatic_colon_positions.rds')
coad_genes = coad_genes$hgnc_symbol %>% unique

# see their fits
lm_by_gene_proteomics %>% 
  filter(gene %in% coad_genes) %>% 
  plot_fit_statistics(., 'estimate') +
  ggtitle('coad CGC genes - Protein vs CNA linear coefficient distribution')

ggsave('res/protein_vs_CNA/fit_coad_genes_protein_vs_cna_linear_coefficient.png', width = 5, height = 4)

lm_by_gene_proteomics %>% 
  filter(gene %in% coad_genes) %>% 
  plot_fit_statistics(., 'r.squared') +
  ggtitle('coad CGC genes - Protein vs CNA R2 distribution')

ggsave('res/protein_vs_CNA/fit_coad_genes_protein_vs_cna_R2.png', width = 5, height = 4)

coad_genes_plot = plot_fit(proteogenomics_data_v2 %>% filter(hgnc_symbol %in% coad_genes), 
                           x = 'tot_cna', 
                           y = 'mean_intensity', 
                           facet = 'hgnc_symbol', 
                           add_fit_line = TRUE) + 
  xlab('Gene ploidy') + 
  ylab('Protein expression (mean of replicates)')

ggsave('res/protein_vs_CNA/coad_genes_cgc_fit_cna_vs_proteins.png', coad_genes_plot, width = 30, height = 20)


proteogenomics_data %>% 
  filter(hgnc_symbol == '') %>% 
  ggplot(aes(x = tot_cna, y = mean_norm_intensity)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~hgnc_symbol, scales = 'free') +
  xlab('total cnas') + 
  scale_color_manual(values = colors) +
  ggtitle('Normalized data, NAs not imputed')


# ----------------------------------------------------------------------------------------------------------------------------------------

# many genes are widely dispersed

proteogenomics_data = proteogenomics_data %>% 
  group_by(proteomics_code, hgnc_symbol) %>% 
  mutate(mean_intensity = mean(norm_intensity)) %>% 
  dplyr::select(everything(), -c(replicate, norm_intensity)) %>% 
  distinct()

sd_genes_proteomics = proteogenomics_data %>% 
  group_by(hgnc_symbol, tot_cna) %>% 
  summarise(sd_intensities = sd(mean_intensity))

sd_genes_proteomics %>% 
  filter(tot_cna < 6) %>% 
  ggplot(aes(sd_intensities)) + 
  geom_histogram(binwidth = 0.01) + 
  facet_wrap(~tot_cna, scales = 'free')

proteogenomics_data_v2 %>% 
  filter(hgnc_symbol == 'KRAS') %>% 
  # mutate(tot_cna = factor(tot_cna)) %>% 
  ggplot(aes(y = mean_intensity, x = tot_cna, label = sample)) +
  geom_point() + 
  theme_bw() + 
  ggrepel::geom_label_repel() +
  geom_smooth(method = 'lm') +
  ggtitle('KRAS, DNA vs Protein')

# ----------------------------------------------------------------------------------------------------------------------------------------
# PLOTS

# plot a bit of statistics over ploidy
genes_cna_status %>% 
  ggplot(aes(CNt, fill = karyotype)) +
  geom_bar(stat = 'count', position = 'dodge') + 
  theme_bw()

# testing a thing --> multiple CNA events on the same gene
drivers_to_check_correct_samples %>% 
  filter(sample == 17, region == 'SRC') %>% 
  dplyr::select(chr, from, to, Major, minor) %>% 
  tidyr::pivot_longer(cols = c(Major, minor), names_to = 'Allele', values_to = 'value') %>% 
  ggplot(aes(
    x = from, 
    xend = to, 
    y = value, 
    yend = value, 
    color = Allele
  )) + 
  geom_segment() + 
  theme_bw() + 
  ggtitle('sample 17, gene SRC, multiple CNAs')

ggsave('res/multiple_cnas.png')



proteomic_genes %>%
  filter(Genes == "TP53") %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(vars(Genes))



colors = setNames(Polychrome::kelly.colors(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)
# colors = setNames(ggthemes::calc_pal()(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)

# x = x %>% 
#   dplyr::mutate(karyotype = factor(karyotype, 
#                                    levels = c("1:0", "1:1", "2:0", "2:1", "3:0", "2:2", "3:1", "4:0", "3:2", "5:0", "3:3")))

x %>% 
  filter(region %in% genes_to_check) %>% 
  ggplot(aes(x = tot_cna, y = value, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~region, scales = 'free') +
  xlab('total cnas') + 
  scale_color_manual(values = colors) +
  ggtitle('Normalized data, NAs not imputed')
# geom_smooth(method=lm , se=TRUE)

x %>% 
  mutate(expr_by_cna = value/tot_cna) %>% 
  ggplot(aes(x = tot_cna, y = expr_by_cna, colour = karyotype)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~region, scales = 'free') +
  xlab('tot_cna') + 
  ylab('expr/tot_cna') +
  scale_color_manual(values = colors) + 
  ggtitle('Normalized data, NAs not imputed')

ggsave('res/cnas_vs_expression_allele.png', width = 15, height = 15)

# x %>% 
#   filter(gene %in% genes_to_check) %>% 
#   # full_join(., cancer_genes, by = join_by('gene' == 'Gene.Symbol')) %>% 
#   # dplyr::rename(Role = Role.in.Cancer) %>% 
#   # filter(driver %in% c("APC", "KRAS", "TP53", "PIK3CA")) %>% 
#   # filter(gene=='KRAS') %>%
#   # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
#   mutate(expr_by_cna = value/tot_cna) %>%
#   ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   facet_wrap(~gene, scales = 'free') +
#   xlab('protein_expr/tot_cna') + 
#   scale_color_manual(values = colors) + 
#   ggtitle('Normalized data, NAs not imputed')


x %>% 
  # filter(gene=='KRAS') %>%
  filter(!is.na(karyotype)) %>% 
  ggplot(aes(value, fill = karyotype)) +
  geom_histogram(binwidth = 0.1) +
  # filter(karyotype %in% c('1:1', '1:0', '2:0', '2:1', '3:0', '2:2', '3:1')) %>%
  # mutate(expr_by_cna = value/tot_cna) %>% 
  # ggplot(aes(x = expr_by_cna, y = tot_cna, colour = karyotype)) +
  # geom_point() +
  theme_bw() +
  facet_wrap(~region, scales = 'free') +
  # scale_color_brewer(palette = 'Paired')
  # ggsci::scale_color_futurama()
  scale_fill_manual(values = colors) + 
  ggtitle('Normalized data, NAs imputed')



