library(tidyverse)
library(ggpmisc)
source('organoids_analysis/R/proteomics/proteomics_utils.R')
source('organoids_analysis/R/functions_utils/fit_plots.R')
source('organoids_analysis/R/scRNA/cna_comparison_utils.R')

transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds') 

cancer_genes = read.table('data/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
  dplyr::as_tibble()

proteogenomics_data = proteogenomics_data %>% 
  group_by(proteomics_code, hgnc_symbol) %>% 
  mutate(mean_intensity = mean(norm_intensity)) %>% 
  dplyr::select(everything(), -c(replicate, norm_intensity)) %>% 
  distinct()

# merge rna and MS

rna_vs_prot = full_join(transcriptomics_data, proteogenomics_data, 
                        by = join_by('sample' == 'sample', 
                                     'hgnc_symbol' == 'hgnc_symbol', 
                                     'karyotype' == 'karyotype', 
                                     'Major' == 'Major', 
                                     'minor' == 'minor', 
                                     'tot_cna' == 'tot_cna', 
                                     'chr' == 'chr', 
                                     'mut_consequence' == 'mut_consequence', 
                                     'driver_label' == 'driver_label', 
                                     'is_mutated' == 'is_mutated')) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(!is.na(proteomics_code)) 
colors = setNames(Polychrome::kelly.colors(length(transcriptomics_data$tot_cna %>% unique)), nm = transcriptomics_data$tot_cna %>% unique)

rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol == 'KRAS') %>% 
  dplyr::mutate(tot_cna = as.character(tot_cna)) %>%
  ggplot(aes(x = value, 
             y = norm_intensity, 
             color = karyotype)) + 
  geom_point() +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~is_mutated) +
  theme_bw() + 
  xlab('RNAseq') + 
  ylab('Proteomics (mean of replicates)')

# look at all the results 
proteomics_genes_fit = readRDS('res/protein_vs_CNA/interesting_genes_fit.rds')
transcriptomics_genes_fit = readRDS('res/RNA_vs_CNA/interesting_genes_fit.rds')

# proteomics_genes_fit_v2 = proteomics_genes_fit[c('good_best_up', 'good_best_down', "cgc_colon", "cgc_somatic", "drivers_coad" )]
# transcriptomics_genes_fit_v2 = transcriptomics_genes_fit[c('good_best_up', 'good_best_down', "cgc_colon", "cgc_somatic", "drivers_coad" )]
# names(proteomics_genes_fit) = gsub('_prot', '', names(proteomics_genes_fit))
# names(transcriptomics_genes_fit) = gsub('_prot', '', names(transcriptomics_genes_fit))

genes_to_select = unique(c(proteogenomics_data$hgnc_symbol, transcriptomics_data$hgnc_symbol))

# by each of the two, then intersection

rna_vs_prot_Pgenes = lapply(proteomics_genes_fit_v2, function(x) {
  plot_all_fit(rna_vs_prot, genes = unique(x$gene), formula = 'mean_intensity ~ value')
})

rna_vs_prot_genes = lapply(transcriptomics_genes_fit_v2, function(x) {
  plot_all_fit(rna_vs_prot, genes = unique(x$gene), formula = 'mean_intensity ~ value')
})


good_prot = proteomics_genes_fit$good_up
good_rna = transcriptomics_genes_fit$good_up

common_up = intersect(good_prot$gene, good_rna$gene)

common_up_fit = plot_all_fit(rna_vs_prot, genes = common_up, formula = 'value ~ mean_intensity')
best_up_all = common_up_fit$fit %>% 
  filter(term == 'mean_intensity') %>% 
  filter(p.value < 0.05) %>% 
  filter(r.squared > 0.6) 

plot_fit(rna_vs_prot %>% 
           filter(hgnc_symbol %in% best_up_all$gene), 
         x = 'value', 
         y = 'mean_intensity', 
         facet = 'hgnc_symbol', 
         add_fit_line = TRUE) +
  xlab('RNA expression') + 
  ylab('Protein expression (mean of replicates)')

ggsave('res/RNA_vs_protein/good_res/best_fits_RNA_vs_prot.png', width = 25, height = 30)

rna_vs_prot_specific_genes = rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol %in% genes_to_select)

lm_by_genes_rna_prot = extract_lm_per_gene(rna_vs_prot_specific_genes, 'mean_intensity ~ value')

plot_fit_statistics(lm_by_genes_rna_prot, 'r.squared') + 
  ggtitle('RNA vs Proteins R2 distribution')
ggsave('res/RNA_vs_protein/RNA_vs_Prot_r_squared.png', width = 5, height = 4)

plot_fit_statistics(lm_by_genes_rna_prot, 'estimate') + 
  ggtitle('RNA vs Protein linear coefficient distribution')
ggsave('res/RNA_vs_protein/RNA_vs_prot_lm_coefficient.png', width = 5, height = 4)

plot_fit(rna_vs_prot_specific_genes, 
         x = 'value', 
         y = 'mean_intensity', 
         facet = 'hgnc_symbol', 
         add_fit_line = TRUE) +
  xlab('RNA expression') + 
  ylab('Protein expression (mean of replicates)')
ggsave('res/RNA_vs_protein/RNA_vs_prot_best_genes_fit.png', width = 30, height = 20)

# fit using the cgc coad genes
coad_genes = readRDS('data/cancer_genes_somatic_colon_positions.rds')
coad_genes = coad_genes$hgnc_symbol %>% unique

rna_vs_prot_cgc_coad = rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol %in% coad_genes)
  
lm_by_coad_genes_rna_prot = extract_lm_per_gene(rna_vs_prot_cgc_coad, 'mean_intensity ~ value')

plot_fit_statistics(lm_by_coad_genes_rna_prot, 'r.squared') + 
  ggtitle('COAD genes - RNA vs Proteins R2 distribution')
ggsave('res/RNA_vs_protein/RNA_vs_Prot_cgc_coad_genes_r_squared.png', width = 5, height = 4)

plot_fit_statistics(lm_by_coad_genes_rna_prot, 'estimate') + 
  ggtitle('COAD genes - RNA vs Protein linear coefficient distribution')
ggsave('res/RNA_vs_protein/RNA_vs_prot_cgc_coad_genes_lm_coefficient.png', width = 5, height = 4)

plot_fit(rna_vs_prot_cgc_coad, 
         x = 'value', 
         y = 'mean_intensity', 
         facet = 'hgnc_symbol', 
         add_fit_line = TRUE) +
  xlab('RNA expression') + 
  ylab('Protein expression (mean of replicates)')
ggsave('res/RNA_vs_protein/RNA_vs_prot_cgc_coad_genes_fit.png', width = 30, height = 20)


# compare fit results 

lm_by_gene_proteomics = readRDS('res/protein_vs_CNA/proteins_vs_cna_fit.rds')
lm_by_gene_transcriptomics = readRDS('res/RNA_vs_CNA/rna_vs_cna_fit.rds')

prot = lm_by_gene_proteomics %>% 
  filter(term == 'tot_cna') 
rna = lm_by_gene_transcriptomics %>% 
  filter(term == 'tot_cna') 

fits = full_join(prot, rna, by = 'gene', suffix = c('_prot', '_rna'))

fits = fits %>% 
  mutate(significance = case_when(
    (p.value_prot < 0.05 & p.value_rna < 0.05) ~ 'both', 
    (p.value_prot < 0.05 & (p.value_rna > 0.05 | is.na(p.value_rna))) ~ 'Proteomics', 
    ((p.value_prot > 0.05 | is.na(p.value_prot)) & p.value_rna < 0.05) ~ 'scRNA', 
    (p.value_prot > 0.05 & p.value_rna > 0.05) ~ 'None', 
    .default = 'None'
  )) %>% 
  mutate(significance = factor(significance, levels = c('None', 'Proteomics', 'scRNA', 'both')))

cols = c(
  'None' = 'gainsboro', 
  'Proteomics' = 'goldenrod', 
  'scRNA' = 'steelblue3', 
  'both' = 'forestgreen'
)

alphas = c(
  'None' = 0.3, 
  'Proteomics' = 1, 
  'scRNA' = 1, 
  'both' = 1
)

percentages = fits %>% 
  group_by(significance) %>% 
  summarise(proportion = (sum(n())/nrow(fits))*100)
vals = percentages %>% 
  mutate(val = paste0(significance, ' (', round(proportion, digits = 2), '%)')) %>% 
  pull(val)

fits %>% 
  ggplot(aes(estimate_rna, estimate_prot, color = significance, alpha = significance)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cols, labels = vals) + 
  scale_alpha_manual(values = alphas, labels = vals) +
  ylab('Proteomics linear coefficients') + 
  xlab('Transcriptomics linear coefficients') + 
  guides(color=guide_legend(title="Significant (pvalue < 0.05)"), 
         alpha=guide_legend(title="Significant (pvalue < 0.05)")) + 
  theme(legend.position = 'bottom') #+ 
  # facet_wrap(~significance, scales = 'free')
ggsave('res/RNA_vs_protein/good_res/comparison_linear_coeffcients.png', width = 8, height = 8)


# lm_by_genes_rna_prot %>% 
#   dplyr::filter(term == 'value') %>% 
#   ggplot(aes(estimate)) + 
#   geom_histogram(binwidth = 0.05) + 
#   theme_bw()
# 
# lm_by_genes_rna_prot %>% 
#   dplyr::filter(term == 'value') %>% 
#   ggplot(aes(r.squared)) + 
#   geom_histogram(binwidth = 0.05) + 
#   theme_bw()

# good_genes = lm_by_genes_rna_prot %>% 
#   dplyr::filter(term == 'value') %>% 
#   dplyr::filter(r.squared > 0.6) %>% 
#   pull(gene)

rna_vs_prot %>%
  dplyr::mutate(tot_cna = factor(tot_cna)) %>% 
  dplyr::filter(hgnc_symbol %in% genes_to_select) %>%
  ggplot(aes(x = value, 
             y = mean_intensity, 
             color = tot_cna)) + 
  geom_point() +
  # scale_color_manual(values = colors) +
  scale_color_brewer(palette = 'Set1') + 
  # facet_wrap(vars(region), scales = 'free') + 
  facet_wrap(vars(hgnc_symbol), scales = 'free') +
  # geom_smooth(method = 'lm', formula = 'y ~ x') + 
  # facet_wrap(vars(tot_cna)) + 
  xlab('RNA') + 
  ylab('Protein (mean of replicates)') + 
  theme_bw()

ggsave('res/proteins_vs_rna_good_genes_in_RNA_fitted.png', width = 20, height = 15)

ggsave('res/proteins_vs_rna_selected_genes.png', width = 10, height = 10)


# checking TSG and oncogenes 
cancer_genes = read.table('data/Cancer_gene_census_tier1.csv', sep = ",", header = T) %>% 
  dplyr::as_tibble()

cancer_genes_somatic_colon = cancer_genes %>% 
  dplyr::filter(Germline != 'yes') %>% 
  dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  dplyr::pull(Gene.Symbol) %>% 
  unique

cancer_genes_somatic_colon = c(cancer_genes_somatic_colon, 'TP53')

rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol %in% cancer_genes_somatic_colon) %>%
  dplyr::mutate(tot_cna = as.factor(tot_cna)) %>% 
  # dplyr::filter(tot_cna <= 6) %>% 
  ggplot(aes(x = value, 
             y = mean_intensity, 
             color = tot_cna)) + 
  geom_point() +
  scale_color_brewer(palette = 'Set1') +
  # facet_wrap(vars(region), scales = 'free') + 
  facet_wrap(vars(hgnc_symbol), scales = 'free') +
  # facet_wrap(vars(tot_cna)) + 
  xlab('RNA') + 
  ylab('Protein (mean of replicates)') + 
  theme_bw() 

ggsave('res/cgs_rna_vs_protein.png', width = 12, height = 10)

rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol %in% genes_to_select) %>%
  # dplyr::mutate(tot_cna = as.factor(tot_cna)) %>% 
  # dplyr::filter(tot_cna <= 6) %>% 
  ggplot(aes(x = tot_cna, 
             y = mean_intensity)) + 
  geom_point() +
  # facet_wrap(vars(region), scales = 'free') + 
  facet_wrap(vars(hgnc_symbol), scales = 'free') +
  # facet_wrap(vars(tot_cna)) + 
  xlab('Ploidy') + 
  ylab('Protein (mean of replicates)') + 
  theme_bw() + 
  geom_smooth(method = 'lm') + 
  ggtitle('Protein expression vs gene ploidy')

ggsave('res/cgs_protein_vs_cna.png', width = 13, height = 8)

rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol %in% cancer_genes_somatic_colon) %>%
  # dplyr::mutate(tot_cna = as.factor(tot_cna)) %>% 
  # dplyr::filter(tot_cna <= 6) %>% 
  ggplot(aes(x = tot_cna, 
             y = value)) + 
  geom_point() +
  # facet_wrap(vars(region), scales = 'free') + 
  facet_wrap(vars(hgnc_symbol), scales = 'free') +
  # facet_wrap(vars(tot_cna)) + 
  xlab('Ploidy') + 
  ylab('RNA expression') + 
  theme_bw() + 
  geom_smooth(method = 'lm') + 
  ggtitle('RNA expression vs gene ploidy')

ggsave('res/cgs_rna_vs_cna.png', width = 13, height = 8)
  

# all genes from CGC
cancer_genes_somatic = cancer_genes %>% 
  dplyr::filter(Somatic == 'yes') %>% 
  # dplyr::filter(grepl('*colorectal*', Tumour.Types.Somatic.)) %>% 
  dplyr::pull(Gene.Symbol) %>% 
  unique