.libPaths()
rm(list=ls())
library(glmnet)
library(tidyverse)
library(ComplexHeatmap)
library(ggpmisc)
library("factoextra")
library(ggExtra)
library(patchwork)

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/comparing_models.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/coefficients_analysis.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/glm_rmse.R')
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/ht_utils.R')
create_annotation = function(x, # data 
                             ann_colors, 
                             position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position, 
                                    annotation_legend_param = list(nrow = 3, width = 12, by_row = T))
  
}

# proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
# transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')
# 
# all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))
# 
# rna_vs_prot = inner_join(transcriptomics_data, proteogenomics_data, 
#                         by = join_by('sample' == 'sample', 
#                                      'hgnc_symbol' == 'hgnc_symbol', 
#                                      'karyotype' == 'karyotype', 
#                                      'Major' == 'Major', 
#                                      'minor' == 'minor', 
#                                      'tot_cna' == 'tot_cna', 
#                                      'chr' == 'chr', 
#                                      'mut_consequence' == 'mut_consequence', 
#                                      'driver_label' == 'driver_label', 
#                                      'is_mutated' == 'is_mutated', 
#                                      'is_driver_intogen' == 'is_driver_intogen', 
#                                      'CGC_role_COAD' == 'CGC_role_COAD', 
#                                      'CGC_role_PANCANCER' == 'CGC_role_PANCANCER')) %>% 
#   dplyr::filter(!is.na(value)) %>% 
#   dplyr::filter(!is.na(proteomics_code)) %>% 
#   dplyr::filter(hgnc_symbol %in% all_genes) %>% 
#   dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
#   dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
#   dplyr::rename(protein_expression = mean_intensity) %>% 
#   dplyr::rename(rna_expression = value) %>% 
#   dplyr::filter(!is.na(protein_expression)) %>% 
#   dplyr::filter(!is.na(rna_expression))

rna_vs_prot = readRDS('data/rna_vs_prot_data.rds')

# fit_multiresponse_model = readRDS('data/multiresponse_glm_lasso_v2.rds')
fit_multiresponse_model = readRDS('data/multiresponse_glm_lasso_multiplicity_v1.rds')

# fit_multiresponse_model = Filter(Negate(is.null), fit_multiresponse_model)
# some genes have not been fitted since they do not have enough observations
fit_multiresponse_model = Filter(function(x) {length(x) > 1}, fit_multiresponse_model) 
names(fit_multiresponse_model) = lapply(fit_multiresponse_model, function(x) {x$gene}) %>% unlist

coefficients_all = lapply(fit_multiresponse_model, get_coefficients_multiresponse)
coefficients_all = lapply(names(coefficients_all), function(x) {
  coefficients_all[[x]] %>% 
    dplyr::mutate(gene = x)
}) %>% 
  bind_rows() %>% 
  # relocate(gene, .before = protein_expression_tot_cna) %>% 
  tibble::column_to_rownames('gene')
saveRDS(coefficients_all, 'data/coefficients_all_multiresponse.rds')
# remove rows that have only zeros for simplicity of plotting
# coefficients_all = coefficients_all[!apply(coefficients_all == 0, 1, all), ]

# plot coefficients based on the effect 

# genes with a role in CGC for coad --> create a function
rna_vs_prot = rna_vs_prot %>% 
  dplyr::mutate(
    CGC_role_COAD = case_when(
      CGC_role_COAD == 'oncogene, fusion' ~ 'oncogene', 
      CGC_role_COAD == 'TSG, fusion' ~ 'TSG', 
      CGC_role_COAD == 'oncogene, TSG, fusion' ~ 'oncogene, TSG', 
      CGC_role_COAD == 'None' ~ NA,
      .default = CGC_role_COAD
    ),
    CGC_role_PANCANCER = case_when(
      CGC_role_PANCANCER == 'oncogene, fusion' ~ 'oncogene', 
      CGC_role_PANCANCER == 'TSG, fusion' ~ 'TSG', 
      CGC_role_PANCANCER == 'oncogene, TSG, fusion' ~ 'oncogene, TSG', 
      CGC_role_PANCANCER == 'None' ~ NA,
      .default = CGC_role_PANCANCER
    )
  )


groups_interesting = c('CGC_role_PANCANCER', 'is_driver_intogen', 'CGC_role_COAD')
genes_by_groups = lapply(groups_interesting, function(x) {
  get_genes_by_group(rna_vs_prot, x)
})
names(genes_by_groups) = groups_interesting

coefficients_groups = lapply(genes_by_groups, function(x) {
  get_groups_coefficients(x, coefficients_all)
})
coefficients_groups$is_driver_intogen = coefficients_groups$is_driver_intogen[-1]

# create annotations -- always the same!
maf = readRDS('data/cohort_maf_drivers.rds')
top_genes = coefficients_all %>% 
  tibble::rownames_to_column('gene') %>% 
  dplyr::filter(gene %in% (maf@data$VEP.SYMBOL %>% unique)) %>% 
  tibble::column_to_rownames('gene') %>% 
  as.matrix()

ann_data = tibble(
  ann = colnames(coefficients_all)
) %>% 
  dplyr::mutate(
    assay = str_extract(ann, 'protein|rna')
  ) %>% 
  dplyr::mutate(predictor = str_extract(ann, 'n_wt|n_low|n_alt|n_trunc')) %>% 
  dplyr::select(-ann)

ann_col = list(
  assay = setNames(
    nm = c('protein', 'rna'), 
    c('orange2', 'lightseagreen')
  ), 
  predictor = setNames(
    nm = c('n_wt', 'n_low', 'n_alt', 'n_trunc'), 
    viridis::viridis(n = 4, direction = -1)
    # nm = c('cna', 'mutation'), 
    # c('mediumvioletred', 'royalblue4', '')
  )
) 

ann = create_annotation(ann_data, 
                  ann_colors = ann_col,
                  position = 'column')

all_heatmaps = lapply(coefficients_groups, function(x) {
  plot_heatmap_coefficients(x, ann)
})

# plot_heatmap_coefficients(list('top_drivers_cohort' = top_genes), ann)

pdf('res/fit_glm_multiresponse/coefficients_heatmap_cgc_coad_genes.pdf')
all_heatmaps$CGC_role_COAD
graphics.off()

pdf('res/fit_glm_multiresponse/coefficients_heatmap_cgc_pancancer_genes.pdf')
all_heatmaps$CGC_role_PANCANCER
graphics.off()

pdf('res/fit_glm_multiresponse/coefficients_heatmap_intogen_driver_genes.pdf')
all_heatmaps$is_driver_intogen
graphics.off()

ht_top_genes = plot_heatmap_coefficients(list('top_mutated_genes' = top_genes), ann)
pdf('res/fit_glm_multiresponse/coefficients_heatmap_top_mutated_driver_genes.pdf', height = 9)
ht_top_genes
graphics.off()

ht = ComplexHeatmap::Heatmap(as.matrix(coefficients_all), 
             col = colors,
             bottom_annotation = ann, 
             name = 'Predictor coefficients', 
             show_column_names = F, 
             column_title = 'all genes fitted', 
             heatmap_legend_param = list(direction = "horizontal"), 
             row_names_gp = gpar(fontsize = 2)
             # show_row_names = F
)

pdf('res/heatmap_coefficients_all_genes.pdf', width = 30, height = 50)
draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


estimate_err = compute_errors(fit_multiresponse_model, 
                              n_obs = 4,
                              data = rna_vs_prot, 
                              response = c('protein_expression', 'rna_expression'),
                              model = as.formula('~ n_low + n_alt + n_trunc + n_wt'), 
                              lambda = "lambda.min")
estimate_err = Filter(Negate(is.null), estimate_err)
names(estimate_err) = lapply(estimate_err, function(x) {x$gene}) %>% unlist

estimate_err_all = bind_rows(estimate_err)

r2 = lapply(fit_multiresponse_model, get_r2) 
names(r2) = names(fit_multiresponse_model)

r2_all = tibble(
  gene = names(r2), 
  r2 = unlist(r2)
)

estimate_err_all %>% 
  ggplot(aes(rna_expression_R2)) + 
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  ggtitle('RNA fit R2')
ggsave('res/fit_glm_multiresponse/R2_rna.png')

estimate_err_all %>% 
  ggplot(aes(protein_expression_R2)) + 
  geom_histogram(binwidth = 0.01) + 
  theme_bw() + 
  ggtitle('Protein fit R2')
ggsave('res/fit_glm_multiresponse/R2_protein.png')
 
saveRDS(estimate_err_all, 'data/errors_estimate_all_multiresponse_fit_v1.rds')

# visualizing the coefficients
couples = c('n_low', 'n_wt', 'n_alt', 'n_trunc')

coefficients_scatters = lapply(couples, function(nn) {
  df = coefficients_all %>% 
    dplyr::select(ends_with(nn)) 
  
  colnames(df) = gsub(paste0('_', nn), '', colnames(df))
  df = df %>% 
    dplyr::filter(protein_expression != 0, rna_expression != 0)
  
  p = df %>% 
    ggplot2::ggplot(aes(x = protein_expression, y = rna_expression)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste0(nn, ' coefficients RNA vs Protein')) + 
    xlab('Protein coefficients') + 
    ylab('RNA coefficients')
  
  p = ggExtra::ggMarginal(p, type="histogram")
  
  ggsave(paste0('res/', nn, '_coefficients.png'), p, width = 8, height = 8)
  
  return(p)
    
})
wrap_plots(coefficients_scatters)
ggsave('res/all_coefficients.png', width = 10, height = 10)
# dbscan::kNNdistplot() 
  
combs_rna = combn(colnames(coefficients_all) %>% grep('rna', ., value = T), 2)
combs_prot = combn(colnames(coefficients_all) %>% grep('protein', ., value = T), 2)

counts_alleles = rna_vs_prot %>% 
  dplyr::group_by(hgnc_symbol) %>% 
  dplyr::summarise(low = sum(n_low), alt = sum(n_alt), trunc = sum(n_trunc), wt = sum(n_wt))


rna_combinations = apply(combs_rna, 2, function(x) {
  print(x)
  
  el = strsplit(x, split = '_') %>% 
    lapply(., function(s) {
      s %>% last
    }) %>% unlist
  
  genes = counts_alleles %>% 
    dplyr::select(hgnc_symbol, any_of(el)) %>% 
    filter(across(all_of(el), ~ .x != 0)) %>% 
    pull(hgnc_symbol)
  
  
  coefficients_all[genes, ] %>%
    dplyr::select(any_of(x)) %>% 
    ggplot(aes(.data[[x[1]]], .data[[x[2]]])) + 
    geom_point() +
    theme_bw()
})
wrap_plots(rna_combinations)

prot_combinations = apply(combs_prot, 2, function(x) {
  
  el = strsplit(x, split = '_') %>% 
    lapply(., function(s) {
      s %>% last
    }) %>% unlist
  
  genes = counts_alleles %>% 
    dplyr::select(hgnc_symbol, any_of(el)) %>% 
    filter(across(all_of(el), ~ .x != 0)) %>% 
    pull(hgnc_symbol)
  
  
  coefficients_all[genes, ] %>%
    dplyr::select(any_of(x)) %>% 
    ggplot(aes(.data[[x[1]]], .data[[x[2]]])) + 
    geom_point() +
    theme_bw()
})
wrap_plots(prot_combinations) / wrap_plots(rna_combinations)
ggsave('data/comparison_betas_all_fits.png', width = 15, height = 15)

(prot_combinations[[3]] + rna_combinations[[3]]) / (prot_combinations[[5]] + rna_combinations[[5]])
ggsave('res/fit_glm_multiresponse/most_interesting_coefficients.png', width = 10, height = 10)
ggsave('res/fit_glm_multiresponse/wt_low_coefficients.png')

prot_combinations[[5]] + rna_combinations[[5]]
ggsave('res/fit_glm_multiresponse/wt_low_coefficients.png')


apply(coefficients_all, 2, max)
apply(coefficients_all, 2, min)

coefficients_all[which.max(coefficients_all$protein_expression_n_alt), ]

smad2_rna = rna_vs_prot %>% 
  filter(hgnc_symbol == 'SMAD2') %>% 
  dplyr::mutate(effect = paste(n_wt, n_low, n_alt, n_trunc, sep = '-')) %>% 
  ggplot(aes(y = rna_expression, fill = effect, x = category)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_fill_viridis_d('viridis') + 
  ylab('RNA expression')

smad2_prot = rna_vs_prot %>% 
  filter(hgnc_symbol == 'SMAD2') %>% 
  dplyr::mutate(effect = paste(n_wt, n_low, n_alt, n_trunc, sep = '-')) %>% 
  ggplot(aes(y = protein_expression, fill = effect, x = category)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_fill_viridis_d('viridis') + 
  ylab('Protein expression')

(smad2_rna +smad2_prot) +
  plot_annotation(title = 'SMAD2 expression') + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave('res/smad2_expression.png', width = 10, height = 5)

rna_vs_prot %>% 
  filter(hgnc_symbol == 'MSH6') %>% 
  dplyr::mutate(effect = paste(n_wt, n_low, n_alt, n_trunc, sep = '-')) %>% 
  ggplot(aes(y = rna_expression, fill = effect, x = category)) + 
  geom_boxplot() + 
  theme_bw()


coefficients_all['SMAD2', ] %>% 
  reshape2::melt() %>% 
  mutate(assay = case_when(grepl('protein', variable) ~ 'Protein', .default = 'RNA')) %>% 
  mutate(variable = gsub('protein_expression_|rna_expression_', '', variable)) %>% 
  ggplot(aes(x = variable, y = value, fill = assay)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() +
  scale_fill_manual(values = setNames(nm = c('Protein', 'RNA'), 
                    c('orange2', 'lightseagreen'))) + 
  xlab('Terms') + 
  ylab('Coefficients') + 
  ggtitle('SMAD2 coefficients')
ggsave('res/SMAD2_coefficients.png')

coefficients_all %>% 
  ggplot(aes(rna_expression_n_wt, rna_expression_n_low)) + 
  geom_point() + 
  theme_bw()

colnames(df) = gsub(paste0('_', nn), '', colnames(df))
df = df %>% 
  dplyr::filter(protein_expression != 0, rna_expression != 0)

p = df %>% 
  ggplot2::ggplot(aes(x = protein_expression, y = rna_expression)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(paste0(nn, ' coefficients RNA vs Protein')) + 
  xlab('Protein') + 
  ylab('RNA')

p = ggExtra::ggMarginal(p, type="histogram")



rna_vs_prot %>% 
  dplyr::group_by(hgnc_symbol) %>% 
  summarise(mean_prot = mean(protein_expression), mean_rna = mean(rna_expression)) %>% 
  ggplot(aes(mean_prot, mean_rna)) + 
  geom_point()

rna_vs_prot %>% 
  # dplyr::group_by(hgnc_symbol) %>%
  # summarise(mean_prot = mean(protein_expression), mean_rna = mean(rna_expression)) %>%
  ggplot(aes(protein_expression, rna_expression)) +
  geom_point() + 
  theme_bw()


coefficients_all %>% 
  dplyr::filter(protein_expression_tot_cna != 0, rna_expression_tot_cna != 0) %>%
  ggplot(aes(y = protein_expression_tot_cna, x = rna_expression_tot_cna)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(method = 'lm') + 
  stat_poly_eq(use_label('eq')) +  
  xlab('RNA cna coefficients') + 
  ylab('Protein cna coefficients')
ggsave('res/glm_multiresponse_coefficients_cna.png')

coefficients_all %>% 
  dplyr::filter(protein_expression_mutation_statusMutated != 0, rna_expression_mutation_statusMutated != 0) %>%
  ggplot(aes(y = protein_expression_mutation_statusMutated, x = rna_expression_mutation_statusMutated)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(method = 'lm') + 
  stat_poly_eq(use_label('eq')) +  
  xlab('RNA mut coefficients') + 
  ylab('Protein mut coefficients')
ggsave('res/glm_multiresponse_coefficients_mut.png')

# trying to define some groups based on the coefficients --> classifying
coefficients_all %>% 
  mutate(prot_effect = case_when(protein_expression_tot_cna > protein_expression_mutation_statusMutated ~ 'cna_protein_stronger', 
                   .default = 'mut_protein_stronger')) %>% 
  pull(prot_effect) %>% 
  table

coefficients_all %>% 
  mutate(prot_effect = case_when(protein_expression_tot_cna > rna_expression_tot_cna ~ 'cna_protein_stronger', 
                                 .default = 'cna_rna_stronger')) %>% 
  pull(prot_effect) %>% 
  table

coefficients_all = coefficients_all %>%
  rename(prot_cna = protein_expression_tot_cna) %>%
  rename(prot_mut = protein_expression_mutation_statusMutated) %>%
  rename(rna_cna = rna_expression_tot_cna) %>%
  rename(rna_mut = rna_expression_mutation_statusMutated)

classification_genes = coefficients_all %>% 
  dplyr::mutate(
    stronger_rna = 
      case_when(
        rna_cna > rna_mut ~ 'CNA', 
        rna_mut > rna_cna ~ 'mut', 
        (rna_mut == 0 & rna_cna == 0) ~ 'no effect'
      ), 
    stronger_prot =
      case_when(
        prot_cna > prot_mut ~ 'CNA', 
        prot_mut > prot_cna ~ 'mut', 
        (prot_cna == 0 & prot_mut == 0) ~ 'no effect'
      ), 
    stronger_cna = 
      case_when(
        prot_cna > rna_cna ~ 'Prot', 
        prot_cna < rna_cna ~ 'RNA', 
        (prot_cna == 0 & rna_cna == 0) ~ 'no effect'
      ), 
    stronger_mut = 
      case_when(
        prot_mut > rna_mut ~ 'Prot', 
        prot_mut < rna_mut ~ 'RNA', 
        (prot_mut == 0 & rna_mut == 0) ~ 'no effect'
      )
  )

genes_by_classes = classification_genes %>% 
  dplyr::mutate(
    effect = 
      case_when(
        stronger_rna == stronger_prot ~ stronger_rna, 
        stronger_rna != stronger_prot ~ 'Different effect'
      )
  ) %>% 
  tibble::rownames_to_column('gene') %>% 
  dplyr::group_by(effect) %>% 
  group_split()
names(genes_by_classes) = lapply(genes_by_classes, function(x) {x$effect %>% unique}) %>% unlist
genes_by_classes = lapply(genes_by_classes, function(x) {
  x$gene %>% unique
})

saveRDS(genes_by_classes, 'data/genes_by_classes.rds')
# 
# coefficients_all %>% 
#   dplyr::mutate(cna_diff = rna_cna - prot_cna) %>% 
#   dplyr::mutate(mut_diff = rna_mut - prot_mut)  
  


# coefficients_all %>% 
#   tibble::rownames_to_column('gene') %>% 
#   reshape2::melt() %>% 
#   dplyr::mutate(
#     assay = case_when(str_detect(variable, 'protein') ~ 'protein', .default = 'rna')
#   ) %>% 
#   dplyr::mutate(
#     predictor = case_when(str_detect(variable, 'tot_cna') ~ 'cna', .default = 'mut_status')) %>% head
  
genes_not_explained_expression = coefficients_all %>% 
  tibble::rownames_to_column('gene') %>%
  as_tibble() %>%
  dplyr::filter(protein_expression_tot_cna == 0, protein_expression_mutation_statusMutated == 0, rna_expression_tot_cna == 0, rna_expression_mutation_statusMutated == 0) %>% 
  dplyr::pull(gene) %>% 
  unique
saveRDS(genes_not_explained_expression, 'data/genes_not_explained_expression.rds')


# visualizations
rna_vs_prot %>% 
  dplyr::filter(hgnc_symbol == 'APC', mut_consequence != 'stop_gained') %>% 
  tidyr::pivot_longer(cols = c(protein_expression, rna_expression), names_to = 'assay', values_to = 'expression') %>%
  ggplot(aes(y = expression, x = tot_cna, color = mut_consequence)) +
  geom_point() + 
  facet_wrap(mutation_status~assay)

# plot the coefficients 
coefficients_all %>% 
  ggplot(aes(rna_expression_tot_cna, protein_expression_tot_cna)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  geom_smooth(method = 'lm') + 
  stat_poly_eq(use_label('eq')) +  
  ylab('Protein cna coefficients') + 
  xlab('RNA cna coefficients') + 
  xlim(c(-1.5, 2.5)) + 
  ylim(c(-1.5, 2.5))
  

# pca_coeffs = prcomp(coefficients_all, scale = TRUE)
# scree_plot = fviz_eig(pca_coeffs, addlabels = T)
# 
# library(plotly)
# library(stats)


# components <- pca_coeffs[["x"]]
# components <- data.frame(components)
# components = components %>% 
#   rownames_to_column('gene')
# # components <- cbind(components, iris$Species)
# components$PC2 <- -components$PC2
# explained_variance <- summary(pca_coeffs)[["sdev"]]
# explained_variance <- explained_variance[1:2]
# comp <- pca_coeffs[["rotation"]]
# comp[,'PC2'] <- - comp[,'PC2']
# loadings <- comp
# for (i in seq(explained_variance)){
#   loadings[,i] <- comp[,i] * explained_variance[i]
# }
# 
# features =colnames(coefficients_all)
# 
# fig <- plot_ly(components, x = ~PC1, y = ~PC2, type = 'scatter', mode = 'markers') %>%
#   layout(
#     legend=list(title=list(text='color')),
#     plot_bgcolor = "#e5ecf6",
#     xaxis = list(
#       title = "0"),
#     yaxis = list(
#       title = "1"))
# for (i in seq(4)){
#   fig <- fig %>%
#     add_segments(x = 0, xend = loadings[i, 1], y = 0, yend = loadings[i, 2], line = list(color = 'black'),inherit = FALSE, showlegend = FALSE) %>%
#     add_annotations(x=loadings[i, 1], y=loadings[i, 2], ax = 0, ay = 0,text = features[i], xanchor = 'center', yanchor= 'bottom')
# }
# 
# fig

