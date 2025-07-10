.libPaths()

library(glmnet)
library(tidyverse)
library(ComplexHeatmap)
library(ggpmisc)
library("factoextra")

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/comparing_models.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/coefficients_analysis.R')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/ht_utils.R')
proteogenomics_data = readRDS('data/proteogenomics_data_all_genes.rds')
transcriptomics_data = readRDS('data/transcriptomics_data_all_genes.rds')

all_genes = intersect(unique(proteogenomics_data$hgnc_symbol), unique(transcriptomics_data$hgnc_symbol))

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
                                     'is_mutated' == 'is_mutated', 
                                     'is_driver_intogen' == 'is_driver_intogen', 
                                     'CGC_role_COAD' == 'CGC_role_COAD', 
                                     'CGC_role_PANCANCER' == 'CGC_role_PANCANCER')) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(!is.na(proteomics_code)) %>% 
  dplyr::filter(hgnc_symbol %in% all_genes) %>% 
  dplyr::mutate(mutation_status = ifelse(is_mutated == TRUE, 'Mutated', 'Wild-type')) %>% 
  dplyr::mutate(mutation_status = factor(mutation_status, levels = c('Wild-type', 'Mutated'))) %>% 
  dplyr::rename(protein_expression = mean_intensity) %>% 
  dplyr::rename(rna_expression = value)

fit_multiresponse_model = readRDS('data/multiresponse_glm_lasso_v2.rds')
fit_multiresponse_model = Filter(Negate(is.null), fit_multiresponse_model)
names(fit_multiresponse_model) = lapply(fit_multiresponse_model, function(x) {x$gene}) %>% unlist

coefficients_all = lapply(fit_multiresponse_model, get_coefficients_multiresponse)
coefficients_all = lapply(names(coefficients_all), function(x) {
  coefficients_all[[x]] %>% 
    dplyr::mutate(gene = x)
}) %>% 
  bind_rows() %>% 
  relocate(gene, .before = protein_expression_tot_cna) %>% 
  tibble::column_to_rownames('gene')
# remove rows that have only zeros for simplicity of plotting
# coefficients_all = coefficients_all[!apply(coefficients_all == 0, 1, all), ]

# plot coefficients based on the effect 

# genes with a role in CGC for coad --> create a function

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
x = readRDS('data/cohort_maf_drivers.rds')
top_genes = coefficients_all %>% 
  tibble::rownames_to_column('gene') %>% 
  dplyr::filter(gene %in% (x@data$VEP.SYMBOL %>% unique)) %>% 
  tibble::column_to_rownames('gene') %>% 
  as.matrix()

ann_data = tibble(
  ann = colnames(coefficients_all)
) %>% 
  dplyr::mutate(
    assay = str_extract(ann, 'protein|rna')
  ) %>% 
  dplyr::mutate(predictor = str_extract(ann, 'cna|mutation')) %>% 
  dplyr::select(-ann)

ann_col = list(
  assay = setNames(
    nm = c('protein', 'rna'), 
    c('orange2', 'lightseagreen')
  ), 
  predictor = setNames(
    nm = c('cna', 'mutation'), 
    c('mediumvioletred', 'royalblue4')
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


colors = circlize::colorRamp2(c(0,max(coefficients_all)+0.2), c('snow', 'dodgerblue4'))

ht = Heatmap(as.matrix(coefficients_all), 
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

# trying to define some groups based on the coefficients
coefficients_all %>% 
  mutate(prot_effect = case_when(protein_expression_tot_cna > protein_expression_mutation_statusMutated ~ 'cna_protein_stronger', 
                   .default = 'mut_protein_stronger')) %>% 
  pull(prot_effect) %>% 
  table









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
  

pca_coeffs = prcomp(coefficients_all, scale = TRUE)
scree_plot = fviz_eig(pca_coeffs, addlabels = T)

library(plotly)
library(stats)


components <- pca_coeffs[["x"]]
components <- data.frame(components)
components = components %>% 
  rownames_to_column('gene')
# components <- cbind(components, iris$Species)
components$PC2 <- -components$PC2
explained_variance <- summary(pca_coeffs)[["sdev"]]
explained_variance <- explained_variance[1:2]
comp <- pca_coeffs[["rotation"]]
comp[,'PC2'] <- - comp[,'PC2']
loadings <- comp
for (i in seq(explained_variance)){
  loadings[,i] <- comp[,i] * explained_variance[i]
}

features =colnames(coefficients_all)

fig <- plot_ly(components, x = ~PC1, y = ~PC2, type = 'scatter', mode = 'markers') %>%
  layout(
    legend=list(title=list(text='color')),
    plot_bgcolor = "#e5ecf6",
    xaxis = list(
      title = "0"),
    yaxis = list(
      title = "1"))
for (i in seq(4)){
  fig <- fig %>%
    add_segments(x = 0, xend = loadings[i, 1], y = 0, yend = loadings[i, 2], line = list(color = 'black'),inherit = FALSE, showlegend = FALSE) %>%
    add_annotations(x=loadings[i, 1], y=loadings[i, 2], ax = 0, ay = 0,text = features[i], xanchor = 'center', yanchor= 'bottom')
}

fig

