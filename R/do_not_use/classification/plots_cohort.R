# visualization and analysis
rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
library(tidyverse)
library(ggExtra)
library(ComplexHeatmap)
library(ggalluvial)
library(patchwork)
library(ggh4x)
source('organoids_analysis/R/classification/utils.R')

test_expr_rna = readRDS('data/test_expr_rna_bootstrap_v4.rds') %>% 
  mutate(Assay = 'RNA')
test_expr_prot = readRDS('data/test_expr_prot_bootstrap_v4.rds') %>% 
  mutate(Assay = 'Protein') 

test_expr_prot = classify_elements(test_expr_prot, pth = 0.1)
test_expr_rna = classify_elements(test_expr_rna, pth = 0.1)

test_all = bind_rows(test_expr_prot, test_expr_rna) %>% 
  mutate(Assay = factor(Assay, levels = c('RNA', 'Protein')))

# define the proportions by gene for both omics
propotion_rna_by_gene = get_props(test_expr_rna, group = 'hgnc_symbol')
propotion_prot_by_gene = get_props(test_expr_prot, group = 'hgnc_symbol')

# test_expr_rna %>% 
#   filter(is.na(cls_dosage))

# summary statistics across cohort ---------------------

## proportions by gene -----

# classes proportions across genes
props_by_gene = full_join(propotion_rna_by_gene, propotion_prot_by_gene, suffix = c('_prot', '_rna'), 
                          by = join_by('hgnc_symbol' == 'hgnc_symbol', 
                                       'class' == 'class'))

props_by_gene_plot = props_by_gene %>%
  pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'assay', values_to = 'prop') %>% 
  mutate(assay = gsub('prop_', '', assay)) %>% 
  mutate(assay = factor(assay, levels = c('rna', 'prot'))) %>% 
  ggplot(aes(x = class, y = prop, fill = assay)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'right', 
        legend.title = element_text(hjust = 0.5, 
                                    vjust = 0)
        ) + 
  xlab('Dosage class') + 
  ylab('Proportion of genes') +
  scale_fill_manual(values = c('rna' = 'steelblue', 'prot' = 'goldenrod'), 
                    labels = c("rna" = "RNA", "prot" = "Protein")) +
  guides(fill = guide_legend(
    title = 'Assay', 
    ncol = 1, 
    title.position = 'top'
  ))
saveRDS(props_by_gene_plot, 'data/plots_rds/props_by_gene_plot.rds')
ggsave('res/props_by_gene_plot.png', props_by_gene_plot, height = 6)
ggsave('res/props_by_gene_plot.pdf', props_by_gene_plot, height = 6)

## props across samples ---- 
propotion_rna_by_sample = get_props(test_expr_rna, group = 'sample', filter = F)
propotion_prot_by_sample = get_props(test_expr_prot, group = 'sample', filter = F)

props_by_sample = full_join(propotion_rna_by_sample, propotion_prot_by_sample, suffix = c('_prot', '_rna'), 
                            by = join_by('sample' == 'sample', 
                                         'class' == 'class'))

props_by_samples_plot = plot_cls_by_sample(props_by_sample, dosage_colors = dosage_colors)
saveRDS(props_by_samples_plot, 'data/plots_rds/props_by_samples_plot.rds')
ggsave('res/props_by_samples_plot.png', props_by_samples_plot, height = 7, width = 12, dpi = 300)
ggsave('res/props_by_samples_plot.pdf', props_by_samples_plot, height = 6)

## prevalent classes ----
rna_prev_class = get_prevalent_class(test_expr_rna) %>% 
  mutate(Assay = 'RNA')
prot_prev_class = get_prevalent_class(test_expr_prot) %>% 
  mutate(Assay = 'Protein')

prev_classes_omics = bind_rows(rna_prev_class, prot_prev_class) %>% 
  mutate(Assay = factor(Assay, levels = c('RNA', 'Protein')))
plot_prev = plot_prevalent_class(prev_classes_omics, classes_cols = dosage_colors, omics_cols = omics_cols)
ggsave('res/plot_prevalent_class.png', plot = plot_prev, width = 15)
ggsave('res/plot_prevalent_class.pdf', plot = plot_prev, width = 15)
saveRDS(plot_prev, 'data/plots_rds/plot_prevalent_classes.rds')

# check differences among tsg and etc

gene_functions = test_all %>% 
  dplyr::select(hgnc_symbol, CGC_role_COAD) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(CGC_role_COAD != 'None') %>% 
  dplyr::filter(!is.na(CGC_role_COAD)) 

prev_classes_omics_by_function = right_join(prev_classes_omics, gene_functions) %>% 
  filter(!is.na(Assay))

# prev_classes_omics_by_function %>% 
#   group_by(hgnc_symbol, Assay) %>% 
#   dplyr::slice_max(prop) %>% 
#   dplyr::group_by(class, Assay) %>% 
#   count() %>% 
#   ggplot(aes(x = reorder(class, +n), y = n, fill = class)) + 
#   geom_bar(stat = 'identity') + 
#   scale_fill_manual(values = dosage_colors) + 
#   theme_bw() + 
#   light_theme() + 
#   labs(x = 'Dosage classes', y = 'Number of genes', 
#        title = 'Prevalent class') +
#   coord_flip() + 
#   facet_wrap(~Assay) + 
#   # facet_wrap2(~Assay, 
#   #             strip = strip_themed(
#   #               background_x = elem_list_rect(fill = omics_cols, alpha = 0.1)
#   #             )) +
#   guides(fill = guide_legend(title = 'Dosage classes', 
#                              ncol = 4, 
#                              title.position = 'top'))

plot_classes_by_function = prev_classes_omics_by_function %>% 
  group_by(hgnc_symbol, Assay) %>% 
  dplyr::slice_max(prop) %>% 
  dplyr::group_by(class, Assay, CGC_role_COAD) %>% 
  count() %>% 
  ggplot(aes(x = reorder(class, +n), y = n, fill = class)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = dosage_colors) + 
  theme_bw() + 
  light_theme() + 
  labs(x = 'Dosage classes', y = 'Number of genes', 
       title = 'Prevalent class') +
  coord_flip() + 
  facet_grid(CGC_role_COAD~Assay) + 
  # facet_wrap2(~Assay, 
  #             strip = strip_themed(
  #               background_x = elem_list_rect(fill = omics_cols, alpha = 0.1)
  #             )) +
  guides(fill = guide_legend(title = 'Dosage classes', 
                             ncol = 3, 
                             title.position = 'top'))
saveRDS(plot_classes_by_function, 'data/plots_rds/plot_classes_by_function.rds')
ggsave('res/plot_classes_by_function.png', plot = plot_classes_by_function, width = 15)
ggsave('res/plot_classes_by_function.pdf', plot = plot_classes_by_function, width = 15)
# check how classes are moving between rna and protein --> sankey plot + visualise how many are in differnt classes

## differences among genes ---- 

class_merge = merging_omics(prot_cls = test_expr_prot, rna_cls = test_expr_rna) 

merge_plot = class_merge %>% 
  filter(!is.na(cls_comparison)) %>% 
  filter(alteration != 'Wild-type, no CNA') %>% 
  group_by(cls_comparison) %>% 
  count(alteration) %>% 
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(
    y = reorder(cls_comparison, -n), x = n, fill = alteration
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw()+ 
  scale_fill_manual(values = alt_colors) + 
  labs(x = 'Proportions', 
       y = 'Classes comparison') + 
  coord_flip() + 
  theme(legend.position = 'bottom', legend.title = element_text(hjust = 0.5, vjust = 1)) + 
  guides(fill = guide_legend(title = 'Alteration types', title.position = 'top'))
  
merge_plot
saveRDS(merge_plot, 'data/plots_rds/cls_merge_plots.rds')
ggsave('res/cls_merge_plot.png', merge_plot, height = 6)
ggsave('res/cls_merge_plot.pdf', merge_plot, height = 6)

## sankey plot ----
# focusing on the genes that have different classes
classes_diff_sankey = class_merge %>% 
  filter(cls_comparison == 'Different classes') %>% 
  group_by(class_RNA, class_Prot) %>% 
  count(name = 'n_genes') %>% 
  mutate(
    alluvium = interaction(class_RNA, class_Prot)
  ) %>% 
  pivot_longer(cols = starts_with('class_'), values_to = 'class', names_to = 'assay') %>%
  mutate(assay = gsub('class_', '', assay)) %>% 
  mutate(assay = factor(assay, levels = c('RNA', 'Prot'))) %>% 
  ggplot(aes(
    x = assay, 
    y = n_genes, 
    # fill = class, 
    stratum = class, 
    alluvium = alluvium
  )) + 
  geom_stratum(aes(fill = class)) +
  geom_flow(color = 'black') +
  # geom_flow() + 
  theme_bw() + 
  scale_fill_manual(values = dosage_colors) + 
  labs(
    x = '', 
    y = 'Number of genes with different classes'
  ) + 
  theme(legend.position = 'bottom', legend.title = element_text(hjust = 0.5, vjust = 1)) + 
  guides(fill = guide_legend(title = 'Classes', title.position = 'top',  ncol = 3)) 
saveRDS(classes_diff_sankey, 'data/plots_rds/classes_diff_sankey_plot.rds')
ggsave(plot = classes_diff_sankey, 'res/classes_diff_sankey_plot_v2.pdf', height = 7, width = 12)
ggsave(plot = classes_diff_sankey, 'res/classes_diff_sankey_plot_v2.png', height = 7, width = 12)

# Assembly 
plt_layout = '
AAAA
BBCC
DDEE
'

wrap_plots(
  list(
    props_by_gene_plot,# + theme(axis.title.y = 
    # element_text(margin = margin(l = -2, unit = 'pt')))),
    props_by_samples_plot,
    plot_prev,
    merge_plot,
    classes_diff_sankey
  ), design = plt_layout
) +
  plot_annotation(tag_levels = 'A')# & theme(legend.box.margin = margin(-5, -5, -5, -5),
#        legend.margin = margin(0, 0, 0, 0))
ggsave('res/classification_fig1.pdf', width = 15, height = 22)



## heatmap ----
drivers = c((test_expr_prot %>% filter(is_driver_intogen == 'True') %>% pull(hgnc_symbol) %>% unique), 
            test_expr_rna %>% filter(is_driver_intogen == 'True') %>% pull(hgnc_symbol) %>% unique) %>% unique

cgc_coad = c((test_expr_prot %>% filter(!is.na(CGC_role_COAD), CGC_role_COAD != 'None') %>% pull(hgnc_symbol) %>% unique), 
             test_expr_rna %>% filter(!is.na(CGC_role_COAD), CGC_role_COAD != 'None') %>% pull(hgnc_symbol) %>% unique) %>% unique

dr = readRDS('data/drivers.rds') %>% pull(gene) %>% unique

dr_heatmap = plot_heatmap(rna = test_expr_rna, 
             prot = test_expr_prot, 
             genes = dr, 
             cols = ht_colors, 
             annotation_colors = annotation_colors, 
             levels_cls = names(dosage_colors), 
             title = 'Cohort drivers'
             )

pdf('res/oncoprint_drivers_heatmap_classes.pdf', width = 12, height = 13)
dr_heatmap
dev.off()


# summary gene specific -----
## get the prelevant class for each sample -----

pdf('res/drivers_classes.pdf', width = 10, height = 8)
lapply(dr, function(x){
  plot_classes_by_gene(props_by_gene, gene = x, dosage_colors)
})
graphics.off()

## plot expression -------

pdf('res/drivers_expression_wt_only.pdf', width = 13, height = 10)
lapply(dr, function(x){
  
  p = tryCatch({
   print( 
     plot_expr_gene_by_cls(test_all, gene = x, ploidy_seq = seq(2:6), color_by_classes = T, filter_wt = FALSE))}, 
    error = function(e) {ggplot()}
   )
  
  # rna = plot_expr_gene_by_cls(test_expr_rna, gene = x, which = 'RNA', ploidy_seq = seq(2:6)) +
  #   ggtitle('RNA')
  # prot = plot_expr_gene_by_cls(test_expr_prot, gene = x, which = 'Protein', ploidy_seq = seq(2:6)) + 
  #   ggtitle('Protein')
  # 
  # (rna / prot) + 
  #   plot_annotation(title = x) +
  #   plot_layout(guides = 'collect') & 
  #   theme(legend.position = 'bottom')
  print(p)
  
})
graphics.off()

# trying putting everything together ############################################################################################

cgc_coad_genes = test_all %>% 
  filter(!is.na(CGC_role_COAD)) %>% 
  filter(CGC_role_COAD != 'None')
  
cgc_coad_genes = cgc_coad_genes %>% 
  filter(mutation_status == "Wild-type") %>% 
  group_by(Assay) %>% 
  nest() %>% 
  deframe()

cgc_coad_genes_scaled = lapply(cgc_coad_genes, function(x) {
  genes_expression = x %>% 
    select(sample, hgnc_symbol, expression) %>% 
    pivot_wider(names_from = hgnc_symbol, values_from = expression) %>% 
    tibble::column_to_rownames('sample')
    
  genes_expression_scaled = apply(genes_expression, 2, function(s) {scale(s, center=T,scale=T)})
  rownames(genes_expression_scaled) = rownames(genes_expression)
  genes_expression_scaled = genes_expression_scaled %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('sample')
   
  x %>% 
    select(-hgnc_symbol, -expression) %>% 
    full_join(., genes_expression_scaled, by = 'sample') %>% 
    pivot_longer(., (x$hgnc_symbol %>% unique), names_to = 'gene', values_to = 'expression')
  
}) %>% 
  bind_rows(., .id = "Assay")

cgc_coad_genes_scaled %>% 
    # mutate(expression_scaled = scale(expression)) %>% 
  ggplot(aes(x = tot_cna,
             y = expression, 
             color = class, 
             shape = mutation_status)) + 
  geom_point() + 
  facet_grid(Assay~CGC_role_COAD, scales = 'free') + 
  theme_bw() + 
  scale_color_manual(values = dosage_colors) + 
  theme(legend.position = 'bottom') + 
  guides(color = guide_legend(ncol = 2, title.position = 'top'), 
         shape = guide_legend(title.position = 'top'))

ggsave('res/test_all.png', width = 10)

x = test_expr_prot %>% 
  filter(hgnc_symbol == 'KRAS') %>% 
  mutate(Assay = 'Protein') %>% 
  bind_rows(., 
            (test_expr_rna %>% filter(hgnc_symbol == 'KRAS') %>% mutate(Assay = 'RNA')))

plot_expr_gene_by_cls(x, gene = 'KRAS', ploidy_seq = seq(2:6)) 

plot_expr_gene_by_cls(test_expr_rna, gene = 'TP53', which = 'RNA', ploidy_seq = seq(2:6)) 
p_prot = plot_expr_gene_by_cls(test_expr_prot, gene = 'CFB', which = 'Protein')




# plot_prevalent_class(rna_prev_class, which = 'RNA', dosage_colors)
# plot_prevalent_class(prot_prev_class, which = 'Protein', dosage_colors)

# rna_prev_class %>% 
#   slice_max(prop) %>% 
#   split(.$cls_dosage) 
#   
# to remove #########################################################################################################
# plot expression of specific genes indicating the class
test_expr_rna %>% 
  filter(cls_dosage == 'Dose_sensitive') %>% 
  pull(hgnc_symbol) %>% 
  unique



(p_rna / p_prot) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
ggsave('res/cfb_expression_class.png', width = 10)

p_rna = plot_expr_gene_by_cls(test_expr_rna, gene = '', which = 'RNA') 
p_prot = plot_expr_gene_by_cls(test_expr_prot, gene = 'MAP2K1', which = 'Protein')

(p_rna / p_prot) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
ggsave('res/msh6_expression_class.png', width = 10)


test_expr_rna %>% 
  filter(cls_dosage == 'Dose_sensitive') %>%
  filter(tot_cna != 2) %>%
  filter(CGC_role_PANCANCER != 'None') %>% 
  pull(hgnc_symbol) %>% unique

