# visualization and analysis
rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
library(tidyverse)
library(ggExtra)
library(ComplexHeatmap)
source('organoids_analysis/R/classification/utils.R')

test_expr_rna = readRDS('data/test_expr_rna_bootstrap_v2.rds') 
test_expr_prot = readRDS('data/test_expr_prot_bootstrap_v2.rds') 
# define the proportions by gene for both omics
propotion_rna_by_gene = get_props(test_expr_rna, group = 'hgnc_symbol')
propotion_prot_by_gene = get_props(test_expr_prot, group = 'hgnc_symbol')

# test_expr_rna %>% 
#   filter(is.na(cls_dosage))

# summary statistics across cohort ---------------------

# classes proportions across genes
props_by_gene = full_join(propotion_rna_by_gene, propotion_prot_by_gene, suffix = c('_prot', '_rna'), 
                          by = join_by('hgnc_symbol' == 'hgnc_symbol', 
                                       'cls_dosage' == 'cls_dosage'))

props_by_gene_plot = props_by_gene %>%
  pivot_longer(cols = c(prop_prot, prop_rna), names_to = 'assay', values_to = 'prop') %>% 
  mutate(assay = gsub('prop_', '', assay)) %>% 
  mutate(assay = factor(assay, levels = c('rna', 'prot'))) %>% 
  ggplot(aes(x = cls_dosage, y = prop, fill = assay)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'bottom', 
        legend.title = element_text(hjust = 0.5, 
                                    vjust = 0)
        ) + 
  xlab('Dosage class') + 
  ylab('Proportion of genes') +
  scale_fill_manual(values = c('rna' = 'steelblue', 'prot' = 'goldenrod'), 
                    labels = c("rna" = "RNA", "prot" = "Protein")) +
  guides(fill = guide_legend(
    title = 'Assay', 
    ncol = 2, 
    title.position = 'top'
  ))
saveRDS(props_by_gene_plot, 'data/plots_rds/props_by_gene_plot.rds')
ggsave('res/props_by_gene_plot.png', props_by_gene_plot, height = 6)
ggsave('res/props_by_gene_plot.pdf', props_by_gene_plot, height = 6)

# props across samples
propotion_rna_by_sample = get_props(test_expr_rna, group = 'sample', filter = F)
propotion_prot_by_sample = get_props(test_expr_prot, group = 'sample', filter = F)

props_by_sample = full_join(propotion_rna_by_sample, propotion_prot_by_sample, suffix = c('_prot', '_rna'), 
                            by = join_by('sample' == 'sample', 
                                         'cls_dosage' == 'cls_dosage'))

props_by_samples_plot = plot_cls_by_sample(props_by_sample, dosage_colors = dosage_colors)
saveRDS(props_by_samples_plot, 'data/plots_rds/props_by_samples_plot.rds')
ggsave('res/props_by_samples_plot.png', props_by_samples_plot, height = 6)
ggsave('res/props_by_samples_plot.pdf', props_by_samples_plot, height = 6)

# heatmap
drivers = c((test_expr_prot %>% filter(is_driver_intogen == 'True') %>% pull(hgnc_symbol) %>% unique), 
            test_expr_rna %>% filter(is_driver_intogen == 'True') %>% pull(hgnc_symbol) %>% unique) %>% unique


# check how classes are moving between rna and protein
class_merge = merging_omics(prot_cls = test_expr_prot, rna_cls = test_expr_rna)

# get the 
plot_expr_gene_by_cls(test_expr_rna, gene = 'CFB', which = 'RNA', ploidy_seq = seq(2:6)) 
p_prot = plot_expr_gene_by_cls(test_expr_prot, gene = 'CFB', which = 'Protein')

rna_prev_class = get_prevalent_class(test_expr_rna)
prot_prev_class = get_prevalent_class(test_expr_prot)

rna_prev_class %>% 
  slice_max(prop) %>% 
  split(.$cls_dosage) 
  
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

