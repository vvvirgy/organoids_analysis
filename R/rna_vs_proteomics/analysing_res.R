rm(list=ls())
.libPaths()
library(tidyverse)
library(gprofiler2)
library(UpSetR)
library(grid)
library(patchwork)
library(ggalluvial)
library(ComplexHeatmap)
library(rcartocolor)
setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj")

source('organoids_analysis/R/plot_utils/dge_utils_and_plots.R')
# source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/dge_utils_and_plots.R')

# multi_omics = readRDS('/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results/lfc_prot_and_rna_bind.rds')
multi_omics = readRDS('data/lfc_prot_and_rna_bind.rds')

# volcano plot of expression
pth = .05
fth = .75

multi_omics = classify_genes(multi_omics)
multi_omics = multi_omics %>% 
  mutate(omic = factor(omic, levels = c('RNA', 'protein')))

volcano = plot_volcano(multi_omics, omic = c('RNA', 'protein'), cols = expr_cols) + 
  xlab('log FC')

# add the protein and rna facet colors
ggsave(plot = volcano, filename = 'res/volcano_rna_prot.png', width = 10,  height = 6)
ggsave(plot = volcano, filename = 'res/volcano_rna_prot.pdf', width = 10,  height = 6)

# summary statistics --> how many DEGs per karyotype and omic?

barplot_multiomics = multi_omics %>%
  filter(!is.na(lfc)) %>%
  filter(significance == 'significant') %>% 
  group_by(omic, karyotype, fc_cls) %>%
  ggplot(aes(fc_cls, fill = fc_cls)) +
  geom_bar(stat = 'count', position = 'dodge') +
  theme_bw() + 
  facet_grid(omic~karyotype, scales = 'free') + 
  scale_fill_manual(values = expr_cols) +
  coord_flip() + 
  guides(fill = guide_legend(title = 'FC class')) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  labs(x = 'FC class', 
       y = 'Number of genes') 
ggsave('res/number_of_degs_omic.png', width = 10, height = 6)
ggsave('res/number_of_degs_omic.png', width = 10, height = 6)

pt_l = 'AAAA
        AAAA
        BBBB'
(volcano / barplot_multiomics) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(design = pt_l)
  
ggsave('res/de_results_summary.png', width = 10, height = 7)
ggsave('res/de_results_summary.pdf', width = 10, height = 7)

# checking how many genes do have the fc for each karyotype

pdf('res/upset_de_multiomics.pdf', width = 10, height = 6)
multi_omics %>% 
  filter(!is.na(lfc)) %>% 
  select(name, karyotype, omic) %>% 
  mutate(data = 1) %>% 
  pivot_wider(names_from = c(omic, karyotype), names_sep = '_',  values_from = data, values_fill = 0) %>% 
  as.data.frame() %>% 
  upset(nsets = 10, order.by = 'freq', 
        sets.x.label = "Genes estimated FC (sign and ns)") 
dev.off()  


# plotting lfc for each gene (then add the CI)
multi_omics_plot = multi_omics %>% 
  filter(!is.na(lfc)) %>% 
  mutate(karyotype = factor(karyotype, levels = c('1:0', '2:0', '2:1', '2:2'))) 

multi_omics_fc_plt(multi_omics_plot, gene = 'PIK3CA', colors = c('RNA' = '#296374', 'protein' = '#D02752'))
multi_omics_fc_plt(multi_omics_plot, gene = 'KRAS', colors = c('RNA' = '#296374', 'protein' = '#D02752'))
multi_omics_fc_plt(multi_omics_plot, gene = 'SMAD2', colors = c('RNA' = '#296374', 'protein' = '#D02752'))
ggsave('res/genes_fc/smad2_fc.png', width = 8, height = 6)
multi_omics_fc_plt(multi_omics_plot, gene = 'AKT1', colors = c('RNA' = '#296374', 'protein' = '#D02752'))
ggsave('res/genes_fc/akt1_fc.png', width = 8, height = 6)

# might add the confidence intervals or a * to indicate the significance

# checking how many results do we have per group

# plot the fc 
multi_omics_cls = multi_omics %>% 
  filter(!is.na(lfc)) %>% 
  # select(name, adj_pval, lfc, karyotype, omic) %>%
  select(name, significance, fc_cls, karyotype, omic, adj_pval, lfc) %>%
  pivot_wider(names_from = omic, values_from = c(adj_pval, lfc, significance, fc_cls)) %>% 
  mutate(fc_cls = case_when(
    (lfc_RNA < 0 & lfc_protein < 0) ~ 'both negative',  
    (lfc_RNA > 0 & lfc_protein > 0) ~ 'both positive', 
    (sign(lfc_RNA) != sign(lfc_protein)) ~ 'opposite', 
    (is.na(lfc_RNA) & !is.na(lfc_protein) & lfc_protein < 0) ~ 'only protein estimated (negative)', 
    (is.na(lfc_RNA) & !is.na(lfc_protein) & lfc_protein > 0) ~ 'only protein estimated (positive)', 
    (!is.na(lfc_RNA) & is.na(lfc_protein) & lfc_RNA > 0) ~ 'only RNA estimated (positive)', 
    (!is.na(lfc_RNA) & is.na(lfc_protein) & lfc_RNA < 0) ~ 'only RNA estimated (negative)'
  )) %>% 
  mutate(sign_RNA = ifelse(adj_pval_RNA <= pth, 'significant', 'ns')) %>% 
  mutate(sign_prot = ifelse(adj_pval_protein <= pth, 'significant', 'ns'))

plot_omics_comparison(multi_omics_cls, filter = T, bg_colors = bg_cols)
ggsave('res/rna_vs_protein_fc_significant_only.png', width = 10, height = 10)  
ggsave('res/rna_vs_protein_fc_significant_only.pdf', width = 10, height = 10) 

plot_omics_comparison(multi_omics_cls, filter = F, bg_colors = bg_cols)
ggsave('res/rna_vs_protein_fc_all.png', width = 10, height = 10)  
ggsave('res/rna_vs_protein_fc_all.pdf', width = 10, height = 10) 

  
# distribution of genes in classes
multi_omics_cls %>%
  mutate(significance = paste0('RNA ', sign_RNA, ', ', 'protein ', sign_prot)) %>% 
  mutate(fc_cls = factor(fc_cls, levels = c('both positive', 'both negative', 'opposite', 
                                            'only RNA estimated (positive)', 'only RNA estimated (negative)', 
                                            'only protein estimated (positive)', 'only protein estimated (negative)') %>% rev)) %>% 
  ggplot(aes(fc_cls, fill = significance)) + 
  geom_bar(stat = 'count') + 
  coord_flip() + 
  # scale_fill_brewer(palette = 'Set3') + 
  ggsci::scale_fill_lancet() + 
  # theme_bw() + 
  theme_light() + 
  labs(
    x = 'FC classes', 
    y = 'number of genes'
  ) + 
  facet_wrap(~karyotype) + 
  theme(legend.position = 'bottom')
ggsave('res/fc_cls_rna_protein_by_karyotype.pdf',width = 12, height = 8)
ggsave('res/fc_cls_rna_protein_by_karyotype.png',width = 12, height = 8)

# considering only genes that are significant for each omic state 
plot_densities(multi_omics, filter = T) 
ggsave('res/fc_significant_density.png', width = 10, height = 6)
ggsave('res/fc_significant_density.pdf', width = 10, height = 6)

plot_barplots(multi_omics, cols = expr_cols)

# fix karyotype and see which are the genes that are significant

res_by_karyo = multi_omics %>% 
  filter(!is.na(lfc)) %>% 
  # mutate(fc_cls = ifelse(is.na(fc_cls), 'not differential', fc_cls)) %>% 
  mutate(fc_cls = case_when(
    is.na(fc_cls) ~ 'not differential', 
    .default = fc_cls
  )) %>% 
  mutate(class = paste(fc_cls, significance, sep = ', ')) %>% 
  select(name, karyotype, omic, class) %>% 
  mutate(data = 1) %>% 
  pivot_wider(names_from = c(omic, class), names_sep = ' ',  values_from = data, values_fill = 0) %>% 
  as.data.frame() %>% 
  # select(name, karyotype, omic) %>% 
  split(.$karyotype)

lapply(res_by_karyo, function(x) {
  
  karyo = gsub(':', '_', x$karyotype %>% unique)
  
  pdf(paste0('res/upsets/karyo_', karyo, '_fixed_among_omics.pdf'), width = 12, height = 7)
  grid.newpage()
  print(x %>% 
    select(-karyotype) %>% 
    UpSetR::upset(nsets = 10, order.by = 'freq', 
                  sets.x.label = paste("Number of genes for", paste0('karyotype ', x$karyotype %>% unique)), 
                  matrix.color = '#4988C4', 
                  shade.color = '#C9CDCF',
                  main.bar.color = '#0F2854', 
                  set_size.show = T))
  
  # grid.text(paste0('karyotype ', x$karyotype %>% unique), 
  #           x = 0.5, y = .95, 
  #           gp = gpar(fontsize = 16, fontface = 'bold'))
  
  graphics.off()
  
})

# fix the omic, and compare the karyotypes

res_by_omic = multi_omics %>% 
  filter(!is.na(lfc)) %>% 
  mutate(fc_cls = case_when(
    is.na(fc_cls) ~ 'not differential', 
    .default = fc_cls
  )) %>% 
  mutate(class = paste(fc_cls, significance, sep = ', ')) %>% 
  select(name, karyotype, omic, class) %>% 
  mutate(data = 1) %>% 
  pivot_wider(names_from = c(class, karyotype), names_sep = ' ',  values_from = data, values_fill = 0) %>% 
  as.data.frame() %>% 
  # select(name, karyotype, omic) %>% 
  split(.$omic)

lapply(res_by_omic, function(x) {
  
  omic = x$omic %>% unique
  
  metadata = tibble(
    sets = colnames(x)[-c(1,2)]
  ) %>% 
    mutate(cls = str_extract(pattern = 'ns|significant', sets))
  
  pdf(paste0('res/upsets/', omic, '_fixed_among_karyotypes.pdf'), width = 12, height = 7)
  
  print(x %>% 
          select(-omic) %>% 
          UpSetR::upset(nsets = 16, order.by = 'freq', 
                        sets.x.label = paste("Number of genes for", omic), 
                        matrix.color = '#4988C4', 
                        shade.color = '#C9CDCF',
                        main.bar.color = '#0F2854', 
                        set_size.show = T, 
                        set.metadata = 
                          list(data = metadata, 
                               plots = list(list(type = 'matrix_rows', 
                                                 column = 'cls', 
                                                 colors = c('significant' = '#9ABF80', 
                                                            'ns' = '#DDAED3')
                                                 )))))
  
  # grid.text(paste0('karyotype ', x$karyotype %>% unique), 
  #           x = 0.5, y = .95, 
  #           gp = gpar(fontsize = 16, fontface = 'bold'))
  
  graphics.off()
  
}) 


# sankey plots
multi_omics_v2 = multi_omics %>% 
  complete(
    name,
    karyotype,
    omic,
    fill = list(
      significance = NA,
      fc_cls = NA
    )
  ) 
  
plot_sankey(multi_omics_v2, strat = 'karyotype', facet = 'omic', cols = sankey_cols)
ggsave('res/omic_sankey_classes.png', width = 10)
ggsave('res/omic_sankey_classes.pdf', width = 10)
plot_sankey(multi_omics, strat = 'karyotype', facet = 'omic', cols = sankey_cols)

plot_sankey(multi_omics_v2, strat = 'omic', facet = 'karyotype', cols = sankey_cols)  
ggsave('res/karyotype_sankey_classes.png', width = 10)
ggsave('res/karyotype_sankey_classes.pdf', width = 10)


# take a look to which genes are significant --> interesting groups
# genes that are significant more or less expressed in all the karyotypes

degs_karyo_omic = multi_omics %>% 
  filter(significance == 'significant', !is.na(fc_cls)) %>% 
  # group_by(karyotype, omic) %>% 
  split(interaction(.$omic, .$karyotype, .$fc_cls, sep = "_"))
  
degs_karyo_omic = lapply(degs_karyo_omic, function(x) {
  x$name %>% unique 
})

# running enrichment

degs_enrich = lapply(degs_karyo_omic, function(x) {
  gprofiler2::gost(
    x, 
    organism = 'hsapiens', 
    sources = c('GO', 'KEGG', 'REAC', 'WP', 'TF', 'MIRNA', 'CORUM'), 
    ordered_query = F, 
    exclude_iea = T, 
    measure_underrepresentation = F, 
    evcodes = T, 
    correction_method = 'fdr', 
    domain_scope = 'annotated', 
    highlight = T
  )
})
saveRDS(degs_enrich, 'data/enrichment_results/degs_enrich.rds')

degs_enrich = readRDS('data/enrichment_results/degs_enrich.rds')

# visualising the results for the real degs

wrap_plots(list(sankey_annotations(degs_enrich, karyo = '1:0', cols = fc_cols),
                sankey_annotations(degs_enrich, karyo = '2:0', cols = fc_cols),
                sankey_annotations(degs_enrich, karyo = '2:1', cols = fc_cols),
                sankey_annotations(degs_enrich, karyo = '2:2', cols = fc_cols))) + 
  plot_annotation(title = 'Annotations among DE genes for each karyotype') + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
ggsave('res/annotation/sankeys_karyotypes_annotations.png', width = 10, height = 8)
ggsave('res/annotation/sankeys_karyotypes_annotations.pdf', width = 10, height = 8)


sankey_annotations_by_omic(degs_enrich, filter_terms = T, karyo = c('1:0', '2:0', '2:1', '2:2'), omic_list = c('RNA'), pth = .05, cols = fc_cols)
sankey_annotations_by_omic(degs_enrich, filter_terms = T, karyo = c('1:0', '2:0', '2:1', '2:2'), omic_list = c('protein'), pth = .05, cols = fc_cols)
ggsave('res/annotation/sankeys_karyotypes_annotations.pdf', width = 10, height = 8)

# pdf('res/annotation/enrichment_degs.pdf', width = 20, height = 30)
# lapply(names(degs_enrich), function(x) {
#   print(plot_enrichment_results(degs_enrich[[x]]$result, 
#                           highlight = F, 
#                           sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
#     ggtitle(gsub('_', ' ', x)))
# })
# graphics.off()

plot_enrichment_results(degs_enrich$`RNA_1:0_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_1:0_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:0_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:0_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:1_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:1_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:2_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:2_up`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)

plot_enrichment_results(degs_enrich$`RNA_1:0_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_1:0_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:0_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:0_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:1_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:1_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`protein_2:2_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)
plot_enrichment_results(degs_enrich$`RNA_2:2_down`$result, highlight = F, sources = c('REAC', 'KEGG', 'WP', 'CORUM', 'GO:BP', 'GO:MF', 'GO:CC')) + 
  facet_wrap(~source, scales = 'free', ncol = 1)

enrichment_barplot(degs_enrich, pth = .05, cols = fc_cols, source_list = c('GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG', 'WP', 'CORUM')) + 
  ggtitle('Only significant terms')
ggsave('res/annotation/enrichment_number_terms_no_tf_mirna.pdf', width = 10, height = 8)
ggsave('res/annotation/enrichment_number_terms_no_tf_mirna.png', width = 10, height = 8)

# pdf('res/annotation/fc_ht/rna_22_up.pdf', width = 4, height = 20)
# plot_fc_heatmap(res = degs_enrich$`RNA_2:2_up`, 
#                 deg = multi_omics, 
#                 # source_list = c('GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG', 'WP', 'CORUM'),
#                 source_list = c('REAC', 'KEGG', 'WP', 'CORUM'),
#                 genes_number = 20,
#                 omics = 'RNA',
#                 karyo = '2:2', 
#                 direction = 'up',
#                 ann_colors = source_colors) %>% 
#   ggplotify::as.ggplot()
# dev.off()
# 
# plot_fc_heatmap(res = degs_enrich$`RNA_1:0_down`, 
#                 deg = multi_omics, 
#                 # source_list = c('GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG', 'WP', 'CORUM'), 
#                 source_list = c('REAC', 'KEGG', 'WP', 'CORUM'), 
#                 genes_number = 10,
#                 omics = 'RNA',
#                 karyo = '1:0', 
#                 direction = 'down',
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`RNA_1:0_up`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG', 'WP', 'CORUM'),
#                 # source_list = c('REAC', 'KEGG', 'WP', 'CORUM'), 
#                 genes_number = 0,
#                 omics = 'RNA',
#                 karyo = '1:0', 
#                 direction = 'down',
#                 pth = .05,
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`protein_2:1_up`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:1', 
#                 omics = 'protein', 
#                 direction = 'up', 
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`protein_2:2_up`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP', 'CORUM'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:2', 
#                 omics = 'protein', 
#                 direction = 'up', 
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`RNA_2:2_down`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP', 'CORUM'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:2', 
#                 omics = 'RNA', 
#                 direction = 'down', 
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`RNA_2:1_down`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP', 'CORUM'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:1', 
#                 omics = 'RNA', 
#                 direction = 'down', 
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`RNA_2:1_up`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP', 'CORUM'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:1', 
#                 omics = 'RNA', 
#                 direction = 'up', 
#                 ann_colors = source_colors)
# 
# plot_fc_heatmap(res = degs_enrich$`protein_2:1_down`, 
#                 deg = multi_omics, 
#                 source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP'), 
#                 pth = .05, 
#                 genes_number = 0, 
#                 karyo = '2:1', 
#                 omics = 'protein', 
#                 direction = 'down', 
#                 ann_colors = source_colors)

plot_terms_fc(res = degs_enrich$`protein_2:1_down`, 
              deg = multi_omics, 
              source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP'), 
              pth = .05, 
              genes_number = 0, 
              karyo = '2:1', 
              omics = 'protein', 
              direction = 'down') 

plot_terms_fc(res = degs_enrich$`RNA_1:0_down`, 
              deg = multi_omics, 
              # source_list = c('GO:MF', 'GO:CC', 'GO:BP', 'REAC', 'KEGG', 'WP'), 
              source_list = 'REAC',
              pth = .05, 
              genes_number = 0, 
              karyo = '1:0', 
              omics = 'RNA') 








# different gene set --------
# genes with the same sign between rna and protein

same_sign_genes = multi_omics_cls %>% 
  filter(sign_RNA =='significant', sign_prot=='significant') %>%  # -- grep the same sign
  filter(fc_cls %in% c('both positive', 'both negative')) %>% 
  filter((lfc_RNA >= fth & lfc_protein >= fth) | (lfc_RNA <= -fth & lfc_protein <= -fth) ) %>% 
  split(interaction(.$fc_cls, .$karyotype, sep = "_"))

same_sign_genes = lapply(same_sign_genes, function(x) {
  x$name %>% unique 
})
  
same_sign_genes_enrich = lapply(same_sign_genes, function(x) {
  gprofiler2::gost(
    x, 
    organism = 'hsapiens', 
    sources = c('GO', 'KEGG', 'REAC', 'WP', 'TF', 'MIRNA', 'CORUM'), 
    ordered_query = F, 
    exclude_iea = T, 
    measure_underrepresentation = F, 
    evcodes = T, 
    correction_method = 'fdr', 
    domain_scope = 'annotated', 
    highlight = T
  )
})
saveRDS(same_sign_genes_enrich, 'data/enrichment_results/same_sign_genes_enrich.rds')

same_sign_genes_enrich = readRDS('data/enrichment_results/same_sign_genes_enrich.rds')

plot_enrichment_results(same_sign_genes_enrich$`both negative_1:0`$result, 
                        sources =  c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP', 'CORUM'), 
                        highlight = F)

plot_enrichment_results(same_sign_genes_enrich$`both positive_1:0`$result,
                        sources =  c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP', 'CORUM'), 
                        highlight = F)

plot_enrichment_results(same_sign_genes_enrich$`both negative_2:2`$result,
                        sources =  c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP', 'CORUM'), 
                        highlight = F)

plot_enrichment_results(same_sign_genes_enrich$`both negative_2:1`$result,
                        sources =  c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP', 'CORUM'), 
                        highlight = F)




# plotting 
# upset plot

# same_sign_genes_enrich = readRDS('data/enrichment_results/same_sign_genes_enrich.rds')
# enrichment_heatmap(same_sign_genes_enrich$`both positive_2:2`, highlight = T, source_list = c('REAC', 'KEGG', 'GO:MF'))


plot_fc_heatmap(same_sign_genes_enrich$`both positive_2:2`, 
                multi_omics, 
                # source_list =  c('GO:MF', 'GO:BP', 'GO:CC', 'KEGG', 'REAC', 'WP', 'CORUM'), 
                genes_number = 0,
                highlight = T, 
                karyo = '2:2', direction = 'up'
)

same_sign_terms_significant = lapply(same_sign_genes_enrich %>% names, function(x) {
  same_sign_genes_enrich[[x]]$result %>% 
    mutate(adj_pval = p.adjust(p_value, method = 'BH')) %>% 
    filter(adj_pval <= pth) %>% 
    filter(!source %in% c('TF', 'MIRNA', 'HP', 'HPA')) %>% 
    select(term_name) %>% 
    mutate(query = x)
}) %>% 
  bind_rows() %>% 
  mutate(query = gsub('_', ' ', query)) %>% 
  rename(karyotype = query) %>% 
  mutate(data = 1) %>% 
  distinct() 


same_sign_terms_significant %>% 
  pivot_wider(names_from = karyotype, values_from = data, values_fill = 0) %>% 
  as.data.frame() %>% 
  upset(nsets = 10, order.by = 'freq')
   




# running some enrichment

gene_list_LOH = multi_omics_cls %>% 
  filter(karyotype == '1:0') %>% 
  filter(fc_cls == 'both negative', sign_RNA == 'significant', sign_prot == 'significant') %>% 
  pull(name)
  
gene_list_LOH_gost = gprofiler2::gost(
  gene_list_LOH, 
  organism = 'hsapiens', 
  ordered_query = F, 
  exclude_iea = T, 
  measure_underrepresentation = F, 
  evcodes = T, 
  correction_method = 'fdr', 
  domain_scope = 'annotated', 
  highlight = T
)

gene_list_LOH_gost$result %>% 
  filter(highlighted == TRUE) %>% 
  select(significant, p_value, term_name, intersection) %>% 
  separate_rows(intersection, sep = ',') %>% 
  rename(gene = intersection) %>% 
 
# genes LOH positive 
gene_list_LOH_pos = multi_omics_cls %>% 
  filter(karyotype == '1:0') %>% 
  filter(fc_cls == 'both positive', sign_RNA == 'significant', sign_prot == 'significant') %>% 
  pull(name)

gene_list_LOH_pos_gost = gprofiler2::gost(
  gene_list_LOH_pos, 
  organism = 'hsapiens', 
  ordered_query = F, 
  exclude_iea = T, 
  measure_underrepresentation = F, 
  evcodes = T, 
  correction_method = 'fdr', 
  domain_scope = 'annotated', 
  highlight = T
)

gene_list_LOH_pos_gost_top = gene_list_LOH_pos_gost$result %>% 
  filter(highlighted) %>% 
  select(intersection, term_name) %>% 
  separate_rows(intersection, sep = ',') %>% 
  rename(gene = intersection)

multi_omics %>% 
  full_join(., gene_list_LOH_pos_gost_top, by = join_by('name' == 'gene')) %>% 
  filter(!is.na(term_name)) %>% 
  ggplot(aes(y = lfc, x = karyotype, color = term_name)) +
  geom_point() + 
  facet_grid(omic~term_name) + 
  theme_light()


# taking a look to the 2:1 negative in both 
gene_list_trisomy = multi_omics_cls %>% 
  filter(karyotype == '2:1') %>% 
  filter(fc_cls == 'both negative', sign_RNA == 'significant', sign_prot == 'significant') %>% 
  pull(name)

gene_list_trisomy_gost = gprofiler2::gost(
  gene_list_trisomy, 
  organism = 'hsapiens', 
  ordered_query = F, 
  exclude_iea = T, 
  measure_underrepresentation = F, 
  evcodes = T, 
  correction_method = 'fdr', 
  domain_scope = 'annotated', 
  highlight = T
)
gene_list_trisomy_gost$result %>% 
  filter(highlighted)

# taking a look to the 2:2 negative in both 
gene_list_tetra = multi_omics_cls %>% 
  filter(karyotype == '2:2') %>% 
  filter(fc_cls == 'both negative', sign_RNA == 'significant', sign_prot == 'significant') %>% 
  pull(name)

gene_list_tetra_gost = gprofiler2::gost(
  gene_list_tetra, 
  organism = 'hsapiens', 
  ordered_query = F, 
  exclude_iea = T, 
  measure_underrepresentation = F, 
  evcodes = T, 
  correction_method = 'fdr', 
  domain_scope = 'annotated', 
  highlight = T
)
gene_list_tetra_gost$result %>% 
  filter(highlighted) %>% 
  view

# taking a look to the 2:2 positive in both 
gene_list_tetra = multi_omics_cls %>% 
  filter(karyotype == '2:2') %>% 
  filter(fc_cls == 'both positive', sign_RNA == 'significant', sign_prot == 'significant') %>% 
  pull(name)

gene_list_tetra_gost = gprofiler2::gost(
  gene_list_tetra, 
  organism = 'hsapiens', 
  ordered_query = F, 
  exclude_iea = T, 
  measure_underrepresentation = F, 
  evcodes = T, 
  correction_method = 'fdr', 
  domain_scope = 'annotated', 
  highlight = T
)

gene_list_tetra_gost$result %>% dim

gene_list_tetra_gost$result %>% 
  filter(highlighted) %>% 
  view

  

# coherent genes
# by only one omic
rna_fc_by_karyo = multi_omics_cls %>% 
  select(name, karyotype, ends_with('RNA'), -sign_RNA) %>% 
  pivot_wider(names_from = karyotype, values_from = c(adj_pval_RNA, lfc_RNA)) 
  
coherent_genes = rna_fc_by_karyo %>% 
  filter(`lfc_RNA_2:1` > 0, `lfc_RNA_2:2` > 0, `lfc_RNA_1:0` < 0 ) %>% 
  filter(`adj_pval_RNA_2:1` <= pth | `adj_pval_RNA_2:2` <= pth) %>% 
  filter(`adj_pval_RNA_1:0` <= pth) %>% 
  select(name, starts_with('lfc')) %>% 
  pull(name)

rna_fc_by_karyo %>% 
  filter(`lfc_RNA_2:1` > 0, `lfc_RNA_2:2` > 0, `lfc_RNA_1:0` < 0 ) %>% 
  filter(`adj_pval_RNA_2:1` <= pth | `adj_pval_RNA_2:2` <= pth) %>% 
  filter(`adj_pval_RNA_1:0` <= pth) %>% 
  select(name, starts_with('lfc')) %>% 
  pivot_longer(cols = starts_with('lfc'), names_to = 'karyotype', values_to = 'lfc') %>% 
  mutate(karyotye = gsub('lfc_RNA_', '', karyotype)) %>% 
  ggplot(aes(x = karyotype, y = lfc)) + 
  geom_point() + 
  theme_bw()

coherent_genes_ann = gprofiler2::gost(coherent_genes, 
                                      organism = 'hsapiens', 
                                      ordered_query = F, 
                                      exclude_iea = T, 
                                      measure_underrepresentation = F, 
                                      evcodes = T, 
                                      correction_method = 'fdr', 
                                      domain_scope = 'annotated', 
                                      highlight = T)

ann = coherent_genes_ann$result %>% 
  filter(significant) %>% 
  
  select(term_name, intersection) %>% 
  



