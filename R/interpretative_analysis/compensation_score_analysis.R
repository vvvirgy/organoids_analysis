rm(list = ls())
.libPaths()
library(tidyverse)
library(dplyr)
library(boot)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(googlesheets4)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/dge_utils_and_plots.R')


# dividing the genes in different groups ---- 
df = readRDS('data/compensation/cs_results.rds')
karyos = c('1:0', '2:0', '2:1', '2:2')

df %>% 
  ggplot(
    aes(
      x = karyotype, 
      y = CS, 
      fill = omic
    )
  ) + 
  geom_violin() + 
  theme_bw() +
  ggpubr::stat_compare_means(comparisons = list(c('1:0', '2:1')), ) + 
  facet_wrap(vars(omic))


# setting the thr karyotype specific 
df_groups = lapply(karyos, function(x) {
  define_th(df, q = .5, karyo = x)  
})

df_groups = df_groups %>% 
  bind_rows()

# checking the effect of compensation
df_groups %>% 
  group_by(name) %>% 
  filter(n() == 4) %>% 
  # filter(name == 'PTGFRN') %>% 
  mutate(karyotype = factor(karyotype), 
         karyo_num = as.numeric(karyotype)) %>% 
  ggplot(aes(
    x = karyo_num, 
    y = RNA, 
    color = reg_group
  )) + 
  geom_point() +
  # geom_line() + 
  scale_color_manual(values = category_colors) +
  theme_bw()
                                                         
# number of genes per karyotype in each reg group
df_groups %>% 
  group_by(karyotype, reg_group) %>% 
  dplyr::count() %>% 
  ggplot(aes(
    x = reg_group, 
    y = n, 
    fill = karyotype
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Pastel2') + 
  labs(
    x = 'Regulatory groups', 
    y = 'Number of genes'
  ) + 
  theme(legend.position = 'bottom') 
ggsave('res/compensation_score/genes_cs_by_karyo.png', width = 8, height = 6)
ggsave('res/compensation_score/genes_cs_by_karyo.pdf', width = 8, height = 6)

# plot the scatter of CS per omic
category_colors <- c(
  "(RNA-prot heavy)" = "#AD002AB2",
  "(Prot-heavy)" = "#E18727B2",
  "(RNA-heavy)" = "#20854Eb2",
  "(RNA-Prot light)" = "#00468BB2",
  "Intermediate/Other" = "gainsboro"
)

df_groups %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") + 
  facet_wrap(~karyotype) + 
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2))
ggsave('res/compensation_score/scatter_cs_rna_protein_groups.png', width = 12, height = 8)


# checking haploinsuffient genes 
haploinsufficients = read.table('data/utilities/haploinsufficiency.bed', sep = '\t', header = F)
haplo_genes = haploinsufficients %>% 
  filter(V11 > .8) %>%   # dplyr::select(V4, V5, V10, V11) %>% 
  pull(V4)

haplo_cs = df_groups %>% 
  filter(name %in% haplo_genes) 

haplo_cs %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS") + 
  facet_wrap(~karyotype) + 
  guides(color = guide_legend(title = 'Regulatory group', 
                              position = 'bottom', 
                              nrow = 2))


# checking CS and effects

df_groups %>% 
  mutate(karyotype = factor(karyotype)) %>% 
  group_by(karyotype) %>% 
  filter(reg_group != 'Intermediate/Other') %>% 
  mutate(freq = row_number()/n()) %>% 
  ggplot(
    aes(x = karyotype, 
        stratum = reg_group, 
        alluvium = name,
        y = freq,
        fill = reg_group
    ) 
  ) +
  geom_flow(alpha = .7) +
  geom_stratum() +
  theme_light() 
  

# identifying the weird genes - those that change the cs or reg class
genes_groups_prop = df_groups %>% 
  group_by(name, reg_group) %>% 
  mutate(n_group = n()) %>%
  group_by(name) %>% 
  mutate(n_gene = n(), 
         prop = n_group/n_gene) 
  
cna_sensitive_genes = genes_groups_prop %>% 
  filter(prop == 1) %>% 
  filter(reg_group == 'Intermediate/Other')
  
sens_genes = cna_sensitive_genes$name %>% unique
df %>% 
  filter(name %in% sens_genes) %>%
  ggplot(aes(
    lfc
  )) + 
  # geom_point() + 
  geom_histogram() +
  facet_wrap(~omic)

sens_pathways = enrichGO(
  gene = sens_genes,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', 
  ont = 'ALL'
)

dotplot(sens_pathways, showCategory = 30)


genes_groups_prop %>% 
  filter(prop == 1) %>% 
  dplyr::select(name, reg_group) %>% 
  distinct() %>% 
  group_by(reg_group) %>% 
  dplyr::count()

rna_prot_heavy_genes = genes_groups_prop %>% 
  filter(prop == 1) %>% 
  filter(reg_group == '(RNA-prot heavy)') %>% 
  pull(name) %>% 
  unique

rna_prot_heavy_pathways = enrichGO(
  gene = rna_prot_heavy_genes,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', 
  ont = 'ALL'
)

rna_prot_heavy_pathways@result %>% view

cls_v2 = enrichplot::pairwise_termsim(rna_prot_heavy_pathways)
enrichplot::treeplot(cls_v2)

# genes that act weird
homozigous_compensating = genes_groups_prop %>% 
  filter(prop < 1) %>% 
  dplyr::select(karyotype, name, reg_group) %>% 
  pivot_wider(names_from = karyotype, values_from = reg_group) %>% 
  filter(`1:0` == '(RNA-Prot light)') %>% 
  pull(name)

homozigous_compensating_pathways = enrichGO(
  gene = homozigous_compensating,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', 
  ont = 'ALL'
)
dotplot(homozigous_compensating_pathways)

  
amp_compensating = genes_groups_prop %>% 
  filter(prop < 1) %>% 
  dplyr::select(karyotype, name, reg_group) %>% 
  pivot_wider(names_from = karyotype, values_from = reg_group) %>% 
  filter(`2:1` %in% c('(RNA-heavy)','(RNA-prot heavy)', '(Prot-heavy)') | `2:2` %in% c('(RNA-heavy)','(RNA-prot heavy)', '(Prot-heavy)')) %>% 
  pull(name)


amp_compensating_pathways = enrichGO(
  gene = amp_compensating,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', 
  ont = 'ALL'
)
dotplot(amp_compensating_pathways)

# Run enrichment -----

#  per groups - ORA ----
# extract the genes per karyo-group
genes_by_group = df_groups %>% 
  dplyr::select(karyotype, reg_group, name) %>% 
  split(interaction(.$karyotype, .$reg_group))


# map gene symbols on entrez
genes_by_group = lapply(genes_by_group, function(x) {
  
  gene_df <- bitr(
    x$name %>% unique,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )  
  
  return(gene_df)
})

# running enrichment on each gene group for karyotype
enrichment_groups = lapply(genes_by_group, function(s) {
  run_ORA(unique(s$ENTREZID), org = 'human', p_th = .1)
})
# saveRDS(enrichment_groups, 'data/compensation/cs_groups_enrichment.rds')

#  per karyotype - GSEA ----
# extract the genes per karyo-group
genes_by_group = df_groups %>% 
  dplyr::select(karyotype, name, RNA, protein) %>% 
  pivot_longer(cols = c(RNA, protein), names_to = 'omic', values_to = 'CS') %>% 
  split(interaction(.$karyotype, .$omic))

# map gene symbols on entrez and order by CS
genes_by_group_ranked = lapply(genes_by_group, function(x) {
  
  gene_df <- bitr(
    x$name %>% unique,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  dd = gene_df %>% 
    full_join(., x, by = join_by('SYMBOL' == 'name')) %>% 
    filter(!is.na(ENTREZID))
  
  setNames(nm = dd$ENTREZID, object = dd$CS) %>% 
    sort(decreasing = T)
  
})

# running enrichment on each gene group for karyotype
enrichment_groups = lapply(genes_by_group_ranked, function(s) {
  print('running')
  run_gsea(s, org = 'human', p_th = .1)
})
saveRDS(enrichment_groups, 'data/compensation/cs_gsea.rds')

# analyse the results ------

enrichment_groups = readRDS('data/compensation/cs_groups_enrichment.rds') 

## visualising results ---- ORA
enrichment_groups_res = lapply(enrichment_groups, function(x) {
  lapply(names(x), function(s) {
    
    x[[s]] %>% 
      setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% 
      as.data.frame() %>% 
      mutate(source = s)
    
  }) %>% 
    bind_rows()
})

enrichment_groups_res = lapply(enrichment_groups_res %>% names, function(g){ 
  
  enrichment_groups_res[[g]] %>% 
    mutate(group = g) %>% 
    tidyr::separate(group, into = c('karyotype', 'reg_group'), sep = '\\.', remove = F, convert = T)
  
}) %>% 
  bind_rows()

enrichment_groups_res %>% 
  filter(p.adjust <= .05) %>% 
  group_by(karyotype, reg_group) %>% 
  dplyr::count() %>% 
  ggplot(aes(
    x = reg_group, 
    y = n, 
    fill = karyotype
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Pastel2') + 
  labs(
    x = 'Regulatory groups', 
    y = 'Number of genes'
  ) +
  guides(fill = guide_legend(title = 'Karyotypes', position = 'bottom'))
ggsave('res/compensation_score/annotations_by_group.png', width = 8, height = 6)
ggsave('res/compensation_score/annotations_by_group.pdf', width = 8, height = 6)

enrichment_groups_res %>% 
  group_by(karyotype, source, reg_group) %>% 
  dplyr::count() %>% 
  view

enrichment_groups_res %>% 
  filter(reg_group == '(RNA-prot heavy)') %>%
  filter(source == 'KEGG') %>%
  filter(p.adjust <= .05) %>% 
  mutate(Gene_Ratio = str_split(GeneRatio, "/") %>%
           map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2]))) %>% 
  group_by(Description, source) %>% 
  # count() %>% 
  # arrange(desc(n))
  ggplot(aes(
    x = karyotype, 
    y = Description, 
    color = p.adjust, 
    size = Gene_Ratio
  )) + 
  geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~source, scales = 'free') + 
  theme_bw() 

ggsave('res/annotation/rna_prot_light_go_ann.png', width = 22, height = 10)
  
# mapping the CS on the pathway
enrichment_cs = enrichment_groups_res %>% 
  separate_rows(geneID, sep = '/') %>% 
  full_join(., df_groups, by = join_by(
    'geneID' == 'name', 
    'reg_group' == 'reg_group', 
    'karyotype' == 'karyotype'
  )) %>% 
  group_by(Description, karyotype, reg_group) %>% 
  mutate(mean_RNA = mean(RNA, na.rm = T), 
         mean_protein = mean(protein, na.rm = T)) %>% 
  relocate(mean_RNA, .after = Description) %>% 
  relocate(mean_protein, .after = Description) %>% 
  relocate(reg_group, .after = Description) %>% 
  filter(!is.na(Description)) %>% 
  mutate(genes = paste(geneID, collapse = '/')) %>% 
  dplyr::select(-c(geneID, RNA, protein, DNA_lfc)) %>% 
  distinct()

write.csv(enrichment_cs, 'res/compensation_score/enrichment_cs.csv', sep = ',', col.names = T, quote = F, row.names = F)
  

view(enrichment_cs)

enrichment_cs %>% 
  dplyr::select(Description, karyotype, reg_group, starts_with('mean')) %>% 
  distinct() %>% 
  pivot_longer(cols = c(mean_RNA, mean_protein), 
               names_to = 'omic', 
               values_to = 'mean_CS') %>% 
  filter(reg_group == "(RNA-Prot light)") %>% 
  group_by(omic) %>%
  slice_min(mean_CS, n = 10) %>% 
  ggplot(aes(
    y = Description, 
    x = omic, 
    color = mean_CS
  )) + 
  geom_point() + 
  scale_color_viridis_c() +
  theme_bw() #+ 
  # facet_wrap(~karyotype, scales = 'free_y')

df_wide = df %>% 
  dplyr::select(-CS) %>% 
  pivot_wider(names_from = omic, values_from = lfc)

enrichment_fc = enrichment_groups_res %>% 
  separate_rows(geneID, sep = '/') %>% 
  full_join(., df_wide, by = join_by(
    'geneID' == 'name', 
    # 'reg_group' == 'reg_group', 
    'karyotype' == 'karyotype'
  )) %>% 
  group_by(Description, karyotype, reg_group) %>% 
  mutate(mean_RNA_FC = mean(RNA, na.rm = T), 
         mean_protein_FC = mean(protein, na.rm = T)) %>% 
  relocate(mean_RNA_FC, .after = Description) %>% 
  relocate(mean_protein_FC, .after = mean_RNA_FC) %>% 
  relocate(reg_group, .after = Description) %>% 
  filter(!is.na(Description)) %>% 
  mutate(genes = paste(geneID, collapse = '/')) %>% 
  dplyr::select(-c(geneID, RNA, protein, DNA_lfc)) %>% 
  distinct()

view(enrichment_fc)
write.csv(enrichment_fc, 'res/compensation_score/enrichment_fc.csv', sep = ',', col.names = T, quote = F, row.names = F)

# compare clusters -----
prot_heavy = genes_by_group[grep('(Prot-heavy)', names(genes_by_group), value = T)] 
prot_heavy = lapply(prot_heavy, function(x) {x$SYMBOL %>% unique})

prot_heavy_clusters <- compareCluster(prot_heavy, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db, 
                              keyType = 'SYMBOL',
                              ont = "ALL")
dotplot(prot_heavy_clusters, showCategory = 20)
# sim_prot_heavy = enrichplot::pairwise_termsim(prot_heavy_clusters)

rna_heavy = genes_by_group[grep('(RNA-heavy)', names(genes_by_group), value = T)] 
rna_heavy = lapply(rna_heavy, function(x) {x$SYMBOL %>% unique})

rna_heavy_clusters <- compareCluster(rna_heavy, 
                                      fun = "enrichGO", 
                                      OrgDb = org.Hs.eg.db, 
                                      keyType = 'SYMBOL',
                                      ont = "ALL")
dotplot(rna_heavy_clusters, showCategory = 5)


homozigous_genes = genes_by_group[grep('1:0', names(genes_by_group), value = T)] 
homozigous_genes = lapply(homozigous_genes, function(x) {x$SYMBOL %>% unique})

homozigous_genes_clusters <- compareCluster(homozigous_genes, 
                                     fun = "enrichGO", 
                                     OrgDb = org.Hs.eg.db, 
                                     keyType = 'SYMBOL',
                                     ont = "ALL")
dotplot(homozigous_genes_clusters, showCategory = 5)


# comparing all clusters 

ggenes = lapply(genes_by_group[1:16], function(x) {x$SYMBOL %>% unique})
names(ggenes)
reg_groups_comparison <- compareCluster(ggenes, 
                                        fun = "enrichGO", 
                                        OrgDb = org.Hs.eg.db, 
                                        keyType = 'SYMBOL',
                                        ont = "ALL")

dotplot(reg_groups_comparison, showCategory = 10) + 
  facet_wrap(~Cluster, scales = 'free')

# cs_cluster_sim = enrichplot::pairwise_termsim(cs_cluster_enrichment)
# enrichplot::treeplot(cs_cluster_sim)

top_10_terms_by_group = reg_groups_comparison@compareClusterResult %>% 
  # filter(p.adjust <= .05) %>% 
  tidyr::separate(Cluster, into = c('karyotype', 'reg_group'), sep = '\\.') %>% 
  group_by(karyotype) %>% 
  slice_min(p.adjust, n = 10) %>% 
  dplyr::select(karyotype, reg_group, Description, p.adjust) %>% 
  split(interaction(.$karyotype, .$reg_group))

test = lapply(top_10_terms_by_group, function(x) {
  tterms = x$Description %>% unique
  
  reg_groups_comparison@compareClusterResult %>% 
    # filter(p.adjust <= .05) %>% 
    tidyr::separate(Cluster, into = c('karyotype', 'reg_group'), sep = '\\.') %>% 
    filter(Description %in% tterms)
})

test$`1:0.(RNA-heavy)` %>% 
  filter(reg_group == '(RNA-heavy)') %>% 
  ggplot(aes(
    x = karyotype, 
    y = Description, 
    color = p.adjust
  )) + 
  geom_point() + 
  theme_bw()

## similarity ----

treeplots_enrichment = lapply(enrichment_groups[1:16], function(x) {
  
  lapply(x, function(s) {
    # s = x$GO
    
    tryCatch({
      df = enrichplot::pairwise_termsim(s, method = 'Jiang')
      enrichplot::treeplot(df, group_color = rcartocolor::carto_pal(n = 5, name = 'Prism')) 
    }, 
    error = function(e) {
      NULL
    })
  })
  
})

pdf('res/compensation_score/annotations_by_group_Jiang_dist.pdf', width = 20, height = 14)
lapply(names(treeplots_enrichment), function(x) {
  print(x)
  print(assebly_plots(enrich = treeplots_enrichment[[x]], cls = x))
})
dev.off()

pdf('res/compensation_score/annotations_by_group_GO_Jiang_dist.pdf', width = 16, height = 8)
lapply(names(treeplots_enrichment), function(x) {
  
  treeplots_enrichment[[x]]$GO + 
    ggtitle(gsub('\\.', ' ', x))
  
  # print(x)
  # print(assebly_plots(enrich = treeplots_enrichment[[x]], cls = x))
})
dev.off()

# x = "2:0.(RNA-Prot light)"

# visualising results - GSEA ----

gsea_res = readRDS('data/compensation/cs_gsea.rds')

go_gsea_RNA = extract_res(gsea_res, db = 'GO', omic = 'RNA')

sapply(go_gsea_RNA, function(x) {
  x %>% 
    # filter(p.adjust <= .05) %>% 
    dim
})

go_gsea_RNA = lapply(go_gsea_RNA %>% names, function(x) {
  go_gsea_RNA[[x]]@result %>% 
    mutate(group = x) %>% 
    tidyr::separate(group, into = c('karyotype', 'omic'), sep = '\\.')
}) %>% 
  bind_rows()

go_gsea_RNA %>%
  filter(p.adjust <= .05) %>% 
  ggplot(aes(
    y = Description, 
    x = NES, 
    fill = p.adjust
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ggsci::scale_fill_gsea() + 
  facet_wrap(~karyotype, scales = 'free')


go_gsea_prot = extract_res(gsea_res, db = 'GO', omic = 'protein')

sapply(go_gsea_prot, function(x) {
  x %>% 
    # filter(p.adjust <= .05) %>% 
    dim
})

go_gsea_prot = lapply(go_gsea_prot %>% names, function(x) {
  go_gsea_prot[[x]]@result %>% 
    mutate(group = x) %>% 
    tidyr::separate(group, into = c('karyotype', 'omic'), sep = '\\.')
}) %>% 
  bind_rows()

go_gsea_prot %>%
  filter(p.adjust <= .05, karyotype == '2:2') %>% 
  ggplot(aes(
    y = Description, 
    x = NES, 
    fill = p.adjust
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ggsci::scale_fill_gsea() + 
  facet_wrap(~karyotype, scales = 'free')


reac_gsea_prot = extract_res(gsea_res, db = 'REAC', omic = 'protein')
reac_gsea_prot = lapply(reac_gsea_prot %>% names, function(x) {
  reac_gsea_prot[[x]]@result %>% 
    mutate(group = x) %>% 
    tidyr::separate(group, into = c('karyotype', 'omic'), sep = '\\.')
}) %>% 
  bind_rows()

go_gsea_prot %>%
  filter(p.adjust <= .05) %>% 
  ggplot(aes(
    y = Description, 
    x = NES, 
    fill = p.adjust
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ggsci::scale_fill_gsea() + 
  facet_wrap(~karyotype, scales = 'free')

kegg_gsea_prot = extract_res(gsea_res, db = 'KEGG', omic = 'protein')
kegg_gsea_prot = lapply(kegg_gsea_prot %>% names, function(x) {
  kegg_gsea_prot[[x]]@result %>% 
    mutate(group = x) %>% 
    tidyr::separate(group, into = c('karyotype', 'omic'), sep = '\\.')
}) %>% 
  bind_rows()

kegg_gsea_prot %>%
  filter(p.adjust <= .05) %>% 
  ggplot(aes(
    y = Description, 
    x = NES, 
    fill = p.adjust
  )) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ggsci::scale_fill_gsea() + 
  facet_wrap(~karyotype, scales = 'free')

# testing on a msidb -- ribosomes/translation ------

cp_gmt = clusterProfiler::read.gmt('data/utilities/c2.cp.v2026.1.Hs.symbols.gmt')

ribo_genes = cp_gmt %>% 
  filter(term %in% c('KEGG_RIBOSOME', 'WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS')) %>% 
  mutate(term = as.character(term)) %>% 
  split(.$term)
  
rr_expr_groups = lapply(ribo_genes, function(x) {
  gg = x$gene
  
  df_groups %>% 
    filter(name %in% gg) %>% 
    mutate(term = unique(x$term))
  
}) %>% 
  bind_rows()
  
rr_expr_groups %>% 
  group_by(term, karyotype, reg_group) %>% 
  dplyr::count() %>% 
  ggplot(aes(
    y = n, 
    x = reg_group, 
    fill  =karyotype
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_brewer(palette = 'Pastel2') + 
  theme_bw() + 
  facet_wrap(~term)



# do not look ###############################################################
# # trying a different method to distinguish the groups -- t-test
# 
# df = df %>% 
#   select(-n) 
# 
# df_test = df %>% 
#   group_by(name) %>%
#   mutate(D_z = (CS - mean(CS)) / sd(CS),
#          close_to_zero = abs(D_z) < 1) 
# 
# # test = df %>% 
# #   split(interaction(.$omic, .$karyotype)) %>% 
# #   lapply(., function(d) {
# #     
# #     t.test(d$CS, mu = 0)    
# #     
# #   })
# #   
# # test$`RNA.1:0`
# 
# 
# # groups = setdiff(unique(df_groups$reg_group), "Intermediate/Other")
# groups = unique(df_groups$reg_group)
# genes_by_cluster = lapply(groups, function(g) {
#   df_groups %>% dplyr::filter(reg_group == g) %>% dplyr::pull(name)
# })
# names(genes_by_cluster) <- groups
# 
# 
# formula_res <- compareCluster(genes_by_cluster, 
#                               fun = "enrichGO", 
#                               OrgDb = org.Hs.eg.db, 
#                               keyType = "SYMBOL",
#                               ont = "ALL")
# 
# # dir.create("results/enrichment", recursive = T)
# # saveRDS(formula_res, "results/enrichment/cluster_comparison.rds")
# formula_res@compareClusterResult %>% view()
# # sheet_write(ss = sheet_url, sheet = "Cluster Enrichment")
# dotplot(formula_res, showCategory = 20)
# 
# formula_res@compareClusterResult %>% 
#   dplyr::group_by(Cluster,ONTOLOGY) %>% 
#   filter(p.adjust <= 0.05) %>% 
#   count()

