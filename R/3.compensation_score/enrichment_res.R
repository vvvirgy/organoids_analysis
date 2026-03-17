library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(msigdbr)
META_PATH = "data/karyotypes_genes_filtered_scrna.rds"

MIN_SAMPLES = 1

# Get gene/karyotypes with at least two samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>% 
  dplyr::group_by(hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= MIN_SAMPLES) %>% 
  dplyr::rename(name = hgnc_symbol)

# load data
df_groups = readRDS('data/compensation_score/sf_psinorm_stable_FALSE/CS_tables/cs_mean_groups.rds')

df_groups %>% 
  filter(reg_group == '(RNA-Prot light)') %>% 
  pull(name)

db = msigdbr(species = "Homo sapiens")
db_c2 = db %>% 
  filter(gs_collection == 'C2') %>% 
  dplyr::select(gs_name, ncbi_gene)

genes_by_group = df_groups %>% 
  dplyr::select(reg_group, name) %>% 
  split(.$reg_group)

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


pathway_enrichment = lapply(genes_by_group[1:4], function(s) {
  enricher(s$ENTREZID, TERM2GENE = db_c2) %>% 
    setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID")
})

dotplot(pathway_enrichment$`(RNA-Prot light)`)

# enrichment with comparecluster
ggenes = lapply(genes_by_group, function(x) {
  x$ENTREZID
})

formula_res <- compareCluster(ggenes, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "ENTREZID",
                              ont = "ALL") %>% 
  setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
# dotplot(formula_res, showCategory = 15)

res_clustering = formula_res@compareClusterResult %>% 
  separate_rows(geneID, sep = '/', convert = T) %>% 
  left_join(., df_groups, by = join_by('geneID' == 'name', 
                                       'Cluster' == 'reg_group')) %>% 
  group_by(Cluster, Description) %>% 
  mutate(RNA_cs = mean(RNA), 
         prot_cs = mean(protein), 
         geneID = paste(geneID, collapse = '/')) %>% 
  dplyr::select(-c(RNA, protein)) %>% 
  distinct() 

res_clustering$Cluster %>% table

plot_res = function(res, 
                    cluster,
                    variable,
                    n_terms) {
  if(variable == 'p.adjust') {
    
    df = res %>% 
      ungroup() %>% 
      dplyr::filter(Cluster == cluster) %>% 
      slice_min(.data[[variable]], n = n_terms)
    
  } else {
    
    df = res %>% 
      ungroup() %>% 
      dplyr::filter(Cluster == cluster) %>% 
      # arrange(desc(prot_cs), desc(RNA_cs)) %>% 
      # arrange(p.adjust) %>% 
      # slice_head(n = n_terms) %>% 
      slice_max(.data[[variable]], n = n_terms)
    
  }
  
  df %>% 
    separate(GeneRatio, into = c('A', 'B'), sep = '/', convert = T) %>% 
    mutate(GeneRatio = A/B, 
           A = NULL, 
           B = NULL) %>% 
    ggplot(
      aes(
        x = GeneRatio, 
        y = reorder(Description, +GeneRatio), 
        color = p.adjust
      )
    ) + 
    geom_point() + 
    theme_bw() + 
    labs(y = 'Biological term') + 
    ggsci::scale_color_gsea() +
    theme(axis.text.y = element_text(size = 10))

}

plot_res(res = res_clustering, cluster = '(RNA-Prot light)', variable = 'p.adjust', n_terms = 30) + 
  ggtitle('RNA-prot heavy associated terms, ranked by p-adjust')

# pdf('res/compensation_score/annotations_rna_prot_heavy.pdf', width = 9, height = 10)
# plot_res(res_clustering, cluster = '(RNA-prot heavy)', variable = 'RNA_cs', n_terms = 30) + 
#   ggtitle('RNA-prot heavy associated terms, ranked by CS RNA')
# plot_res(res = res_clustering, 
#          cluster = '(RNA-prot heavy)', 
#          variable = 'prot_cs', 
#          n_terms = 30) + 
#   ggtitle('RNA-prot heavy associated terms, ranked by CS protein')
# plot_res(res_clustering, cluster = '(RNA-prot heavy)', variable = 'p.adjust', n_terms = 30) + 
#   ggtitle('RNA-prot heavy associated terms, ranked by p-adjust')
# dev.off()

plot_res(res_clustering, cluster = '(RNA-prot heavy)', variable = 'RNA_cs', n_terms = 30) + 
  ggtitle('RNA-prot heavy associated terms, ranked by CS RNA')
ggsave('res/compensation_score/annotations_rna_prot_heavy_rank_rna.png', width = 11, height = 10)

plot_res(res_clustering, cluster = '(RNA-prot heavy)', variable = 'prot_cs', n_terms = 30) + 
  ggtitle('RNA-prot heavy associated terms, ranked by CS Protein')
ggsave('res/compensation_score/annotations_rna_prot_heavy_rank_prot.png', width = 9, height = 10)


pdf('res/compensation_score/annotations_prot_heavy.pdf', width = 9, height = 10)
plot_res(res_clustering, cluster = '(Prot-heavy)', variable = 'p.adjust', n_terms = 30) + 
  ggtitle('RNA-prot heavy associated terms, ranked by p-adjust')
dev.off()

pdf('res/compensation_score/annotations_rna_prot_light.pdf', width = 9, height = 10)
plot_res(res_clustering, cluster = '(RNA-Prot light)', variable = 'p.adjust', n_terms = 30) + 
  ggtitle('RNA-prot heavy associated terms, ranked by p-adjust')
dev.off()


plot_res(res_clustering, cluster = "Intermediate/Other", variable = 'p.adjust', n_terms = 30) + 
  ggtitle('Intermediate terms')
ggsave('res/compensation_score/annotations_intermediate.png', width = 9, height = 10)


cs = '(RNA-prot heavy)'
res_clustering %>% 
  ungroup() %>% 
  filter(Cluster == cs) %>% 
  slice_max(RNA_cs, n = 30) %>% 
  dplyr::select(Description, geneID, p.adjust) %>% 
  separate_rows(geneID, sep = '/') %>% 
  right_join(., karyotypes_df_good, by = join_by('geneID' == 'name')) %>% 
  filter(!is.na(Description)) %>% 
  ggplot(aes(
    reorder(Description, +p.adjust), 
    fill = karyotype
  )) + 
  geom_bar(stat = 'count', position = 'dodge')+
  coord_flip()

cs = '(RNA-prot heavy)'
res_clustering %>% 
  ungroup() %>% 
  filter(Cluster == cs) %>% 
  slice_max(prot_cs, n = 30) %>% 
  dplyr::select(Description, geneID, p.adjust) %>% 
  separate_rows(geneID, sep = '/') %>% 
  right_join(., karyotypes_df_good, by = join_by('geneID' == 'name')) %>% 
  filter(!is.na(Description)) %>% 
  ggplot(aes(
    reorder(Description, +p.adjust), 
    fill = karyotype
  )) + 
  geom_bar(stat = 'count', position = 'dodge')+
  coord_flip()
  

# visualise per karyotype results 

enrichment_groups = readRDS('data/compensation/cs_groups_enrichment_GO.rds')
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

# upset to see if there is anything in common

enrich_list = enrichment_groups_res %>%
  dplyr::select(Description, reg_group, karyotype) %>% 
  split(interaction(.$karyotype, .$reg_group))
enrich_list = lapply(enrich_list, function(x) {
  x$Description %>% unique
})

UpSetR::upset(fromList(enrich_list), nsets = 16, order.by = 'freq')

both_light_top_terms = enrichment_groups_res %>% 
  filter(reg_group == '(RNA-Prot light)') %>% 
  slice_min(p.adjust, by = karyotype, n = 10) %>% 
  pull(Description)

enrichment_groups_res %>% 
  filter(reg_group == '(RNA-Prot light)') %>% 
  filter(Description %in% both_light_top_terms) %>% 
  ggplot(aes(
    y = Description, 
    x = karyotype, 
    color = p.adjust
)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw()
ggsave('res/compensation_score/top10_by_karyo_both_light.pdf', width = 10, height = 10)

both_heavy_top_terms = enrichment_groups_res %>% 
  filter(reg_group == '(RNA-prot heavy)') %>% 
  slice_min(p.adjust, by = karyotype, n = 10) %>% 
  pull(Description)

enrichment_groups_res %>% 
  filter(reg_group == '(RNA-prot heavy)') %>% 
  filter(Description %in% both_heavy_top_terms) %>% 
  ggplot(aes(
    y = Description, 
    x = karyotype, 
    color = p.adjust
  )) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw()
ggsave('res/compensation_score/top10_by_karyo_both_light.pdf', width = 10, height = 10)
ggsave('res/compensation_score/top10_by_karyo_both_light.png', width = 10, height = 10)

rna_heavy_top_terms = enrichment_groups_res %>% 
  filter(reg_group == '(RNA-heavy)') %>% 
  slice_min(p.adjust, by = karyotype, n = 10) %>% 
  pull(Description)

enrichment_groups_res %>% 
  filter(reg_group == '(RNA-heavy)') %>% 
  filter(Description %in% rna_heavy_top_terms) %>% 
  ggplot(aes(
    y = Description, 
    x = karyotype, 
    color = p.adjust
  )) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw()
ggsave('res/compensation_score/top10_by_karyo_both_light.pdf', width = 10, height = 10)



  
  


# res_clustering %>% 
#   ungroup() %>% 
#   dplyr::filter(Cluster == '(RNA-prot heavy)') %>% 
#   # arrange(desc(prot_cs), desc(RNA_cs)) %>% 
#   arrange(p.adjust) %>% 
#   slice_head(n = 30) %>% 
#   # slice_max(prot_cs, n = 30) %>% 
#   separate(GeneRatio, into = c('A', 'B'), sep = '/', convert = T) %>% 
#   mutate(GeneRatio = A/B, 
#          A = NULL, 
#          B = NULL) %>% 
#   ggplot(
#     aes(
#       x = GeneRatio, 
#       y = reorder(Description, +GeneRatio), 
#       color = p.adjust
#     )
#   ) + 
#   geom_point() + 
#   theme_bw() + 
#   labs(y = 'Biological term') + 
#   ggsci::scale_color_gsea()






