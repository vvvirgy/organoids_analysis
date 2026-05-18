library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/constants.R')
df_groups = readRDS('data/compensation/cs_classification_gene_labels.rds')
genes_groups = df_groups %>% 
  split(.$reg_group)

genes_groups = lapply(genes_groups, function(x) {
  x$name %>% unique
})

# running enrichment on all genes
analysis_functions_GO = compareCluster(genes_groups, 
                                       fun = 'enrichGO', 
                                       keyType = 'SYMBOL', 
                                       ont = 'ALL', 
                                       OrgDb = org.Hs.eg.db)

saveRDS(analysis_functions_GO, 'data/enrichment_results/ORA_GO_categories_genes.rds')
write_csv(analysis_functions_GO@compareClusterResult, file = 'data/enrichment_results/ORA_GO_categories_genes.csv', col_names = T)

dotplot(analysis_functions_GO)

# running enrichment on specific pathways
m_df <- msigdbr(db_species = 'HS', species = 'human')

m_df$gs_collection %>% unique
collection_db = c('H', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C9')

# convert gene names to entrez
convert_names = function(x, from, to, org) {
  
  # convert
  gene_df <- bitr(
    x,
    fromType = from,
    toType   = to,
    OrgDb    = org
  )  
  
  x <- gene_df[,to]
  
  return(x)
}

genes_groups_entrez = lapply(genes_groups, function(x) {
  convert_names(x, from = 'SYMBOL', 'ENTREZID', org = 'org.Hs.eg.db')
})



enrichment_msigbd = lapply(collection_db, function(M) {
  m_t2g <- m_df  %>% 
    dplyr::filter(gs_collection == M)  %>% 
    dplyr::select(gs_name, ncbi_gene)
  
  lapply(genes_groups_entrez, function(gg) {
    enricher(gg, TERM2GENE = m_t2g)
  }) #%>%
    # setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  # compareCluster(genes_groups, 
  #                fun = 'enricher', 
  #                TERM2GENE = m_t2g, 
  #                # keyType = 'ENTREZID', 
  #                # OrgDb = org.Hs.eg.db
  #                ) 
  # lapply
  # enricher(rank_padj_entrez, TERM2GENE = m_t2g, pAdjustMethod = 'BH') %>% 
    # setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID")
})

names(enrichment_msigbd) = collection_db

# compare clusters
enrichment_msigbd_comparison = lapply(collection_db, function(M) {
  m_t2g <- m_df  %>% 
    dplyr::filter(gs_collection == M)  %>% 
    dplyr::select(gs_name, ncbi_gene)

  compareCluster(genes_groups_entrez,
                 fun = 'enricher',
                 TERM2GENE = m_t2g
                 ) 
 
})
names(enrichment_msigbd_comparison) = collection_db
saveRDS(enrichment_msigbd_comparison, 'data/enrichment_results/ORA_msigndb_categories_genes.rds')

dotplot(enrichment_msigbd_comparison$H)
dotplot(enrichment_msigbd_comparison$C1)
dotplot(enrichment_msigbd_comparison$C2)
dotplot(enrichment_msigbd_comparison$C3)
dotplot(enrichment_msigbd_comparison$C4)
dotplot(enrichment_msigbd_comparison$C6)

hallmark_results = enrichment_msigbd_comparison$H@compareClusterResult

chr_position_results = enrichment_msigbd_comparison$C1@compareClusterResult

pathway_results = enrichment_msigbd_comparison$C2@compareClusterResult
terms_pathway_best = pathway_results %>% 
  group_by(Cluster) %>% 
  slice_min(p.adjust, n = 20) %>% 
  pull(Description) %>% 
  unique
pathway_results_v2 = pathway_results %>% 
  filter(Description %in% terms_pathway_best)

oncogenic_sign_results = enrichment_msigbd_comparison$C6@compareClusterResult

  
# plotting function
plot_enrichment_results = function(df, cols, pth = .05) {
  df %>% 
    filter(p.adjust <= pth) %>% 
    tidyr::separate(GeneRatio, into = c('A', 'B'), sep = '/', convert = T, remove = T) %>% 
    mutate(GeneRatio = A/B) %>% 
    ggplot(aes(
      y = Description, 
      x = RichFactor, 
      fill = Cluster
    )) + 
    geom_bar(stat = 'identity', alpha = 1) + 
    facet_wrap(~Cluster, scales = 'free_x') + 
    theme_bw() + 
    scale_fill_manual(values = cols)
}

plot_enrichment_results(hallmark_results, cols = category_colors)
p2 = plot_enrichment_results(chr_position_results, cols = category_colors)
p2@facet$params$nrow = 1
p3 = plot_enrichment_results(pathway_results_v2, cols = category_colors) 
p3@facet$params$nrow = 1
p3
plot_enrichment_results(oncogenic_sign_results, cols = category_colors)

df = lfc_prot_and_rna_bind %>% filter(name %in% c('MYC', 'PIK3CA', 'ELL', 'EGFR')) 
p1 = df %>% 
  mutate(sig = ifelse(adj_pval <= 0.05, 'significant', 'ns')) %>% 
  ggplot(aes(
    x = karyotype, 
    y = lfc, 
    fill = omic,
    color = sig
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_color_manual(values = c('significant' = 'black', 'ns' = NA)) +
  facet_wrap(~name) + 
  theme_bw() 

dotplot(enrichment_msigbd_comparison$C9) + p1

# plot distribution of genes in the results
df_groups %>% 
  mutate(class_gene = ifelse(label != '', class, NA)) %>% 
  filter(!is.na(class_gene)) %>% 
  ggplot()