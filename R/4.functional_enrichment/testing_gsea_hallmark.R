library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/orfeo/scratch/cdslab/vgazziero/organoids_prj/')
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/jovoniR/constants.R')
source('organoids_analysis/R/plot_utils/dge_utils_and_plots.R')


df_groups = readRDS('data/compensation/cs_classification.rds')
genes_groups = df_groups %>% 
  split(.$reg_group)

m_df <- msigdbr(species = "Homo sapiens")

m_t2g <- m_df  %>% 
  dplyr::filter(gs_collection == "H")  %>% 
  dplyr::select(gs_name, ncbi_gene)

genes_by_group_ranked = lapply(genes_groups, function(x) {
  
  gene_df <- bitr(
    x$name %>% unique,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  dd = gene_df %>% 
    full_join(., x, by = join_by('SYMBOL' == 'name')) %>% 
    filter(!is.na(ENTREZID)) %>% 
    mutate(mean_cs = (RNA + Protein)/2)
  
  setNames(nm = dd$ENTREZID, object = dd$mean_cs) %>% 
    sort(decreasing = T)
  
})

formula_res = compareCluster(genes_by_group_ranked, 
                             fun = 'GSEA',
                             TERM2GENE = m_t2g)
dotplot(formula_res)

# running enrichment on each gene group for karyotype
enrichment_groups = lapply(genes_by_group_ranked, function(s) {
  
  clusterProfiler::enricher(gene = names(s), TERM2GENE = m_t2g) %>% 
    setReadable(., OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% 
    as.data.frame()
  
})

# trying enrichment in another way





