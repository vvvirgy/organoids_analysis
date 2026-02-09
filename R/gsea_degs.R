library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(enrichplot)
library(org.Hs.eg.db)

multi_omics = readRDS('data/lfc_prot_and_rna_bind.rds')

# volcano plot of expression
pth = .05
fth = .75

multi_omics = classify_genes(multi_omics)
multi_omics = multi_omics %>% 
  mutate(omic = factor(omic, levels = c('RNA', 'protein')))

degs_karyo_omic = multi_omics %>% 
  filter(!is.na(fc_cls)) %>% 
  # group_by(karyotype, omic) %>% 
  split(interaction(.$omic, .$karyotype, sep = "_"))

degs_karyo_omic_ranked = lapply(degs_karyo_omic, function(x) {
  
  x = x %>% 
    mutate(rank = pval*sign(lfc))
  
  setNames(nm = x$name,
           object = x$rank) %>%
    sort(., decreasing = T)
    
})

degs_karyo_omic_ranked_ids = lapply(degs_karyo_omic_ranked, function(x) {
  
  # convert
  gene_df <- bitr(
    names(x),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )  
  
  x <- x[gene_df$SYMBOL]
  names(x) <- gene_df$ENTREZID
  x <- sort(x, decreasing = TRUE)
  
})

gsea_res = lapply(degs_karyo_omic_ranked_ids, run_gsea)

plot_nes(gsea_res$`RNA_1:0`$GO)
plot_nes(gsea_res$`protein_1:0`$GO)
plot_nes(gsea_res$`protein_2:2`$GO)
plot_nes(gsea_res$`RNA_2:1`$GO)
plot_nes(gsea_res$`protein_2:0`$GO)



