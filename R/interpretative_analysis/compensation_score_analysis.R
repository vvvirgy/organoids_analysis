rm(list = ls())
.libPaths()
library(tidyverse)
# library(tidyplots)
# library(googlesheets4)
library(dplyr)
library(boot)
library(clusterProfiler)
library(ReactomePA)

# BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
source('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/plot_utils/dge_utils_and_plots.R')

DATA_PATH = '/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_proj/results'
SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_karyo_all_organoids_filt.rds"
META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

MIN_SAMPLES = 1

# Get gene/karyotypes with at least two samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>% 
  dplyr::group_by(hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= MIN_SAMPLES) %>% 
  dplyr::rename(name = hgnc_symbol)

df_dna = readRDS(paste(DATA_PATH, "DNA_lfc.rds", sep = '/')) %>% dplyr::rename(DNA_lfc = lfc)

df = readRDS(paste(DATA_PATH, "sf_psinorm_stable_FALSE/lfc_prot_and_rna_bind.rds", sep = '/'))  

df = df %>% 
  dplyr::group_by(karyotype, name) %>% 
  dplyr::filter(n() == 2)

df = df %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

df_dna = df_dna %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

df %>% 
  dplyr::select(name, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise("n_distinct_karyotypes" = n()) %>% 
  dplyr::count(n_distinct_karyotypes) %>% 
  dplyr::rename(n_genes = n) 

df %>% 
  dplyr::select(name, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() %>% 
  dplyr::select(karyotype) %>% 
  dplyr::count(karyotype) 

df = df %>% 
  dplyr::filter(omic != "DNA")

df = df %>% 
  ungroup() %>% 
  # dplyr::left_join(df_dna %>% dplyr::select(!omic), by = join_by("name" == "name", "karyotype" == "karyotype")) %>% 
  dplyr::left_join(df_dna) %>% 
  dplyr::select(-n) 

df = df %>%
  # mutate(CS = DNA_lfc - lfc) #%>%
  dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc))

df_mean_CS = df %>% 
  dplyr::group_by(name, omic) %>% 
  dplyr::summarise(mean_CS = mean(CS)) %>% 
  tidyr::pivot_wider(values_from = mean_CS, names_from = omic)

df %>% 
  ggplot(aes(
    x = karyotype, 
    y = CS,
    fill = omic
  )) + 
  geom_violin()

df %>% 
  dplyr::select(karyotype, name, CS, omic, DNA_lfc) %>% 
  tidyr::pivot_wider(values_from = CS, names_from = omic) %>% 
  ggplot(mapping = aes(x = RNA, y = protein)) +
  geom_point() + 
  facet_wrap(~karyotype) + 
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red')

saveRDS(df, 'data/compensation/cs_results.rds')

# dividing the genes in different groups

karyos = c('1:0', '2:0', '2:1', '2:2')

df_groups = lapply(karyos, function(x) {
  define_th(df, q = .5, karyo = x)  
})

df_groups = df_groups %>% 
  bind_rows()

df_groups %>% 
  group_by(karyotype, reg_group) %>% 
  count() %>% 
  ggplot(aes(
    x = reg_group, 
    y = n, 
    fill = karyotype
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Pastel2')

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
  facet_wrap(~karyotype)

genes_by_group = df_groups %>% 
  dplyr::select(karyotype, reg_group, name) %>% 
  split(interaction(.$karyotype, .$reg_group))

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

enrichment_groups = readRDS('data/compensation/cs_groups_enrichment.rds') 

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
  count() %>% 
  ggplot(aes(
    x = reg_group, 
    y = n, 
    fill = karyotype
  )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Pastel2')

enrichment_groups_res %>% 
  # filter(reg_group == "(RNA-Prot light)") %>% 
  # filter(reg_group == '(RNA-Prot light)') %>% 
  filter(reg_group == '(RNA-heavy)') %>% 
  filter(source %in% c('GO', 'REAC')) %>% 
  # filter(source == 'GO') %>% 
  # filter(source == 'REAC') %>%
  # filter(source != 'REAC') %>%
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
enrichment_groups_res %>% 
  separate_rows(geneID, sep = '/') %>% 
  full_join(., df_groups, by = join_by(
    'geneID' == 'name', 
    'reg_group' == 'reg_group', 
    'karyotype' == 'karyotype'
  )) %>% 
  group_by(Description, karyotype, reg_group) %>% 
  mutate(mean_RNA = mean(RNA, na.rm = T), 
         mean_protein = mean(protein, na.rm = T)) %>% 
  filter(reg_group == '(RNA-prot heavy)') %>% 
  filter(source == 'KEGG') %>% 
  ggplot(aes(
    y = Description, 
    x = karyotype, 
    color = mean_RNA
  )) + 
  geom_point() + 
  scale_color_viridis_c() +
  theme_bw()
  




# trying a different method to distinguish the groups -- t-test

df = df %>% 
  select(-n) 

df_test = df %>% 
  group_by(name) %>%
  mutate(D_z = (CS - mean(CS)) / sd(CS),
         close_to_zero = abs(D_z) < 1) 

# test = df %>% 
#   split(interaction(.$omic, .$karyotype)) %>% 
#   lapply(., function(d) {
#     
#     t.test(d$CS, mu = 0)    
#     
#   })
#   
# test$`RNA.1:0`


# old way


category_colors <- c(
  "(RNA-prot heavy)" = "#AD002AB2",
  "(Prot-heavy)" = "#E18727B2",
  "(RNA-heavy)" = "#20854Eb2",
  "(RNA-Prot light)" = "#00468BB2",
  "Intermediate/Other" = "gainsboro"
)

reg_groups_plot = df_groups %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS")

# groups = setdiff(unique(df_groups$reg_group), "Intermediate/Other")
groups = unique(df_groups$reg_group)
genes_by_cluster = lapply(groups, function(g) {
  df_groups %>% dplyr::filter(reg_group == g) %>% dplyr::pull(name)
})
names(genes_by_cluster) <- groups


formula_res <- compareCluster(genes_by_cluster, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "SYMBOL",
                              ont = "ALL")

# dir.create("results/enrichment", recursive = T)
# saveRDS(formula_res, "results/enrichment/cluster_comparison.rds")
formula_res@compareClusterResult %>% view()
# sheet_write(ss = sheet_url, sheet = "Cluster Enrichment")
dotplot(formula_res, showCategory = 20)

formula_res@compareClusterResult %>% 
  dplyr::group_by(Cluster,ONTOLOGY) %>% 
  filter(p.adjust <= 0.05) %>% 
  count()

library(GOSemSim)
library(igraph)

TOP_K = 20

go_df <- formula_res@compareClusterResult %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::arrange(p.adjust) %>% 
  filter(ONTOLOGY == 'BP') %>% 
  dplyr::slice_head(n = TOP_K)

# ---- semantic similarity matrix (Wang)
hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
go_ids <- unique(go_df$ID)

S <- mgoSim(go_ids, go_ids, semData = hsGO, measure = "Wang", combine = NULL)
S[is.na(S)] <- 0
D <- as.dist(1 - S)

# ---- 2D embedding (MDS)
xy <- cmdscale(D, k = 2) %>% as.data.frame()
colnames(xy) <- c("dim1", "dim2")
xy$ID <- go_ids

plot_df <- go_df %>%
  left_join(xy, by = "ID") %>%
  mutate(size = -log10(p.adjust))

make_mst_edges <- function(df) {
  n <- nrow(df)
  if (n < 3) return(tibble(x=double(), y=double(), xend=double(), yend=double()))
  
  coords <- as.matrix(df[, c("dim1","dim2")])
  distmat <- as.matrix(dist(coords))
  
  g <- graph_from_adjacency_matrix(distmat, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  mst_g <- mst(g, weights = E(g)$weight)
  el <- as.data.frame(as_edgelist(mst_g))
  colnames(el) <- c("i","j")
  el$i <- as.integer(el$i); el$j <- as.integer(el$j)
  
  tibble(
    x    = df$dim1[el$i],
    y    = df$dim2[el$i],
    xend = df$dim1[el$j],
    yend = df$dim2[el$j]
  )
}

edges_df <- plot_df %>%
  group_by(Cluster) %>%
  group_modify(~ make_mst_edges(.x)) %>%
  ungroup()

# ---- Base semantic plot (overlay by condition)
semantic_plot_base <- ggplot(plot_df, aes(dim1, dim2)) +
  geom_segment(
    data = edges_df,
    aes(x=x, y=y, xend=xend, yend=yend, color=Cluster),
    alpha = 0.25, linewidth = 0.6
  ) +
  geom_point(aes(color = Cluster, size = size), alpha = 0.75) +
  #facet_grid(~Cluster) +
  #facet_wrap(~condition, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = category_colors) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    x = "GO semantic dimension 1",
    y = "GO semantic dimension 2",
    color = "Cluster",
    size  = expression(-log[10](adjP))
  )

n_terms_to_plot = TOP_K
label_df = plot_df %>% 
  dplyr::slice_sample(n = n_terms_to_plot)

# ---------------------------
# 4) Two final semantic plots
# ---------------------------

# A) Overlay (labels once per condition, across methods) — uses distinct labels w/out method
semantic_plot_overlay <- semantic_plot_base +
  ggrepel::geom_text_repel(
    data = label_df %>% dplyr::select(Description, Cluster, dim1, dim2) %>% dplyr::distinct(),
    aes(label = Description),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  )

ggsave("img/semantic_enrichment.pdf", plot = semantic_plot_overlay, width = 10, height = 8, units = "in")
ggsave("img/semantic_enrichment.png", plot = semantic_plot_overlay, width = 10, height = 8, units = "in", dpi = 450)