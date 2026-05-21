
rm(list = ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(analogue)
library(dendextend)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

df = readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/CS_scores_prot_and_rna.rds") %>% 
  dplyr::filter(!karyotype %in% c("2:0"))

# --- 1. Reshape ---
wide <- df %>%
  dplyr::select(name, omic, karyotype, CS) %>%
  pivot_wider(names_from = c(omic, karyotype), values_from = CS, names_sep = "_") %>%
  tibble::column_to_rownames("name")

wide <- wide[complete.cases(wide), ]
mat  <- as.matrix(wide)

# --- 2. Row z-score, then drop any zero-variance rows ---
mat_scaled <- t(scale(t(mat)))

mat_scaled_plot <- t(apply(mat, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)   # leave constant rows alone
  x / s
}))
rownames(mat_scaled_plot) <- rownames(mat)
mat_scaled_plot <- mat_scaled_plot[complete.cases(mat_scaled), ]

# --- 3. Cluster ---
d  <- dist(mat_scaled, method = "euclidean")
hc <- hclust(d, method = "complete")

k <- 4
clusters <- cutree(hc, k = k)
gene_clusters <- tibble::tibble(name = names(clusters), cluster = clusters)
print(table(gene_clusters$cluster))

# --- 4. Annotations ---
annot_row <- data.frame(cluster = factor(clusters))
rownames(annot_row) <- names(clusters)

col_meta <- do.call(rbind, strsplit(colnames(mat_scaled), "_", fixed = TRUE))
annot_col <- data.frame(omic = col_meta[, 1], karyotype = col_meta[, 2])
rownames(annot_col) <- colnames(mat_scaled)

# Robust cluster palette
n_clust <- length(unique(clusters))
cluster_palette <- if (n_clust <= 8) {
  brewer.pal(max(3, n_clust), "Set2")[seq_len(n_clust)]
} else {
  colorRampPalette(brewer.pal(8, "Set2"))(n_clust)
}

annot_colors <- list(
  cluster   = setNames(cluster_palette, levels(annot_row$cluster)),
  omic      = setNames(brewer.pal(max(3, length(unique(annot_col$omic))), "Dark2"),
                       unique(annot_col$omic)),
  karyotype = setNames(brewer.pal(max(3, length(unique(annot_col$karyotype))), "Pastel1"),
                       unique(annot_col$karyotype))
)

# --- 5. Column ordering ---
karyo_levels <- c("1:0", "2:0", "2:1", "2:2")
col_order <- order(annot_col$omic,
                   factor(annot_col$karyotype, levels = karyo_levels))
mat_scaled_plot <- mat_scaled_plot[, col_order]
annot_col  <- annot_col[col_order, , drop = FALSE]

# --- 6. Breaks based on SCALED matrix (the one being plotted) ---
# breaks <- seq(-max(abs(mat_scaled_plot), na.rm = TRUE),
#               max(abs(mat_scaled_plot), na.rm = TRUE),
#               length.out = 101)
c <- 0.0
half_range <- max(abs(range(mat_scaled, na.rm = TRUE) - c))
breaks <- seq(c - half_range, c + half_range, length.out = 101)
cs_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Optional: visual gaps between omics
gaps_col <- cumsum(rle(as.vector(annot_col$omic))$lengths)
gaps_col <- gaps_col[-length(gaps_col)]

# --- 7. Heatmap ---
pheatmap(
  mat_scaled_plot, 
  color             = cs_colors,
  breaks            = breaks,
  cluster_rows      = hc,
  cluster_cols      = FALSE,
  cutree_rows       = k,
  gaps_col          = gaps_col,
  annotation_row    = annot_row,
  annotation_col    = annot_col,
  annotation_colors = annot_colors,
  show_rownames     = nrow(mat_scaled) < 80,
  show_colnames     = TRUE,
  fontsize_col      = 8,
  main              = "CS (row z-score) by omic × karyotype"
)

library(clusterProfiler)
library(org.Hs.eg.db)   # human; swap for org.Mm.eg.db etc. if needed
library(dplyr)
library(tidyr)
library(purrr)

# --- 1. Convert gene symbols to ENTREZ IDs (required by enrichKEGG) ---
gene_clusters <- gene_clusters %>%
  mutate(entrez = mapIds(org.Hs.eg.db,
                         keys      = name,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")) %>%
  filter(!is.na(entrez))

# Universe = all genes that went into clustering (background for enrichment)
universe_entrez <- unique(gene_clusters$entrez)

# --- 2. KEGG enrichment per cluster ---
kegg_per_cluster <- gene_clusters %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(entrez)), .groups = "drop") %>%
  mutate(kegg = map(genes, ~ enrichKEGG(
    gene          = .x,
    universe      = universe_entrez,
    organism      = "hsa",         # 'mmu' for mouse, etc.
    keyType       = "kegg",
    pvalueCutoff  = 1,             # keep everything; we'll filter on FDR
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )))

kegg_results <- kegg_per_cluster %>%
  mutate(df = map(kegg, ~ if (!is.null(.x)) as.data.frame(.x) else NULL)) %>%
  dplyr::select(cluster, df) %>%
  unnest(df) %>%
  filter(p.adjust < 0.1)

print(kegg_results)

diploid_noise = readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise.RDS") %>% 
  dplyr::filter(omic == "Protein")

k_to_plot = 1
df %>% 
  dplyr::filter(karyotype != "1:0") %>% 
  dplyr::filter(name %in% (gene_clusters %>% dplyr::filter(cluster == k_to_plot) %>% dplyr::pull(name))) %>% 
  ggplot(mapping = aes(x = CS)) +
  geom_histogram(binwidth = .1) +
  facet_grid(omic~karyotype) +
  geom_vline(xintercept = diploid_noise$mean - 3 * diploid_noise$sigma) +
  geom_vline(xintercept = diploid_noise$mean + 3 * diploid_noise$sigma) +
  theme_light()
