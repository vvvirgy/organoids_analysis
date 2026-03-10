library(GOSemSim)
library(igraph)

formula_res = reg_groups_comparison
TOP_K = 10

go_df <- formula_res@compareClusterResult %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::arrange(p.adjust) %>% 
  dplyr::slice_head(n = TOP_K) %>% 
  tidyr::separate(Cluster, into = c('karyotype', 'reg_group'), sep = '\\.', remove = T)

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
  group_by(karyotype, reg_group) %>%
  group_modify(~ make_mst_edges(.x)) %>%
  ungroup()

# ---- Base semantic plot (overlay by condition)
semantic_plot_base <- ggplot(plot_df, aes(dim1, dim2)) +
  geom_segment(
    data = edges_df,
    aes(x=x, y=y, xend=xend, yend=yend, color=reg_group),
    alpha = 0.25, linewidth = 0.6
  ) +
  geom_point(aes(color = reg_group, size = size), alpha = 0.75) +
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
  ) + 
  facet_wrap(~karyotype, scales = 'free')

n_terms_to_plot = TOP_K
label_df = plot_df %>% 
  dplyr::slice_sample(n = n_terms_to_plot)

# ---------------------------
# 4) Two final semantic plots
# ---------------------------

# A) Overlay (labels once per condition, across methods) — uses distinct labels w/out method
semantic_plot_overlay <- semantic_plot_base +
  ggrepel::geom_text_repel(
    data = label_df %>% dplyr::select(Description, karyotype, reg_group, dim1, dim2) %>% dplyr::distinct(),
    aes(label = Description),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  )

semantic_plot_overlay
