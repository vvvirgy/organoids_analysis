rm(list=ls())
library(tidyverse)
library(dplyr)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

args = NULL

# Default values if not provided
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))

basepath = '/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data'
SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_karyo_all_organoids_filt.rds"
META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

IMG_PATH = paste0("img/sf_", sf_method, "_stable_", use_stable)
RES_PATH = paste0("data/compensation_score/sf_", sf_method, "_stable_", use_stable)
RNA_PATH = paste0("data/results/RNA/lfc_res_",sf_method,"_stable_",use_stable,".rds")
# gs4_auth(email = "santacatterinagiovanni@gmail.com")
# ss = gs4_create(paste0("organoids_", sf_method, "_", use_stable))

dir.create(IMG_PATH, recursive = T)
dir.create(RES_PATH, recursive = T)

# Parameters
#sheet_url <- "https://docs.google.com/spreadsheets/d/1uoiYWQg9EemHiriYU9EzFPBf1fF2KNQYsiWfurJnID0/"
alpha = .05
MIN_SAMPLES = 1

# Get gene/karyotypes with at least two samples
karyotypes_df_all = readRDS(META_PATH)
karyotypes_df_good = karyotypes_df_all %>% 
  dplyr::group_by(hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= MIN_SAMPLES) %>% 
  dplyr::rename(name = hgnc_symbol)

res_rna = readRDS(RNA_PATH) %>% dplyr::mutate(omic = "RNA") %>% 
  dplyr::mutate(lfc = ifelse(abs(lfc) > 10, sign(lfc) * 10, lfc))

res_prot = readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/fc_tb_clean_v2.rds") %>% 
  dplyr::rename(coef = condition, name = PG.Genes, lfc = diff, pval = p.val, adj_pval = p.adj) %>% 
  dplyr::mutate(omic = "protein") %>% 
  tidyr::separate(coef, sep = "X", into = c(".", "karyotype")) %>% 
  dplyr::select(lfc, name, karyotype, omic) %>% 
  mutate(karyotype = str_replace(karyotype, "\\.", ":")) %>% 
  dplyr::group_by(karyotype)

df_dna = readRDS("data/results/DNA_lfc.rds") %>% dplyr::rename(DNA_lfc = lfc)

df = dplyr::bind_rows(res_rna, res_prot) %>% na.omit()
df$omic = factor(df$omic, levels = c("RNA", "protein"))

df %>% saveRDS(file.path(RES_PATH, "lfc_prot_and_rna_bind.rds"))

df = readRDS(file.path(RES_PATH, "lfc_prot_and_rna_bind.rds"))

df = df %>% 
  dplyr::group_by(karyotype, name) %>% 
  dplyr::filter(n() == 2)

df = df %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

df_dna = df_dna %>% 
  dplyr::left_join(karyotypes_df_good) %>% 
  na.omit()

# df %>% dplyr::filter(name == "AASS")
genes_with_4_karyotypes = df %>%
  dplyr::select(name, karyotype) %>%
  dplyr::distinct() %>%
  dplyr::group_by(name) %>%
  dplyr::summarise("n_distinct_karyotypes" = n()) %>% 
  dplyr::filter(n_distinct_karyotypes == 4) %>% 
  dplyr::pull(name) %>% unique()

df %>%
  dplyr::select(name, karyotype) %>%
  dplyr::distinct() %>%
  dplyr::group_by(name) %>%
  dplyr::summarise("n_distinct_karyotypes" = n()) %>%
  dplyr::count(n_distinct_karyotypes) %>%
  dplyr::rename(n_genes = n) %>% 
  saveRDS(file.path(RES_PATH, "karyotypes_distribution.rds"))

# corum_genes = read.delim("data/CORUM_gene_list.txt") %>% dplyr::pull(GeneSym)
# df = df %>% dplyr::mutate(group = ifelse(name %in% corum_genes, "complex", "non-complex"))
# df %>% ungroup() %>% dplyr::select(name, group) %>% dplyr::distinct() %>% dplyr::pull(group) %>% table()

df = df %>% dplyr::filter(omic != "DNA")

# df = df %>% 
#   tidyr::separate(karyotype, into = c("A", "B"), sep = ":", remove = FALSE, convert = TRUE) %>% 
#   dplyr::mutate(copies = A + B, DNA_lfc = log2(copies / 2)) %>% 
#   dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc))

df = df %>% 
  ungroup() %>% 
  # dplyr::left_join(df_dna %>% dplyr::select(!omic), by = join_by("name" == "name", "karyotype" == "karyotype")) %>% 
  dplyr::left_join(df_dna) %>% 
  dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc)) 

# computing a mean CS

df_mean_CS = df %>% 
  dplyr::group_by(name, omic) %>% 
  dplyr::summarise(mean_CS = mean(CS)) %>% 
  tidyr::pivot_wider(values_from = mean_CS, names_from = omic)

df_mean_CS %>% 
  ggplot(mapping = aes(x = RNA, y = protein)) +
  geom_point()

q = .5

rna_cuts = df_mean_CS %>%
  dplyr::mutate(s = sign(RNA)) %>%
  dplyr::group_by(s) %>%
  dplyr::summarise(q = quantile(RNA, q)) %>%
  dplyr::pull(q) %>%
  sort()

prot_cuts = df_mean_CS %>%
  dplyr::mutate(s = sign(protein)) %>%
  dplyr::group_by(s) %>%
  dplyr::summarise(q = quantile(protein, q)) %>%
  dplyr::pull(q) %>%
  sort()

# 3. Assign the genes to the specific groups from your text
df_groups <- df_mean_CS %>%
  mutate(reg_group = case_when(
    # Group 1: High RNA CS (>65th) and Low Protein CS (<35th)
    RNA > rna_cuts[2] & protein < prot_cuts[1] ~ "(RNA-heavy)",
    RNA > rna_cuts[2] & protein > prot_cuts[2] ~ "(RNA-prot heavy)",
    
    # Group 2: Low RNA CS (<35th) and High Protein CS (>65th)
    RNA < rna_cuts[1] & protein > prot_cuts[2] ~ "(Prot-heavy)",
    RNA < rna_cuts[1] & protein < prot_cuts[1] ~ "(RNA-Prot light)",
    
    TRUE ~ "Intermediate/Other"
  ))
saveRDS(df_groups, 'data/compensation_score/sf_psinorm_stable_FALSE/CS_tables/cs_mean_groups.rds')

category_colors <- c(
  "(RNA-prot heavy)" = "#AD002AB2",
  "(Prot-heavy)" = "#E18727B2",
  "(RNA-heavy)" = "#20854Eb2",
  "(RNA-Prot light)" = "#00468BB2",
  "Intermediate/Other" = "gainsboro"
)

# df_groups %>% 
#   dplyr::group_by(reg_group) %>% 
#   dplyr::count(reg_group) %>% 
#   sheet_write(ss = ss, sheet = "Regulation groups distribution")

# df %>% 
#   dplyr::filter(name == "ENPP4")
# 
# df %>% 
#   dplyr::filter(name %in% c("APOBEC3C","CDA","ENPP4","GDA","LCMT2","MACROD1"))
# df_groups %>% 
#   dplyr::filter(name %in% c("APOBEC3C","CDA","ENPP4","GDA","LCMT2","MACROD1"))

reg_groups_plot = df_groups %>% 
  ggplot(mapping = aes(x = RNA, y = protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS")

reg_groups_plot

ggsave(plot = reg_groups_plot, filename = file.path(IMG_PATH, "reg_groups_plot.pdf"), width = 7, height = 5)
ggsave(plot = reg_groups_plot, filename = file.path(IMG_PATH, "reg_groups_plot.png"), width = 7, height = 5, dpi = 450, units = "in")

groups = setdiff(unique(df_groups$reg_group), "Intermediate/Other")
genes_by_cluster = lapply(groups, function(g) {
  df_groups %>% dplyr::filter(reg_group == g) %>% dplyr::pull(name)
})
names(genes_by_cluster) <- groups

sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-prot heavy)`) / length(genes_by_cluster$`(RNA-prot heavy)`)
sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-Prot light)`) / length(genes_by_cluster$`(RNA-Prot light)`)
sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-prot heavy)`) / length(genes_with_4_karyotypes)

formula_res <- compareCluster(genes_by_cluster, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "SYMBOL",
                              ont = "BP")


dir.create(file.path(RES_PATH, "enrichment"), recursive = T)
saveRDS(formula_res, file.path(RES_PATH, "enrichment", "cluster_comparison.rds"))
# formula_res@compareClusterResult %>%
#   sheet_write(ss = ss, sheet = "Cluster Enrichment")
dotplot(formula_res, showCategory = 20)

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


res_clustering %>% 
  ungroup() %>% 
  dplyr::filter(Cluster == '(RNA-Prot light)') %>% 
  arrange(desc(prot_cs)) %>% 
  slice_head(n = 30) %>% 
  # slice_max(prot_cs, n = 30) %>% 
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
  ggsci::scale_color_gsea()




library(GOSemSim)
library(igraph)

TOP_K = 10

go_df <- res_clustering %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::arrange(desc(RNA_cs)) %>% 
  dplyr::slice_head(n = TOP_K)

row = go_df %>% 
  dplyr::filter(Cluster == "(RNA-Prot light)") %>% 
  dplyr::slice_head(n = 1)

classify_term_by_genetics <- function(gene_string, ref_df) {
  # Split the geneID string into individual genes
  genes <- unlist(str_split(gene_string, "/"))
  
  # Look up classes for each gene
  gene_classes <- sapply(genes, function(g) {
    # Get karyotypes for this gene from your reference df
    karyotypes <- ref_df %>% 
      filter(name == g) %>% 
      pull(karyotype) %>% 
      unique()
    
    if (length(karyotypes) == 0) return("Unknown")
    
    # Your original logic refined
    if (length(karyotypes) == 4) {
      return("all")
    } else if (all(karyotypes %in% c("2:1", "2:2"))) {
      return("Gain")
    } else if (all(karyotypes %in% c("2:1", "2:2", "2:0"))) {
      return("Gain/Neutral")
    } else if (all(karyotypes %in% c("1:0", "2:0"))) {
      return("Loss/Neutral")
    } else if (all(karyotypes %in% c("1:0"))) {
      return("Loss")
    } else if (any(karyotypes %in% c("1:0")) & any(karyotypes %in% c("2:1", "2:2"))) {
      return("Gain/Loss")
    } else {
      return("Neutral")
    }
  })
  
  # 2. Summarize the whole term (The Assessment)
  class_counts <- table(gene_classes)
  
  # Logic to determine the "Main" driver
  if (length(class_counts) == 0) return("No Data")
  
  # Calculate proportions (excluding Neutral if desired)
  total <- sum(class_counts)
  gain_sum <- sum(class_counts[c("Gain", "Gain/Neutral")], na.rm = TRUE)
  loss_sum <- sum(class_counts[c("Loss", "Loss/Neutral")], na.rm = TRUE)
  
  if (gain_sum > (total * 0.5)) {
    return("Mainly Gain")
  } else if (loss_sum > (total * 0.5)) {
    return("Mainly Loss")
  } else if (gain_sum > 0 & loss_sum > 0) {
    return("Mixed (Both)")
  } else {
    return("Balanced/Neutral")
  }
}

go_df_classified <- go_df %>%
  rowwise() %>%
  mutate(
    term_classification = classify_term_by_genetics(geneID, df)
  ) %>%
  ungroup()

print(table(go_df_classified$term_classification))

# definying which pathway is more or less impacted by the cna

go_df_classified %>% 
  dplyr::select(Cluster, Description, RNA_cs, prot_cs, geneID, term_classification) %>% 
  group_by(Cluster) %>% 
  dplyr::arrange(desc(RNA_cs)) %>% 
  dplyr::slice_head(n = TOP_K) %>% 
  pivot_longer(ends_with('cs'), names_to = 'omic', values_to = 'CS') %>% 
  ggplot(aes(
    y = Description, 
    x = CS, 
    fill = term_classification
  )) + 
  geom_bar(stat = 'identity') + 
  facet_grid(Cluster~omic, scales = 'free' ) + 
  theme_bw()



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

semantic_plot_overlay

ggsave(file.path(IMG_PATH, "semantic_enrichment.pdf"), plot = semantic_plot_overlay, width = 10, height = 8, units = "in")
ggsave(file.path(IMG_PATH, "semantic_enrichment.png"), plot = semantic_plot_overlay, width = 10, height = 8, units = "in", dpi = 450)

# Classify based on DNA to RNA reg and with RNA to prot reg
# df_RNA_to_DNA = df %>% 
#   dplyr::filter(omic == "RNA") %>% 
#   dplyr::mutate(CS_dna_rna = (DNA_lfc - lfc) * sign(DNA_lfc))
# 
# df_prot_to_RNA = df %>% 
#   dplyr::select(omic, karyotype, name, group, lfc) %>% 
#   tidyr::pivot_wider(names_from = omic, values_from = lfc) %>% 
#   dplyr::mutate(CS_rna_prot = (RNA - protein) * sign(RNA))
# 
# df_regulation = df_RNA_to_DNA %>% 
#   dplyr::left_join(df_prot_to_RNA) %>% 
#   dplyr::select(karyotype, name, group, CS_dna_rna, CS_rna_prot)
# 
# df_regulation %>% 
#   ggplot(mapping = aes(x = CS_dna_rna, y = CS_rna_prot)) +
#   geom_point() +
#   facet_wrap(~karyotype)

