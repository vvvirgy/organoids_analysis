
library(tidyverse)
library(googlesheets4)
library(dplyr)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

# Default values if not provided
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))
source("utils.R")
source("constants.R")

IMG_PATH = paste0("img/sf_", sf_method, "_stable_", use_stable)
RES_PATH = paste0("results/sf_", sf_method, "_stable_", use_stable)
RNA_PATH = paste0("results/RNA/lfc_res_",sf_method,"_stable_",use_stable,".rds")
# gs4_auth(email = "santacatterinagiovanni@gmail.com")
# ss = gs4_create(paste0("organoids_", sf_method, "_", use_stable))

dir.create(IMG_PATH, recursive = T)
dir.create(RES_PATH, recursive = T)

# Parameters
#sheet_url <- "https://docs.google.com/spreadsheets/d/1uoiYWQg9EemHiriYU9EzFPBf1fF2KNQYsiWfurJnID0/"
alpha = .05
use_only_common_genes = TRUE
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
  dplyr::mutate(omic = "Protein") %>% 
  tidyr::separate(coef, sep = "X", into = c(".", "karyotype")) %>% 
  dplyr::select(lfc, name, karyotype, omic, pval, adj_pval) %>% 
  mutate(karyotype = str_replace(karyotype, "\\.", ":")) %>% 
  dplyr::group_by(karyotype)

df_dna = readRDS("results/DNA_lfc.rds") %>% dplyr::rename(DNA_lfc = lfc)

df = dplyr::bind_rows(res_rna, res_prot) %>% 
  dplyr::filter(!is.na(lfc))
df$omic = factor(df$omic, levels = c("RNA", "Protein"))

df %>% saveRDS(file.path(RES_PATH, "lfc_prot_and_rna_bind.rds"))

if (use_only_common_genes) {
  df = df %>% 
    dplyr::group_by(karyotype, name) %>% 
    dplyr::filter(n() == 2)
}

df = df %>% 
  dplyr::left_join(karyotypes_df_good)

df_dna = df_dna %>% 
  dplyr::left_join(karyotypes_df_good)

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

corum_genes = read.delim("data/CORUM_gene_list.txt") %>% dplyr::pull(GeneSym)
df = df %>% dplyr::mutate(group = ifelse(name %in% corum_genes, "complex", "non-complex"))
df %>% ungroup() %>% dplyr::select(name, group) %>% dplyr::distinct() %>% dplyr::pull(group) %>% table()

df = df %>% dplyr::filter(omic != "DNA")

# df = df %>% 
#   tidyr::separate(karyotype, into = c("A", "B"), sep = ":", remove = FALSE, convert = TRUE) %>% 
#   dplyr::mutate(copies = A + B, DNA_lfc = log2(copies / 2)) %>% 
#   dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc))

df = df %>% 
  ungroup() %>% 
  #dplyr::left_join(df_dna %>% dplyr::select(!omic), by = join_by("name" == "name", "karyotype" == "karyotype")) %>% 
  dplyr::left_join(df_dna %>% dplyr::select(DNA_lfc, name, karyotype)) %>% 
  dplyr::mutate(CS = ifelse(DNA_lfc > 0, DNA_lfc - lfc, lfc - DNA_lfc)) 

# Volcano plot

df_volcano_and_bar = df %>% 
  dplyr::select(pval, lfc, name, karyotype, omic) %>% 
  dplyr::group_by(omic, karyotype) %>% 
  dplyr::mutate(adj_pval = p.adjust(pval, "BH") + 1e-300) %>% 
  dplyr::mutate(signif = ifelse(adj_pval <= .05, "< 0.05", "ns")) %>% 
  dplyr::mutate(FC_class = ifelse((adj_pval <= .05) & (lfc > 0.75), "Up", ifelse((adj_pval <= .05) & (lfc < -0.75), "Down", "Not differential"))) %>% 
  dplyr::mutate(FC_class = factor(FC_class, levels = c("Not differential", "Up", "Down"))) %>% 
  dplyr::mutate(signif = factor(signif, levels = c("ns", "< 0.05"))) %>% 
  dplyr::mutate(karyotype = karyotype_mapping[karyotype]) %>% 
  dplyr::mutate(karyotype = factor(karyotype, levels = karyotype_mapping))
  
p_volcano = df_volcano_and_bar %>% 
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval), col = FC_class, alpha = signif)) +
  geom_point(size = .75) +
  facet_grid(omic~karyotype, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "darkblue", "Not differential" = "gray50")) +
  labs(x = "logFC", y = "-log10(adj-p-value)", col = "DE class", alpha = "Adjusted p-value")
p_volcano

ggsave(file.path(IMG_PATH, "volcano_full.pdf"), width = 10, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "volcano_full.png"), width = 10, height = 5, units = "in", dpi = 450)
saveRDS(p_volcano, file.path(IMG_PATH, "volcano_full.rds"))

p_bar = df_volcano_and_bar %>% 
  dplyr::filter(adj_pval <= .05) %>% 
  dplyr::group_by(FC_class, karyotype, omic) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(mapping = aes(x = "", y = n, fill = FC_class)) +
  geom_col(position = "dodge2") +
  facet_grid(omic~karyotype) +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "firebrick", "Down" = "darkblue", "Not differential" = "gray50")) +
  theme_bw() +
  labs(x = "", y = "Number of genes", fill = "DE class") +
  theme(axis.ticks.y = element_blank())
p_bar

ggsave(file.path(IMG_PATH, "deg_barplot_full.pdf"), width = 10, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "deg_barplot_full.png"), width = 10, height = 5, units = "in", dpi = 450)
saveRDS(p_bar, file.path(IMG_PATH, "deg_barplot_full.rds"))


p_volcano_sub = df_volcano_and_bar %>% 
  dplyr::filter(karyotype %in% c("LOH")) %>% 
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval), col = FC_class, alpha = signif)) +
  geom_point(size = .75) +
  facet_wrap(~omic, scales = "free_y", nrow = 2) +
  theme_bw() +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "darkblue", "Not differential" = "gray50")) +
  labs(x = "logFC", y = "-log10(adj-p-value)", col = "DE class", alpha = "Adjusted p-value")
p_volcano_sub

ggsave(file.path(IMG_PATH, "volcano_sub.pdf"), width = 6, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "volcano_sub.png"), width = 6, height = 5, units = "in", dpi = 450)
saveRDS(p_volcano_sub, file.path(IMG_PATH, "volcano_sub.rds"))

p_bar_sub = df_volcano_and_bar %>% 
  dplyr::mutate(FC_class = factor(FC_class, levels = c("Down", "Not differential", "Up"))) %>% 
  dplyr::filter(karyotype %in% c("LOH")) %>% 
  dplyr::filter(adj_pval <= .05) %>% 
  dplyr::group_by(FC_class, karyotype, omic) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(mapping = aes(x = "", y = n, fill = FC_class)) +
  geom_col(position = "dodge2") +
  facet_wrap(~omic, nrow = 2) +
  scale_fill_manual(values = c("Up" = "firebrick", "Down" = "darkblue", "Not differential" = "gray50")) +
  theme_bw() +
  labs(x = "", y = "Number of genes", fill = "DE class") +
  theme(axis.ticks.x = element_blank())
p_bar_sub

ggsave(file.path(IMG_PATH, "dge_barplot_sub.pdf"), width = 6, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "dge_barplot_sub.png"), width = 6, height = 5, units = "in", dpi = 450)
saveRDS(p_bar_sub, file.path(IMG_PATH, "dge_barplot_sub.rds"))
  
  
# # Perform CS computation per omic, karyotype, and group ####
# dir.create(file.path(RES_PATH, "CS_tables"), recursive = T)
# 
# N_BOOT = 5000
# CS_df = dplyr::tibble()
# for (k in unique(df$karyotype)) {
#   for (o in unique(df$omic)) {
#     for (g in unique(df$group)) {
#       df_sub = df %>% dplyr::filter(karyotype == k, omic == o, group == g) %>% dplyr::ungroup()
#       wilcox.test(df_sub$CS, mu = 0, alternative = "greater")
#       
#       CS_boot = lapply(1:N_BOOT, function(i) {
#         indexes = sample(1:nrow(df_sub), size = 50)
#         df_sub$CS[indexes] %>% median(na.rm = TRUE)
#       }) %>% unlist()
#       
#       p_val = mean(CS_boot <= 0)
#       
#       # wilcox.test(CS_boot, mu = 0, alternative = "greater")
#       # result <- t.test(CS_boot, mu = 0, alternative = "greater")
#       
#       CS_df = dplyr::bind_rows(
#         CS_df, 
#         dplyr::tibble(
#           group = g, karyotype = k, omic = o, median_CS = median(CS_boot), CI_low = quantile(CS_boot, probs = .05), CI_high = quantile(CS_boot, probs = .95), p_val = p_val
#         )  
#       )
#     }
#   }
# }
# df %>% dplyr::select(name, group) %>% 
#   dplyr::distinct() %>% 
#   dplyr::count(group)
# CS_df$adj_pval = p.adjust(CS_df$p_val, method = "BH")
# saveRDS(CS_df, file.path(RES_PATH, "CS_tables", "CS_table.RDS"))
# # sheet_write(CS_df, ss = ss, sheet = "Compensation scores")
# 
# # Perform CS between complex and non-complex ####
# N_BOOT = 5000
# CS_comp_v_noncomp_df = dplyr::tibble()
# for (k in unique(df$karyotype)) {
#   for (o in unique(df$omic)) {
#     df_sub = df %>% dplyr::filter(karyotype == k, omic == o) %>% dplyr::ungroup()
#     df_sub_complex = df_sub %>% dplyr::filter(group == "complex")
#     df_sub_non_complex = df_sub %>% dplyr::filter(group == "non-complex")
#     CS_boot = lapply(1:N_BOOT, function(i) {
#       complex_indexes = sample(1:nrow(df_sub_complex), size = 50)
#       non_complex_indexes = sample(1:nrow(df_sub_non_complex), size = 50)
#       complex_CS = df_sub_complex$CS[complex_indexes] %>% median(na.rm = T)
#       non_complex_CS = df_sub_non_complex$CS[non_complex_indexes] %>% median(na.rm = T)
#       complex_CS - non_complex_CS
#     }) %>% unlist()
#     
#     # Effect size + CI (bootstrap)
#     median_CS = median(CS_boot)
#     CI_low = quantile(CS_boot, probs = alpha)
#     CI_high = quantile(CS_boot, probs = 1 - alpha)
#     
#     # p-value (Wilcoxon)
#     p_val = wilcox.test(df_sub_complex$CS, df_sub_non_complex$CS, 
#                         alternative = "greater", na.action = na.omit)$p.value
#     
#     CS_comp_v_noncomp_df = dplyr::bind_rows(
#       CS_comp_v_noncomp_df, 
#       dplyr::tibble(
#         karyotype = k, omic = o, median_CS = median_CS, CI_low = CI_low, CI_high = CI_high, p_val = p_val
#       )  
#     )
#   }
# }
# CS_comp_v_noncomp_df$adj_pval = p.adjust(CS_comp_v_noncomp_df$p_val, method = "BH")
# saveRDS(CS_comp_v_noncomp_df, file.path(RES_PATH, "CS_tables", "CS_comp_vs_noncomp_table.RDS"))
# # sheet_write(CS_comp_v_noncomp_df, ss = ss, sheet = "Compensation scores (complex vs non-complex)")

# DNA depth ratio plot ####
df_expected_ratio = dplyr::tibble(
  karyotype = c("1:0", "2:0", "2:1", "2:2"),
  value = c(-1.0, 0, 0.5, 1)
)
df_expected_ratio$karyotype = karyotype_mapping[df_expected_ratio$karyotype]

p_dna = df_dna %>%
  dplyr::mutate(karyotype = karyotype_mapping[karyotype]) %>% 
  dplyr::mutate(karyotype = factor(karyotype, levels = karyotype_mapping)) %>% 
  #dplyr::mutate(group = ifelse(name %in% corum_genes, "complex", "non-complex")) %>% 
  ggplot(mapping = aes(x = karyotype, y = DNA_lfc)) +
  #geom_violin() +
  geom_boxplot(outliers = F, fill = omic_colors["DNA"]) +
  #geom_boxplot(outliers = F, width = .1, fill = omic_colors["DNA"]) +
  
  #scale_fill_manual(values = c("non-complex" = "#E69F00", "complex" = "#56B4E9")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(df_expected_ratio, mapping = aes(x = karyotype, y = value), inherit.aes = F, col = "firebrick3") +
  theme_bw() +
  labs(y = "DNA log-ratio") +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right"
  )
p_dna

ggsave(file.path(IMG_PATH, "DNA_ratio.pdf"), width = 5, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "DNA_ratio.png"), width = 5, height = 5, units = "in", dpi = 450)
saveRDS(p_dna, file.path(IMG_PATH, "DNA_ratio.rds"))

# RNA vs DNA plot of compensation fraction ####
compare_omics_boot <- function(sub_df, n_genes = 500, n_boot = 1000) {
  
  rna_labels <- sub_df$CS_class[sub_df$omic == "RNA"]
  prot_labels <- sub_df$CS_class[sub_df$omic == "Protein"]
  
  if(length(rna_labels) == 0 | length(prot_labels) == 0) return(NULL)
  
  boot_results <- replicate(n_boot, {
    
    s_rna  <- sample(rna_labels, size = n_genes, replace = (length(rna_labels) < n_genes))
    s_prot <- sample(prot_labels, size = n_genes, replace = (length(prot_labels) < n_genes))
    
    f_rna  <- sum(s_rna == "CS>0") / n_genes
    f_prot <- sum(s_prot == "CS>0") / n_genes
    
    c(f_rna = f_rna, f_prot = f_prot, diff = f_prot - f_rna)
    
  }) %>% t() %>% as.data.frame()
  
  p_val <- sum(boot_results$diff <= 0) / n_boot
  
  data.frame(
    omic = c("RNA", "Protein"),
    f_mean = c(mean(boot_results$f_rna), mean(boot_results$f_prot)),
    low = c(quantile(boot_results$f_rna, 0.025), quantile(boot_results$f_prot, 0.025)),
    high = c(quantile(boot_results$f_rna, 0.975), quantile(boot_results$f_prot, 0.975)),
    p_prot_gt_rna = p_val
  )
}

# 2. Run analysis
plot_data <- df %>%
  dplyr::mutate(CS_class = ifelse(CS > 0, "CS>0", "CS<0")) %>%
  dplyr::group_by(karyotype) %>%
  nest() %>%
  dplyr::mutate(stats = map(data, ~compare_omics_boot(.x, n_genes = 500))) %>%
  unnest(stats) %>%
  dplyr::select(-data) %>%
  mutate(sig_50 = ifelse(low > 0.5, "*", "")) %>%
  mutate(
    p_prot_gt_rna = ifelse(
      p_prot_gt_rna <= 0.001, "**",
      ifelse(
        p_prot_gt_rna <= 0.05, "*",
        "ns"
      )
    )
  ) %>%
  mutate(karyotype = karyotype_mapping[karyotype]) %>%
  mutate(karyotype = factor(karyotype, levels = karyotype_mapping))

# 3. Create bracket annotation dataset
pvals_df <- plot_data %>%
  group_by(karyotype) %>%
  summarise(
    y_position = max(high) + 0.075,
    p_label = unique(p_prot_gt_rna),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = as.numeric(karyotype) - 0.2,
    xmax = as.numeric(karyotype) + 0.2,
    annotation = paste0(p_label)
  )

# 4. Plot
p_compensated_genes_fraction <- ggplot(plot_data, aes(x = karyotype, y = f_mean, fill = omic)) +
  
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # Star if >50%
  geom_text(
    aes(label = sig_50, y = high + 0.02),
    position = position_dodge(width = 0.8),
    size = 6
  ) +
  
  # Brackets with p-values
  geom_signif(
    data = pvals_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      annotations = annotation,
      y_position = y_position, group = karyotype
    ),
    manual = TRUE,
    textsize = 3.5,
    tip_length = 0.01, inherit.aes = F
  ) +
  
  theme_bw() +
  scale_fill_manual(values = omic_colors) +
  labs(
    y = "Fraction of compensated genes",
    fill = "",
    x = "Karyotype"
  )

# p_compensated_genes_fraction +
#   coord_flip()

ggsave(file.path(IMG_PATH, "compensated_frac.pdf"), width = 5, height = 4, units = "in", plot = p_compensated_genes_fraction)
ggsave(file.path(IMG_PATH, "compensated_frac.png"), width = 5, height = 4, units = "in", dpi = 450, plot = p_compensated_genes_fraction)
saveRDS(p_compensated_genes_fraction, file.path(IMG_PATH, "compensated_frac.rds"))

# # 1. The Main Boxplot
# p_dna = df_dna %>%
#   dplyr::mutate(group = ifelse(name %in% corum_genes, "complex", "non-complex")) %>% 
#   ggplot(mapping = aes(x = karyotype, y = DNA_lfc)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("non-complex" = "#E69F00", "complex" = "#56B4E9")) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#   theme_bw() +
#   labs(y = "DNA lfc") +
#   theme(
#     strip.background = element_rect(fill = "white", color = "white"),
#     strip.text = element_text(face = "bold", size = 11),
#     legend.position = "right"
#   )
# 
# p_main <- df %>% 
#   ggplot(aes(x = karyotype, y = lfc, fill = group)) +
#   geom_boxplot(outliers = FALSE, alpha = 0.8, width = 0.7) +
#   facet_wrap(~omic) +
#   scale_fill_manual(values = c("non-complex" = "#E69F00", "complex" = "#56B4E9")) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white", color = "white"),
#     strip.text = element_text(face = "bold", size = 11),
#     legend.position = "right"
#   )
# 
# # 2. The Compensation Tile (Middle)
# p_comp <- CS_df %>%
#   mutate(signif_label = ifelse(adj_pval < 0.05, "*", "")) %>%
#   ggplot(aes(x = karyotype, y = group)) + # Simplified x to match p_main
#   geom_tile(aes(fill = median_CS), color = "white", size = 0.5) +
#   geom_text(aes(label = signif_label), vjust = 0.8) +
#   facet_wrap(~omic) +
#   scale_fill_gradientn(
#     colours = c("#f6f6f6", "#b1eeec", "#00aaa9"),
#     name = "Comp. Score"
#   ) +
#   theme_void() +
#   theme(
#     strip.text = element_blank(),
#     axis.text.y = element_text(size = 8), # Keep Y text for context
#     plot.margin = margin(b = 2)
#   )
# 
# # 3. The Significance Bar (Top)
# p_sig <- CS_comp_v_noncomp_df %>%
#   dplyr::mutate(signif_label = ifelse(adj_pval < alpha, "^", "")) %>%
#   dplyr::mutate(omic = factor(omic, levels = c("RNA", "Protein"))) %>% 
#   ggplot(aes(x = karyotype, y = 1)) +
#   geom_text(aes(label = signif_label), fontface = "bold") +
#   facet_wrap(~omic) +
#   theme_void() +
#   theme(
#     strip.text = element_blank(),
#     plot.margin = margin(b = 0)
#   )
# 
# # --- THE MAGIC ALIGNMENT ---
# # Wrap all plots and collect legends to prevent shifting
# final_plot <- (p_sig / p_comp / p_main) + 
#   plot_layout(heights = c(0.5, 1, 10), guides = "collect") & 
#   theme(
#     # Force the left margin to be the same across all plots
#     plot.margin = margin(l = 50, r = 5, t = 5, b = 5)
#   )
# final_plot
# 
# ggsave(plot = final_plot, filename = file.path(IMG_PATH, "RNA_protein_CS.pdf"), width = 10, height = 7)
# ggsave(plot = final_plot, filename = file.path(IMG_PATH, "RNA_protein_CS.png"), width = 10, height = 7, dpi = 450, units = "in")

# RNA v PROTEIN correlation
# df_long = df %>% 
#   dplyr::mutate(DNA = copies) %>% 
#   dplyr::select(lfc, name, karyotype, omic, group, DNA) %>% 
#   tidyr::pivot_wider(names_from = omic, values_from = lfc)

# df_long = df %>% 
#   dplyr::select(lfc, name, karyotype, omic, group) %>% 
#   tidyr::pivot_wider(names_from = omic, values_from = lfc)
# 
# annotations <- df_long %>%
#   group_by(karyotype, group) %>%
#   filter(!is.na(RNA) & !is.na(protein)) %>%
#   summarise(
#     rho = cor(RNA, protein, method = "spearman", use = "pairwise.complete.obs"),
#     slope = coef(lm(protein ~ RNA))[2], # Extracts the 'RNA' coefficient
#     .groups = "drop"
#   ) %>%
#   # Create a label string
#   mutate(label = paste0("rho == ", round(rho, 2), "\nslope == ", round(slope, 2))) %>% 
#   group_by(karyotype) %>%
#   mutate(y_pos = seq(max(df_long$protein, na.rm=T), by = -0.5, length.out = n()))
# 
# rna_prot_corr <- df_long %>% 
#   ggplot(mapping = aes(x = RNA, y = protein, color = group)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(method = "lm", se = TRUE) + 
#   geom_text(data = annotations, 
#             aes(label = label, x = -Inf, y = y_pos, color = group),
#             hjust = -0.1, vjust = 1.1,
#             parse = TRUE,              
#             show.legend = FALSE) + 
#   facet_wrap(~karyotype, nrow = 1) +
#   scale_color_manual(values = c("non-complex" = "#E69F00", "complex" = "#56B4E9")) +
#   theme_bw() +
#   labs(
#     x = "logFC RNA",
#     y = "logFC Protein"
#   )
# rna_prot_corr
# 
# ggsave(plot = rna_prot_corr, filename = file.path(IMG_PATH, "rna_prot_corr.pdf"), width = 14, height = 7)
# ggsave(plot = rna_prot_corr, filename = file.path(IMG_PATH, "rna_prot_corr.png"), width = 14, height = 7, dpi = 450, units = "in")
# correlation_table <- df_long %>%
#   group_by(karyotype, group) %>%
#   summarise(
#     spearman_rho = cor(RNA, protein, method = "spearman"),
#     p_val = cor.test(RNA, protein, method = "spearman")$p.value,
#     n_genes = n(),
#     .groups = "drop"
#   ) %>% 
#   dplyr::mutate(p_adj = p.adjust(p_val, "BH"))
# correlation_table
# 
# dir.create(file.path(RES_PATH, "RNA_protein_corr"), recursive = T)
# saveRDS(correlation_table, file.path(RES_PATH, "RNA_protein_corr", "CS_comp_vs_noncomp_table.RDS"))
# # sheet_write(correlation_table, ss = ss, sheet = "RNA-protein correlations")
# 
# correlation_spearman_barplot = correlation_table %>% 
#   ggplot(mapping = aes(x = karyotype, y = spearman_rho, fill = group)) +
#   geom_col(position = "dodge", width = .5) +
#   scale_fill_manual(values = c("non-complex" = "#E69F00", "complex" = "#56B4E9")) +
#   theme_bw(base_size = 15) +
#   labs(x = "Karyotype", y = "Spearman correlation")
# correlation_spearman_barplot
# 
# ggsave(plot = correlation_spearman_barplot, filename = file.path(IMG_PATH, "rna_prot_spearman_barplot.pdf"), width = 7, height = 5)
# ggsave(plot = correlation_spearman_barplot, filename = file.path(IMG_PATH, "rna_prot_spearman_barplot.png"), width = 7, height = 5, dpi = 450, units = "in")

df_mean_CS = df %>% 
  dplyr::group_by(name, omic) %>% 
  dplyr::summarise(mean_CS = mean(CS)) %>% 
  tidyr::pivot_wider(values_from = mean_CS, names_from = omic)

df_mean_CS %>% 
  ggplot(mapping = aes(x = RNA, y = Protein)) +
  geom_point()

q = .5

rna_cuts = df_mean_CS %>%
  dplyr::mutate(s = sign(RNA)) %>%
  dplyr::group_by(s) %>%
  dplyr::summarise(q = quantile(RNA, q)) %>%
  dplyr::pull(q) %>%
  sort()

prot_cuts = df_mean_CS %>%
  dplyr::mutate(s = sign(Protein)) %>%
  dplyr::group_by(s) %>%
  dplyr::summarise(q = quantile(Protein, q)) %>%
  dplyr::pull(q) %>%
  sort()

# 3. Assign the genes to the specific groups from your text
df_groups <- df_mean_CS %>%
  mutate(reg_group = case_when(
    # Group 1: High RNA CS (>65th) and Low Protein CS (<35th)
    RNA > rna_cuts[2] & Protein < prot_cuts[1] ~ "RNA-comp.",
    RNA > rna_cuts[2] & Protein > prot_cuts[2] ~ "Compensated",
    
    # Group 2: Low RNA CS (<35th) and High Protein CS (>65th)
    RNA < rna_cuts[1] & Protein > prot_cuts[2] ~ "Prot-comp.",
    RNA < rna_cuts[1] & Protein < prot_cuts[1] ~ "Hyper-responders",
    
    TRUE ~ "Intermediate/Other"
  ))

df_groups %>%
  dplyr::group_by(reg_group) %>%
  dplyr::count(reg_group)


df_beeswarm = dplyr::bind_rows(
  df %>% 
    dplyr::select(lfc, karyotype, name, omic) %>% 
    dplyr::left_join(df_groups %>% dplyr::select(name, reg_group)),  
  df_dna %>% 
    dplyr::mutate(omic = "DNA") %>% 
    dplyr::rename(lfc = DNA_lfc) %>% 
    dplyr::select(lfc, karyotype, name, omic) %>% 
    dplyr::left_join(df_groups %>% dplyr::select(name, reg_group))
) %>% 
  # na.omit() %>% 
  dplyr::mutate(karyotype = factor(karyotype, levels = c("1:0", "2:0", "2:1", "2:2"))) %>% 
  dplyr::mutate(karyotype = karyotype_mapping[karyotype]) %>% 
  dplyr::mutate(karyotype = factor(karyotype, levels = karyotype_mapping))

pd = position_dodge(width = 0.4)
p_omics_trend_by_reg_groups = df_beeswarm %>%
  dplyr::filter(reg_group != "Intermediate/Other") %>% 
  group_by(omic, reg_group, karyotype) %>%
  summarise(
    m = mean(lfc, na.rm = TRUE), 
    q75 = quantile(lfc, 0.75, na.rm = TRUE), 
    q25 = quantile(lfc, 0.25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(mapping = aes(
    x = as.factor(as.numeric(karyotype)), 
    y = m, 
    ymin = q25, 
    ymax = q75, 
    col = omic,
    group = omic  # Essential for geom_line to dodge correctly
  )) +
  # Use the position_dodge object defined above
  geom_line(position = pd, linewidth = 1) +
  geom_pointrange(position = pd) +
  # Clean up the look
  facet_wrap(~reg_group, nrow = 1) +
  scale_color_manual(values = omic_colors) + # High contrast for RNA vs Protein
  theme_bw() +
  labs(
    x = "Karyotype",
    y = "Log2 FC",
    color = "Omic"
  ) +
  scale_x_discrete(labels = c("LOH", "CNLOH", "Trisomy", "Tetrasomy")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_omics_trend_by_reg_groups

ggsave(file.path(IMG_PATH, "omics_trends.pdf"), width = 12, height = 3.5, units = "in", plot = p_omics_trend_by_reg_groups)
ggsave(file.path(IMG_PATH, "omics_trends.png"), width = 12, height = 3.5, units = "in", dpi = 450, plot = p_omics_trend_by_reg_groups)
saveRDS(p_omics_trend_by_reg_groups, file.path(IMG_PATH, "omics_trends.rds"))

# p_karyo_abundance_smooth = ggplot() +
#   #geom_boxplot(outliers = F) +
#   geom_boxplot(df_beeswarm, mapping = aes(x = as.factor(as.numeric(karyotype)), y = lfc, col = omic), outliers = F) +
#   #ggbeeswarm::geom_quasirandom(df_beeswarm %>% dplyr::group_by(omic, karyotype, reg_group) %>% dplyr::slice_sample(n = 50, replace = T), mapping = aes(x = as.factor(as.numeric(karyotype)), y = lfc, col = omic), bandwidth = .1, size = 0.5) +
#   geom_smooth(df_beeswarm, method = "lm", inherit.aes = F, mapping = aes(x = as.numeric(karyotype), y = lfc, col = omic), se = TRUE, alpha = .2) +
#   facet_wrap(~reg_group, nrow = 1) +
#   scale_fill_manual(values = omic_colors) +
#   scale_color_manual(values = omic_colors) +
#   theme_bw() +
#   lims(y = c(-2, 2))
# p_karyo_abundance_smooth


df_groups %>% 
  dplyr::left_join(df)

# df %>% 
#   dplyr::filter(name == "ENPP4")
# 
# df %>% 
#   dplyr::filter(name %in% c("APOBEC3C","CDA","ENPP4","GDA","LCMT2","MACROD1"))
# df_groups %>% 
#   dplyr::filter(name %in% c("APOBEC3C","CDA","ENPP4","GDA","LCMT2","MACROD1"))

df_groups = df_groups %>% dplyr::mutate(label = ifelse(name %in% genes_to_plot, name, ""))
reg_groups_plot = df_groups %>% 
  ggplot(mapping = aes(x = RNA, y = Protein, col = reg_group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  labs(x = "RNA CS", y = "Protein CS", col = "")


reg_groups_plot
ggsave(plot = reg_groups_plot, filename = file.path(IMG_PATH, "reg_groups_plot.pdf"), width = 5, height = 3.5)
ggsave(plot = reg_groups_plot, filename = file.path(IMG_PATH, "reg_groups_plot.png"), width = 5, height = 3.5, dpi = 450, units = "in")
saveRDS(reg_groups_plot, file.path(IMG_PATH, "reg_groups_plot.rds"))

reg_groups_plot_w_annotation = reg_groups_plot + ggrepel::geom_label_repel(mapping = aes(label = label), show.legend = F, max.overlaps = Inf)

ggsave(plot = reg_groups_plot_w_annotation, filename = file.path(IMG_PATH, "reg_groups_plot_annotated.pdf"), width = 5, height = 3.5)
ggsave(plot = reg_groups_plot_w_annotation, filename = file.path(IMG_PATH, "reg_groups_plot_annotated.png"), width = 5, height = 3.5, dpi = 450, units = "in")
saveRDS(reg_groups_plot_w_annotation, file.path(IMG_PATH, "reg_groups_plot_annotated.rds"))

# og and tsg by reg groups (excluding intermediate)

cols = c('Hugo_Symbol', 	
         'Entrez_Gene_ID', 
         'GRCh37_Isoform', 	
         'GRCh37_RefSeq',	
         'GRCh38_Isoform', 	
         'GRCh38_RefSeq', 	
         'Gene_Type', 
         'number_occurrences', 
         'OncoKB_Annotated',
         'MSK_IMPACT', 
         'MSK_HEME',
         'FOUNDATION_ONE',
         'FOUNDATION_ONE_HEME',
         'Vogelstein',
         'COSMIC_CGC',
         'Gene_Aliases')

cancer_genes = read.table('/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/utilities/cancerGeneList.tsv', sep = "\t", header = F, skip = 1, col.names = cols) %>%
  dplyr::as_tibble() %>% 
  filter(Gene_Type %in% c('ONCOGENE', 'TSG', 'ONCOGENE_AND_TSG')) %>% 
  dplyr::select(Hugo_Symbol, Gene_Type) %>% 
  distinct() %>% 
  mutate(Gene_Type = case_when(
    Gene_Type == 'ONCOGENE' ~ 'Oncogene', 
    Gene_Type == 'ONCOGENE_AND_TSG' ~ 'Oncogene/TSG', 
    .default = Gene_Type
  )) %>% 
  mutate(Gene_Type = factor(Gene_Type, levels = c('TSG', 'Oncogene', 'Oncogene/TSG')))

gene_roles_groups = df_groups %>% 
  dplyr::right_join(., cancer_genes, by = join_by('name' == 'Hugo_Symbol')) %>% 
  dplyr::rename(role = 'Gene_Type') %>% 
  filter(!is.na(reg_group)) 

tsg_og_reg_group_dist = gene_roles_groups %>% 
  group_by(role,reg_group) %>% 
  summarise(n_genes = n()) %>% 
  group_by(role) %>% 
  mutate(tot = sum(n_genes), 
         frac = n_genes/tot) %>% 
  filter(reg_group != 'Intermediate/Other') %>% 
  ggplot(aes(
    y = role, 
    x = frac,
    fill = reg_group
  )) + 
  geom_bar(stat = 'identity', 
           # position = 'dodge',
           position = position_dodge2(preserve = "single", width = 1),
           alpha = 1, 
           # width = 
  ) + 
  theme_bw()  +
  scale_fill_manual(values = category_colors) + 
  coord_flip() +
  labs(
    y = 'Gene role', 
    x = '% genes'
  ) + 
  guides(fill = guide_legend(title = '')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        # legend.position = 'bottom',
        # legend.key.size = unit(.2, 'cm'), 
        legend.key.height = unit(.2, 'cm'),
        legend.key.width = unit(.2, 'cm'),
        legend.text = element_text(size=8)) + 
  scale_x_continuous(labels = scales::percent) 

tsg_og_reg_group_dist

ggsave(plot = tsg_og_reg_group_dist, filename = file.path(IMG_PATH, "tsg_og_reg_group_dist.png"), height = 5, width = 4, dpi = 450, units = "in")
ggsave(plot = tsg_og_reg_group_dist, filename = file.path(IMG_PATH, "tsg_og_reg_group_dist.pdf"), height = 5, width = 4)
saveRDS(tsg_og_reg_group_dist, file.path(IMG_PATH, "tsg_og_reg_group_dist.rds"))


# wrap everything 

pt_l = '
AAAA
BBCC
'

p_omics_trend_by_reg_groups = readRDS(file.path(IMG_PATH, "omics_trends.rds"))

wrap_plots(list(
  p_omics_trend_by_reg_groups
))



# annotation
groups = setdiff(unique(df_groups$reg_group), "Intermediate/Other")
genes_by_cluster = lapply(groups, function(g) {
  df_groups %>% dplyr::filter(reg_group == g) %>% dplyr::pull(name)
})
names(genes_by_cluster) <- groups

# sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-prot heavy)`) / length(genes_by_cluster$`(RNA-prot heavy)`)
# sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-Prot light)`) / length(genes_by_cluster$`(RNA-Prot light)`)
# sum(genes_with_4_karyotypes %in% genes_by_cluster$`(RNA-prot heavy)`) / length(genes_with_4_karyotypes)

formula_res <- compareCluster(genes_by_cluster, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db, 
                              keyType = "SYMBOL",
                              ont = "BP")


dir.create(file.path(RES_PATH, "enrichment"), recursive = T)
saveRDS(formula_res, file.path(RES_PATH, "enrichment", "cluster_comparison.rds"))
# formula_res@compareClusterResult %>%
#   sheet_write(ss = ss, sheet = "Cluster Enrichment")
dotplot(formula_res)


library(GOSemSim)
library(igraph)

TOP_K = 10

go_df <- formula_res@compareClusterResult %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::arrange(p.adjust) %>% 
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

ggsave(file.path(IMG_PATH, "semantic_enrichment.pdf"), plot = semantic_plot_overlay, width = 7, height = 5, units = "in")
ggsave(file.path(IMG_PATH, "semantic_enrichment.png"), plot = semantic_plot_overlay, width = 7, height = 5, units = "in", dpi = 450)
saveRDS(semantic_plot_overlay, file.path(IMG_PATH, "semantic_enrichment.rds"))

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
  
