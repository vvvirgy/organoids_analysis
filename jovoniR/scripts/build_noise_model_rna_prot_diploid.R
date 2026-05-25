
rm(list = ls())
library(tidyverse)

# ── Paths ──────────────────────────────────────────────────────────────────
PROT_BOOT_PATH <- "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/noise_model/protein/fc_tb_clean.rds"
RNA_BOOT_PATH  <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/diploid_bootstrap.rds"
OUT_PATH       <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise.RDS"
IMG_PATH       <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/img/"

# ── Load bootstrap data ────────────────────────────────────────────────────
rna_df <- readRDS(RNA_BOOT_PATH) %>%
  dplyr::filter(min_mean_expr > 0.05, non_zero_percent > 1, n_samples > 5) %>%
  dplyr::select(gene, iteration, lfc) %>%
  dplyr::rename(lfc_rna = lfc)

prot_df <- readRDS(PROT_BOOT_PATH) %>%
  dplyr::select(PG.Genes, iteration, B_vs_A_diff) %>%
  dplyr::rename(gene = PG.Genes, lfc_prot = B_vs_A_diff) %>%
  na.omit()

# Matched pairs for the buffering null: CS_Buffering = RNA_lfc - Protein_lfc
# requires draws from the same gene × iteration.
df_matched <- dplyr::inner_join(rna_df, prot_df, by = c("gene", "iteration")) %>%
  dplyr::mutate(lfc_buff = lfc_rna - lfc_prot)

# ── Fit noise models ───────────────────────────────────────────────────────
# A new observation x is tested as: z = (x - mu) / sigma, p = 2*pnorm(-|z|).
diploid_noise <- dplyr::tibble(
  quantity = c("CS_RNA", "CS_Protein", "CS_Buffering"),
  mu       = c(mean(rna_df$lfc_rna), mean(prot_df$lfc_prot), mean(df_matched$lfc_buff)),
  sigma    = c(sd(rna_df$lfc_rna),   sd(prot_df$lfc_prot),   sd(df_matched$lfc_buff)),
  n        = c(nrow(rna_df),         nrow(prot_df),           nrow(df_matched))
)

print(diploid_noise)
saveRDS(diploid_noise, OUT_PATH)

# ── Diagnostic plot ────────────────────────────────────────────────────────
df_density <- dplyr::bind_rows(
  dplyr::tibble(lfc = rna_df$lfc_rna,        quantity = "CS_RNA"),
  dplyr::tibble(lfc = prot_df$lfc_prot,      quantity = "CS_Protein"),
  dplyr::tibble(lfc = df_matched$lfc_buff,   quantity = "CS_Buffering")
) %>%
  dplyr::mutate(quantity = factor(quantity, levels = c("CS_RNA", "CS_Protein", "CS_Buffering")))

sigma_lines <- diploid_noise %>%
  dplyr::mutate(quantity = factor(quantity, levels = c("CS_RNA", "CS_Protein", "CS_Buffering")))

p <- ggplot(df_density, aes(x = lfc, fill = quantity)) +
  geom_density(alpha = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_vline(data = sigma_lines, aes(xintercept =  sigma), linetype = "dotted", colour = "firebrick3") +
  geom_vline(data = sigma_lines, aes(xintercept = -sigma), linetype = "dotted", colour = "firebrick3") +
  facet_wrap(~quantity, nrow = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "Diploid LFC (null)", y = "Density", subtitle = "Red dotted lines = ±1σ") +
  theme(legend.position = "none")

ggsave(file.path(IMG_PATH, "noise_model_null_densities.pdf"), p, width = 9, height = 3.5)
ggsave(file.path(IMG_PATH, "noise_model_null_densities.png"), p, width = 9, height = 3.5, dpi = 450)
