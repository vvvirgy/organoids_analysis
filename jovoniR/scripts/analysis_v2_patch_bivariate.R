# ===========================================================================
# analysis_v3.R  —  Bivariate noise model variant of analysis_v2.R
#
# What changed vs v2:
#   • Loads the bivariate diploid null (μ, Σ, ρ) from diploid_noise_bivariate.RDS
#   • Adds p_joint (Mahalanobis on (RNA, Protein)) used as a global "is this
#     gene's joint behaviour unusual?" gate before the 3-bit class is trusted.
#   • Adds p_Buffer_cond: conditional test Prot | RNA, which is the proper,
#     more powerful buffering signal (uses σ_P√(1−ρ²)).
#   • Marginal p_RNA / p_Protein / p_Buffer kept for interpretability; their
#     σ's now come from Σ so they are internally consistent.
#   • class_label is gated by joint significance; class_label_raw kept for
#     comparison.
#   • A new diagnostic section visualises the (RNA, Protein) cloud per
#     karyotype against the diploid 1σ/2σ Mahalanobis ellipses.
#   • Buffering biology sections (CORUM, chr 17 mito) are run twice, once with
#     marginal comp_Buffer and once with comp_Buffer_cond, so any power gain
#     from the conditional test is auditable.
# ===========================================================================

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/utils.R")
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/constants.R")
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/getters.R")

IMG_PATH <- paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/img/sf_", sf_method, "_stable_", use_stable, "_biv")
RES_PATH <- paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/sf_", sf_method, "_stable_", use_stable, "_biv")


dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)

dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
MIN_N_PER_GROUP <- 4   # per-karyotype minimum sample size filter

# ===========================================================================
# 1. Load data and BIVARIATE noise model
# ===========================================================================

df_dna <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/DNA_lfc.rds") %>%
  dplyr::rename(DNA_lfc = lfc)

df_raw <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/CS_scores_prot_and_rna.rds")

# Bivariate null: list(mu, Sigma, Sigma_inv, rho, cond_PgivenR, n)
noise_biv <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise_bivariate.RDS")

mu_R   <- noise_biv$mu[1];  sd_R <- sqrt(noise_biv$Sigma[1, 1])
mu_P   <- noise_biv$mu[2];  sd_P <- sqrt(noise_biv$Sigma[2, 2])
rho    <- noise_biv$rho
Sig    <- noise_biv$Sigma
Sinv   <- noise_biv$Sigma_inv

# Marginal nulls derived from Σ
mu_buf  <- mu_R - mu_P
sd_buf  <- sqrt(Sig[1, 1] + Sig[2, 2] - 2 * Sig[1, 2])

# Conditional Prot | RNA
beta_PgivenR  <- noise_biv$cond_PgivenR$beta
intc_PgivenR  <- noise_biv$cond_PgivenR$intercept
sd_PgivenR    <- sqrt(noise_biv$cond_PgivenR$sigma2)

cat(sprintf("Bivariate null: μ=(%.3f, %.3f), σ=(%.3f, %.3f), ρ=%.3f\n",
            mu_R, mu_P, sd_R, sd_P, rho))
cat(sprintf("Conditional Prot|RNA SD: %.3f  (vs marginal σ_P %.3f → ratio %.2f)\n",
            sd_PgivenR, sd_P, sd_PgivenR / sd_P))

# ===========================================================================
# 2. Build the unified per-(gene, karyotype) table
# ===========================================================================
df_wide <- df_raw %>%
  dplyr::select(name, karyotype, omic, lfc, pval) %>%
  tidyr::pivot_wider(names_from = omic, values_from = c(lfc, pval)) %>%
  dplyr::rename(RNA_lfc = lfc_RNA, Protein_lfc = lfc_Protein,
                RNA_pval = pval_RNA, Protein_pval = pval_Protein) %>%
  dplyr::left_join(df_dna %>% dplyr::select(name, karyotype, DNA_lfc),
                   by = c("name", "karyotype")) %>%
  dplyr::filter(!is.na(RNA_lfc), !is.na(Protein_lfc), !is.na(DNA_lfc)) %>%
  dplyr::mutate(
    CS_RNA     = DNA_lfc - RNA_lfc,
    CS_Protein = DNA_lfc - Protein_lfc,
    CS_Buffer  = RNA_lfc - Protein_lfc
  )
saveRDS(df_wide, file.path(RES_PATH, "df_wide_three_CS.rds"))

# ===========================================================================
# 3. Per-gene classification — joint + conditional + marginal
# ===========================================================================

# Mahalanobis squared on (RNA_lfc, Protein_lfc)
xy   <- as.matrix(df_wide[, c("RNA_lfc", "Protein_lfc")])
ctr  <- sweep(xy, 2, noise_biv$mu, "-")
mah2 <- rowSums((ctr %*% Sinv) * ctr)            # χ²₂ under the null

df_classified <- df_wide %>%
  dplyr::mutate(
    mahal2     = mah2,
    p_joint    = stats::pchisq(mah2, df = 2, lower.tail = FALSE),
    # Marginal one-sided tests (compensation = positive CS)
    p_RNA      = pnorm(CS_RNA,     mean = mu_R,   sd = sd_R,   lower.tail = FALSE),
    p_Protein  = pnorm(CS_Protein, mean = mu_P,   sd = sd_P,   lower.tail = FALSE),
    p_Buffer   = pnorm(CS_Buffer,  mean = mu_buf, sd = sd_buf, lower.tail = FALSE),
    # Conditional buffering: protein below what RNA predicts?
    Prot_resid     = Protein_lfc - (intc_PgivenR + beta_PgivenR * RNA_lfc),
    p_Buffer_cond  = pnorm(Prot_resid / sd_PgivenR, lower.tail = TRUE)
  ) %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(
    padj_joint        = p.adjust(p_joint,       method = "BH"),
    padj_RNA          = p.adjust(p_RNA,         method = "BH"),
    padj_Protein      = p.adjust(p_Protein,     method = "BH"),
    padj_Buffer       = p.adjust(p_Buffer,      method = "BH"),
    padj_Buffer_cond  = p.adjust(p_Buffer_cond, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    comp_RNA         = (CS_RNA     > 0) & (padj_RNA     <= ALPHA),
    comp_Protein     = (CS_Protein > 0) & (padj_Protein <= ALPHA),
    comp_Buffer      = (CS_Buffer  > 0) & (padj_Buffer  <= ALPHA),
    comp_Buffer_cond = (Prot_resid < 0) & (padj_Buffer_cond <= ALPHA),
    joint_sig        = padj_joint <= ALPHA
  ) %>%
  dplyr::mutate(
    class_3bit = paste0(
      as.integer(comp_RNA), as.integer(comp_Protein), as.integer(comp_Buffer)
    ),
    class_label_raw = dplyr::case_when(
      class_3bit == "000" ~ "Full dosage",
      class_3bit == "110" ~ "RNA-compensated",
      class_3bit == "011" ~ "Buffered (post-transcriptional)",
      class_3bit == "111" ~ "Compensated at both levels",
      class_3bit == "010" ~ "Protein-only compensated",
      class_3bit == "001" ~ "Buffering only",
      class_3bit == "100" ~ "RNA-only (unstable)",
      class_3bit == "101" ~ "RNA-comp + amplified protein",
      TRUE                 ~ "Other"
    ),
    # Joint-gated label: collapse to "Full dosage" if the gene's joint
    # (RNA, Protein) didn't escape the diploid bivariate cloud.
    class_label = ifelse(joint_sig, class_label_raw, "Full dosage")
  )

# Report how many genes the joint gate moved to "Full dosage"
gate_audit <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    n_total            = dplyr::n(),
    n_raw_compensated  = sum(class_label_raw != "Full dosage"),
    n_gated_compensated = sum(class_label != "Full dosage"),
    pct_demoted = round(100 * (n_raw_compensated - n_gated_compensated) /
                          pmax(n_raw_compensated, 1), 1),
    .groups = "drop"
  )
saveRDS(gate_audit, file.path(RES_PATH, "joint_gate_audit.rds"))
print(gate_audit)

# ---------------------------------------------------------------------------
# 3b. Up / Down / Not-differential classification per omic
# ---------------------------------------------------------------------------

classify_direction <- function(lfc, pval_adj, mu, sigma, alpha) {
  z          <- (lfc - mu) / sigma
  p_noise    <- 2 * pnorm(-abs(z))
  noise_ok   <- p_noise <= alpha
  de_ok      <- pval_adj <= alpha
  dplyr::case_when(
    !noise_ok | !de_ok ~ "Not differential",
    lfc > 0            ~ "Up",
    lfc < 0            ~ "Down",
    TRUE               ~ "Not differential"
  )
}

df_classified <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(
    RNA_padj     = p.adjust(RNA_pval,     method = "BH"),
    Protein_padj = p.adjust(Protein_pval, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    dir_RNA     = classify_direction(RNA_lfc,     RNA_padj,     mu_R, sd_R, ALPHA),
    dir_Protein = classify_direction(Protein_lfc, Protein_padj, mu_P, sd_P, ALPHA)
  )

saveRDS(df_classified, file.path(RES_PATH, "df_classified_three_CS.rds"))

# ---------------------------------------------------------------------------
# 3c. Up / Down summary per karyotype × omic
# ---------------------------------------------------------------------------

dir_counts <- df_classified %>%
  dplyr::select(name, karyotype, dir_RNA, dir_Protein) %>%
  tidyr::pivot_longer(c(dir_RNA, dir_Protein),
                      names_to = "omic", values_to = "direction") %>%
  dplyr::mutate(omic = sub("dir_", "", omic)) %>%
  dplyr::group_by(karyotype, omic, direction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype, omic) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()

dir_asymmetry <- dir_counts %>%
  dplyr::filter(direction %in% c("Up", "Down")) %>%
  dplyr::select(karyotype, omic, direction, n) %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  dplyr::mutate(
    n_diff    = Up + Down,
    asymmetry = ifelse(n_diff == 0, NA_real_, (Up - Down) / n_diff),
    binom_p   = purrr::map2_dbl(Up, Down, function(u, d) {
      if ((u + d) == 0) return(NA_real_)
      stats::binom.test(u, u + d, p = 0.5)$p.value
    })
  ) %>%
  dplyr::group_by(omic) %>%
  dplyr::mutate(binom_padj = p.adjust(binom_p, method = "BH")) %>%
  dplyr::ungroup()

saveRDS(dir_counts,    file.path(RES_PATH, "dir_counts_per_karyotype_omic.rds"))
saveRDS(dir_asymmetry, file.path(RES_PATH, "dir_asymmetry_per_karyotype_omic.rds"))

# ---------------------------------------------------------------------------
# 3d. Plots — Up/Down fractions and asymmetry
# ---------------------------------------------------------------------------

dir_counts_plot <- dir_counts %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    direction     = factor(direction, levels = c("Down", "Not differential", "Up")),
    omic          = factor(omic, levels = c("RNA", "Protein"))
  )

p_updown_frac <- ggplot(dir_counts_plot,
                        aes(x = karyotype_lab, y = frac, fill = direction)) +
  geom_col(position = "stack") +
  facet_wrap(~omic) +
  scale_fill_manual(values = c("Up"               = "firebrick",
                               "Down"             = "darkblue",
                               "Not differential" = "grey80")) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes", fill = "Direction")

ggsave(file.path(IMG_PATH, "up_down_fraction_per_karyotype.pdf"), p_updown_frac, width = 8, height = 4)
saveRDS(p_updown_frac, file.path(IMG_PATH, "up_down_fraction_per_karyotype.rds"))

asym_plot <- dir_asymmetry %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    omic          = factor(omic, levels = c("RNA", "Protein")),
    sig_label     = dplyr::case_when(
      is.na(binom_padj)      ~ "",
      binom_padj <= 0.001    ~ "***",
      binom_padj <= 0.01     ~ "**",
      binom_padj <= 0.05     ~ "*",
      TRUE                   ~ "ns"
    )
  )

p_asymmetry <- ggplot(asym_plot, aes(x = karyotype_lab, y = asymmetry, fill = omic)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sig_label,
                y = asymmetry + sign(asymmetry) * 0.04),
            position = position_dodge(width = 0.8), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "Karyotype",
       y = "(Up - Down) / (Up + Down)",
       fill = "Omic",
       caption = paste0("Binomial test vs 50/50, BH-adjusted. ",
                        "Positive = Up-skewed, Negative = Down-skewed."))

ggsave(file.path(IMG_PATH, "up_down_asymmetry.pdf"), p_asymmetry, width = 7, height = 4)
saveRDS(p_asymmetry, file.path(IMG_PATH, "up_down_asymmetry.rds"))

# ===========================================================================
# 3e. NEW: bivariate diagnostic — (RNA, Protein) cloud per karyotype with
#         diploid 1σ / 2σ Mahalanobis ellipses.
# ===========================================================================

make_ellipse <- function(mu, Sigma, k, npts = 200) {
  theta  <- seq(0, 2 * pi, length.out = npts)
  circle <- rbind(cos(theta), sin(theta))
  L      <- chol(Sigma)
  pts    <- t(mu + k * t(L) %*% circle)
  dplyr::tibble(x = pts[, 1], y = pts[, 2], level = paste0(k, "σ"))
}
ell_df <- dplyr::bind_rows(
  make_ellipse(noise_biv$mu, Sig, 1),
  make_ellipse(noise_biv$mu, Sig, 2)
)

p_joint_per_karyotype <- df_classified %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)) %>%
  ggplot(aes(RNA_lfc, Protein_lfc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = joint_sig), alpha = 0.35, size = 0.5) +
  geom_path(data = ell_df, aes(x, y, group = level, linetype = level),
            colour = "firebrick3", inherit.aes = FALSE, linewidth = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "darkorange"),
                      name = "Joint p < α") +
  scale_linetype_manual(values = c(`1σ` = "solid", `2σ` = "dashed"),
                        name = "Diploid null") +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC",
       subtitle = sprintf("Bivariate diploid null ρ = %.2f", rho))

ggsave(file.path(IMG_PATH, "joint_cloud_per_karyotype.pdf"),
       p_joint_per_karyotype, width = 9, height = 7)
saveRDS(p_joint_per_karyotype, file.path(IMG_PATH, "joint_cloud_per_karyotype.rds"))

# ===========================================================================
# 4. Sanity-check plot: distribution of LFCs vs DNA expectation
# ===========================================================================

df_long_lfc <- df_classified %>%
  dplyr::select(name, karyotype, DNA_lfc, RNA_lfc, Protein_lfc) %>%
  tidyr::pivot_longer(c(RNA_lfc, Protein_lfc), names_to = "omic", values_to = "lfc") %>%
  dplyr::mutate(omic = sub("_lfc", "", omic)) %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

p_lfc_density <- ggplot(df_long_lfc, aes(x = lfc, fill = omic, colour = omic)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~karyotype_lab, scales = "free_y") +
  scale_fill_manual(values = omic_colors) +
  scale_colour_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "log2 fold-change vs diploid", y = "Density")

ggsave(file.path(IMG_PATH, "lfc_density_vs_dosage.pdf"), p_lfc_density, width = 9, height = 6)
saveRDS(p_lfc_density, file.path(IMG_PATH, "lfc_density_vs_dosage.rds"))

# ===========================================================================
# 5. Genome-wide fraction compensated per karyotype, with bootstrap CIs
#    Now reported for FOUR axes: RNA, Protein, Buffer (marginal), Buffer_cond.
# ===========================================================================

boot_fraction <- function(x, n_boot = 1000) {
  if (length(x) == 0) return(c(mean = NA, low = NA, high = NA))
  reps <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(reps),
    low  = unname(quantile(reps, 0.025)),
    high = unname(quantile(reps, 0.975)))
}

frac_by_karyotype <- df_classified %>%
  dplyr::select(karyotype, comp_RNA, comp_Protein, comp_Buffer, comp_Buffer_cond) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein, comp_Buffer, comp_Buffer_cond),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(level = sub("comp_", "", level)) %>%
  dplyr::group_by(karyotype, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_comp)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)) %>%
  dplyr::mutate(level = factor(level,
                               levels = c("RNA", "Protein", "Buffer", "Buffer_cond")))

p_frac_compensated <- ggplot(frac_by_karyotype,
                             aes(x = karyotype_lab, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of compensated genes",
       fill = "Compensation level",
       caption = "Buffer = marginal R−P test; Buffer_cond = conditional Prot|RNA test.")

ggsave(file.path(IMG_PATH, "fraction_compensated_three_levels.pdf"), p_frac_compensated, width = 8, height = 4)
saveRDS(p_frac_compensated, file.path(IMG_PATH, "fraction_compensated_three_levels.rds"))

# ===========================================================================
# 6. RNA-LFC vs Protein-LFC scatter (buffering visualization)
#    Now with the conditional buffering call highlighted alongside the
#    marginal one, and the diploid 2σ ellipse for reference.
# ===========================================================================

p_rna_vs_prot <- df_classified %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
                buffer_status = dplyr::case_when(
                  comp_Buffer & comp_Buffer_cond ~ "Both",
                  comp_Buffer                    ~ "Marginal only",
                  comp_Buffer_cond               ~ "Conditional only",
                  TRUE                           ~ "Neither"
                )) %>%
  ggplot(aes(x = RNA_lfc, y = Protein_lfc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_path(data = ell_df %>% dplyr::filter(level == "2σ"),
            aes(x, y, group = level), inherit.aes = FALSE,
            colour = "firebrick3", linewidth = 0.4) +
  geom_point(aes(colour = buffer_status), alpha = 0.4, size = 0.6) +
  facet_wrap(~karyotype_lab) +
  scale_colour_manual(values = c("Neither"          = "grey75",
                                 "Marginal only"    = "skyblue3",
                                 "Conditional only" = "goldenrod3",
                                 "Both"             = "darkorange3")) +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC", colour = "Buffering call")

ggsave(file.path(IMG_PATH, "rna_vs_protein_lfc.pdf"), p_rna_vs_prot, width = 9, height = 6)
saveRDS(p_rna_vs_prot, file.path(IMG_PATH, "rna_vs_protein_lfc.rds"))

# ===========================================================================
# 7. Per-chromosome compensation analysis
# ===========================================================================

mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_chrom <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "chromosome_name"),
  filters    = "hgnc_symbol",
  values     = unique(df_classified$name),
  mart       = mart
) %>%
  dplyr::as_tibble() %>%
  dplyr::rename(name = hgnc_symbol, chr = chromosome_name) %>%
  dplyr::filter(chr %in% c(as.character(1:22), "X", "Y")) %>%
  dplyr::distinct(name, .keep_all = TRUE)

n_total  <- length(unique(df_classified$name))
n_mapped <- nrow(gene_chrom)
message(sprintf("Chromosome mapping: %d / %d genes mapped (%.1f%%)",
                n_mapped, n_total, 100 * n_mapped / n_total))

if (nrow(gene_chrom) == 0) stop("No genes mapped to chromosomes; aborting section 7.")

chrom_levels <- c(as.character(1:22), "X", "Y")

# ---- 7.1 Per-chromosome coverage ------------------------------------------
chrom_coverage <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(!is.na(chr)) %>%
  dplyr::group_by(chr, karyotype) %>%
  dplyr::summarise(
    n_genes           = dplyr::n(),
    n_buffer          = sum(comp_Buffer),
    n_buffer_cond     = sum(comp_Buffer_cond),
    n_rna             = sum(comp_RNA),
    n_prot            = sum(comp_Protein),
    frac_Buffer       = mean(comp_Buffer),
    frac_Buffer_cond  = mean(comp_Buffer_cond),
    frac_RNA          = mean(comp_RNA),
    frac_Protein      = mean(comp_Protein),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_levels))

saveRDS(chrom_coverage, file.path(RES_PATH, "chrom_coverage.rds"))

p_chrom_coverage <- chrom_coverage %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)) %>%
  ggplot(aes(x = chr, y = karyotype_lab, fill = n_genes)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_genes), size = 2.5) +
  scale_fill_gradient(low = "white", high = "darkgreen", name = "n genes") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Coverage: number of (gene, karyotype) combinations per chromosome")

ggsave(file.path(IMG_PATH, "chrom_coverage_heatmap.pdf"), p_chrom_coverage, width = 11, height = 4)
saveRDS(p_chrom_coverage, file.path(IMG_PATH, "chrom_coverage_heatmap.rds"))

# ---- 7.2 Long-form fractions for heatmap ----------------------------------
frac_by_chrom <- chrom_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::select(karyotype, chr,
                frac_RNA, frac_Protein, frac_Buffer, frac_Buffer_cond) %>%
  tidyr::pivot_longer(c(frac_RNA, frac_Protein, frac_Buffer, frac_Buffer_cond),
                      names_to = "level", values_to = "frac") %>%
  dplyr::mutate(level = sub("frac_", "", level),
                level = factor(level,
                               levels = c("RNA", "Protein", "Buffer", "Buffer_cond")),
                karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

# ---- 7.3 Enrichment vs karyotype-wide baseline ----------------------------
genome_wide <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    genome_RNA         = mean(comp_RNA),
    genome_Protein     = mean(comp_Protein),
    genome_Buffer      = mean(comp_Buffer),
    genome_Buffer_cond = mean(comp_Buffer_cond),
    .groups = "drop"
  )

enrichment_by_chrom <- chrom_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::left_join(genome_wide, by = "karyotype") %>%
  dplyr::mutate(
    enr_RNA         = log2((frac_RNA         + 0.01) / (genome_RNA         + 0.01)),
    enr_Protein     = log2((frac_Protein     + 0.01) / (genome_Protein     + 0.01)),
    enr_Buffer      = log2((frac_Buffer      + 0.01) / (genome_Buffer      + 0.01)),
    enr_Buffer_cond = log2((frac_Buffer_cond + 0.01) / (genome_Buffer_cond + 0.01))
  ) %>%
  dplyr::select(karyotype, chr,
                enr_RNA, enr_Protein, enr_Buffer, enr_Buffer_cond) %>%
  tidyr::pivot_longer(c(enr_RNA, enr_Protein, enr_Buffer, enr_Buffer_cond),
                      names_to = "level", values_to = "log2_enr") %>%
  dplyr::mutate(level = sub("enr_", "", level),
                level = factor(level,
                               levels = c("RNA", "Protein", "Buffer", "Buffer_cond")),
                karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

saveRDS(enrichment_by_chrom, file.path(RES_PATH, "chrom_enrichment.rds"))

# ---- 7.4 Fisher tests per (karyotype, chromosome) -------------------------
fisher_one <- function(a, b, a_other, b_other) {
  ft <- suppressWarnings(
    fisher.test(matrix(c(a, b, a_other, b_other), nrow = 2))
  )
  tibble::tibble(p = ft$p.value, OR = unname(ft$estimate))
}

fisher_per_chrom <- function(level_col) {
  df_classified %>%
    dplyr::left_join(gene_chrom, by = "name") %>%
    dplyr::filter(!is.na(chr)) %>%
    dplyr::mutate(is_comp = .data[[level_col]]) %>%
    dplyr::group_by(karyotype, chr) %>%
    dplyr::filter(dplyr::n() >= 10) %>%
    dplyr::summarise(a = sum(is_comp),
                     b = sum(!is_comp),
                     .groups = "drop") %>%
    dplyr::group_by(karyotype) %>%
    dplyr::mutate(a_other = sum(a) - a,
                  b_other = sum(b) - b) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(fisher_one(a, b, a_other, b_other)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(karyotype) %>%
    dplyr::mutate(fisher_padj = p.adjust(p, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(level = sub("comp_", "", level_col))
}

chrom_fisher <- dplyr::bind_rows(
  fisher_per_chrom("comp_RNA"),
  fisher_per_chrom("comp_Protein"),
  fisher_per_chrom("comp_Buffer"),
  fisher_per_chrom("comp_Buffer_cond")
) %>%
  dplyr::mutate(
    level = factor(level, levels = c("RNA", "Protein", "Buffer", "Buffer_cond")),
    chr   = factor(chr,   levels = chrom_levels),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

saveRDS(chrom_fisher, file.path(RES_PATH, "chrom_fisher_per_level.rds"))

# ---- 7.5–7.7 Heatmaps -----------------------------------------------------
p_chrom_heatmap_raw <- ggplot(frac_by_chrom,
                              aes(x = chr, y = karyotype_lab, fill = frac)) +
  geom_tile(colour = "white") +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "white", high = "darkblue",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Frac.\ncompensated") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Raw fraction of compensated genes per chromosome")

ggsave(file.path(IMG_PATH, "chrom_compensation_heatmap.pdf"), p_chrom_heatmap_raw, width = 11, height = 7)
saveRDS(p_chrom_heatmap_raw, file.path(IMG_PATH, "chrom_compensation_heatmap.rds"))

p_chrom_heatmap_enr <- ggplot(enrichment_by_chrom,
                              aes(x = chr, y = karyotype_lab, fill = log2_enr)) +
  geom_tile(colour = "white") +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0,
                       name = "log2 enr.\nvs karyotype\nbaseline") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Per-chromosome enrichment for compensation (above/below karyotype baseline)")

ggsave(file.path(IMG_PATH, "chrom_enrichment_heatmap.pdf"), p_chrom_heatmap_enr, width = 11, height = 7)
saveRDS(p_chrom_heatmap_enr, file.path(IMG_PATH, "chrom_enrichment_heatmap.rds"))

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ ""
  )
}

p_chrom_heatmap_fisher <- chrom_fisher %>%
  dplyr::mutate(
    log2_OR = log2(pmin(pmax(OR, 0.1), 10)),
    star    = sig_label(fisher_padj)
  ) %>%
  ggplot(aes(x = chr, y = karyotype_lab, fill = log2_OR)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 3) +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "log2(OR)\nvs rest of\ngenome") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Fisher enrichment per chromosome (stars = BH-adjusted p)")

ggsave(file.path(IMG_PATH, "chrom_fisher_heatmap.pdf"), p_chrom_heatmap_fisher, width = 11, height = 7)
saveRDS(p_chrom_heatmap_fisher, file.path(IMG_PATH, "chrom_fisher_heatmap.rds"))

# ---- 7.8 Chromosome-17 deep dive ------------------------------------------
df_chr17 <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::mutate(is_chr17 = !is.na(chr) & chr == "17")

# ===========================================================================
# 8. Gain vs loss stratification
# ===========================================================================

df_classified <- df_classified %>%
  tidyr::separate(karyotype, into = c("A", "B"), sep = ":",
                  remove = FALSE, convert = TRUE) %>%
  dplyr::mutate(
    total_copies = A + B,
    cn_direction = dplyr::case_when(
      total_copies < 2  ~ "Loss",
      total_copies == 2 ~ "Neutral (LOH)",
      total_copies > 2  ~ "Gain"
    )
  ) %>%
  dplyr::select(-A, -B)

frac_by_direction <- df_classified %>%
  dplyr::select(cn_direction,
                comp_RNA, comp_Protein, comp_Buffer, comp_Buffer_cond) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein, comp_Buffer, comp_Buffer_cond),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(level = sub("comp_", "", level)) %>%
  dplyr::group_by(cn_direction, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_comp)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(level = factor(level,
                               levels = c("RNA", "Protein", "Buffer", "Buffer_cond")))

p_frac_direction <- ggplot(frac_by_direction,
                           aes(x = cn_direction, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  theme_bw() +
  labs(x = "Copy-number change", y = "Fraction compensated",
       fill = "Level")

ggsave(file.path(IMG_PATH, "fraction_compensated_gain_vs_loss.pdf"),
       p_frac_direction, width = 7, height = 4)

# ===========================================================================
# 9. 3-bit class distribution per karyotype (joint-gated)
#    Side-by-side comparison: raw classes vs gated classes.
# ===========================================================================

class_dist_raw <- df_classified %>%
  dplyr::group_by(karyotype, class_label_raw) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
                gate          = "Raw (marginal only)") %>%
  dplyr::rename(class_label = class_label_raw)

class_dist_gated <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
                gate          = "Joint-gated")

class_dist_combined <- dplyr::bind_rows(class_dist_raw, class_dist_gated) %>%
  dplyr::mutate(gate = factor(gate, levels = c("Raw (marginal only)", "Joint-gated")))

p_class_stack <- ggplot(class_dist_combined,
                        aes(x = karyotype_lab, y = frac, fill = class_label)) +
  geom_col() +
  facet_wrap(~gate, ncol = 1) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes",
       fill = "Compensation class",
       caption = "Joint-gated: genes failing the bivariate Mahalanobis test are collapsed to Full dosage.") +
  theme(legend.position = "right")

ggsave(file.path(IMG_PATH, "class_distribution_per_karyotype.pdf"),
       p_class_stack, width = 9, height = 7)

# ===========================================================================
# 10. CORUM complex enrichment among buffered genes
#     Run with BOTH marginal and conditional buffering calls.
# ===========================================================================

corum_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/CORUM_gene_list.txt"
if (file.exists(corum_path)) {
  corum_genes <- read.delim(corum_path) %>% dplyr::pull(GeneSym)
  
  df_classified <- df_classified %>%
    dplyr::mutate(in_complex = name %in% corum_genes)
  
  run_corum_fisher <- function(buffer_col) {
    df_classified %>%
      dplyr::group_by(karyotype) %>%
      dplyr::summarise(
        p  = fisher.test(table(in_complex, .data[[buffer_col]]))$p.value,
        OR = fisher.test(table(in_complex, .data[[buffer_col]]))$estimate,
        n_buffered = sum(.data[[buffer_col]]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(padj = p.adjust(p, method = "BH"),
                    buffer_def = buffer_col)
  }
  
  fisher_buffered <- dplyr::bind_rows(
    run_corum_fisher("comp_Buffer"),
    run_corum_fisher("comp_Buffer_cond")
  )
  
  complex_enrichment <- df_classified %>%
    dplyr::group_by(karyotype, in_complex) %>%
    dplyr::summarise(
      frac_buffered       = mean(comp_Buffer),
      frac_buffered_cond  = mean(comp_Buffer_cond),
      frac_rna_comp       = mean(comp_RNA),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  saveRDS(complex_enrichment, file.path(RES_PATH, "complex_enrichment.rds"))
  saveRDS(fisher_buffered,    file.path(RES_PATH, "complex_fisher_buffered.rds"))
  
  # Plot: marginal vs conditional buffering OR per karyotype
  fisher_plot <- fisher_buffered %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
      buffer_def    = dplyr::recode(buffer_def,
                                    "comp_Buffer"      = "Marginal (R−P)",
                                    "comp_Buffer_cond" = "Conditional (Prot|RNA)"),
      log2_OR       = log2(pmax(OR, 0.1)),
      star          = sig_label(padj)
    )
  
  p_complex_or <- ggplot(fisher_plot,
                         aes(x = karyotype_lab, y = log2_OR, fill = buffer_def)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = star, y = log2_OR + sign(log2_OR) * 0.1),
              position = position_dodge(width = 0.8), size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw() +
    labs(x = "Karyotype", y = "log2(OR) for CORUM membership among buffered genes",
         fill = "Buffering definition",
         caption = "Higher OR = complex subunits more enriched in buffered set.")
  
  ggsave(file.path(IMG_PATH, "complex_buffering_marginal_vs_conditional.pdf"),
         p_complex_or, width = 8, height = 4)
  
  # Original-style bar plot using the marginal call (preserved for continuity)
  complex_plot_df <- complex_enrichment %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
    )
  
  pvals_complex <- fisher_buffered %>%
    dplyr::filter(buffer_def == "comp_Buffer") %>%
    dplyr::left_join(
      complex_plot_df %>%
        dplyr::group_by(karyotype_lab, karyotype) %>%
        dplyr::summarise(y_top = max(frac_buffered), .groups = "drop"),
      by = "karyotype"
    ) %>%
    dplyr::mutate(
      y_position = y_top + 0.05,
      xmin       = as.numeric(karyotype_lab) - 0.2,
      xmax       = as.numeric(karyotype_lab) + 0.2,
      annotation = paste0(sig_label(padj),
                          "  OR=", formatC(OR, digits = 2, format = "f"))
    )
  
  p_complex <- ggplot(complex_plot_df,
                      aes(x = karyotype_lab, y = frac_buffered, fill = in_complex)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_signif(
      data = pvals_complex,
      aes(xmin = xmin, xmax = xmax,
          annotations = annotation,
          y_position  = y_position,
          group       = karyotype_lab),
      manual = TRUE, textsize = 3.2, tip_length = 0.01, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "steelblue"),
                      labels = c("Non-complex", "CORUM complex")) +
    theme_bw() +
    labs(x = "Karyotype", y = "Fraction buffered (post-transcriptional)",
         fill = "",
         caption = "Marginal buffer call. Fisher's exact test, BH-adjusted.")
  
  ggsave(file.path(IMG_PATH, "complex_buffering.pdf"), p_complex, width = 7, height = 4)
}

# ===========================================================================
# 11. GO enrichment per compensation class (joint-gated classes)
# ===========================================================================

genes_by_class <- df_classified %>%
  dplyr::filter(class_label %in% c("RNA-compensated",
                                   "Buffered (post-transcriptional)",
                                   "Compensated at both levels",
                                   "Full dosage")) %>%
  dplyr::distinct(name, class_label) %>%
  split(.$class_label) %>%
  lapply(function(x) unique(x$name))

if (length(genes_by_class) > 1) {
  enrich_res <- compareCluster(
    geneClusters = genes_by_class,
    fun          = "enrichGO",
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    pAdjustMethod = "BH"
  )
  
  dir.create(file.path(RES_PATH, "enrichment"),
             recursive = TRUE, showWarnings = FALSE)
  saveRDS(enrich_res, file.path(RES_PATH, "enrichment", "compensation_class_GO.rds"))
  
  p_enrich <- dotplot(enrich_res, showCategory = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(IMG_PATH, "GO_by_class.pdf"), p_enrich, width = 10, height = 8)
}

# ===========================================================================
# 12. Chr 17 deep dive: GO, MitoCarta, LFC illustration
#     Now run twice with the two buffering definitions for the gene-set fisher
#     tests; GO uses the conditional set (more powerful) but exports both.
# ===========================================================================

# ---- 12.0 Gene sets per buffering definition ------------------------------
chr17_buffered_marg <- df_chr17 %>%
  dplyr::filter(chr == "17", comp_Buffer) %>%
  dplyr::pull(name) %>% unique()

chr17_buffered_cond <- df_chr17 %>%
  dplyr::filter(chr == "17", comp_Buffer_cond) %>%
  dplyr::pull(name) %>% unique()

chr17_all <- df_chr17 %>%
  dplyr::filter(chr == "17") %>%
  dplyr::pull(name) %>% unique()

all_measured <- unique(df_classified$name)

message(sprintf("chr 17 buffered (marginal): %d | (conditional): %d | chr 17 all: %d",
                length(chr17_buffered_marg),
                length(chr17_buffered_cond),
                length(chr17_all)))

# Save gene sets
writeLines(chr17_buffered_marg, file.path(RES_PATH, "chr17_buffered_marginal.txt"))
writeLines(chr17_buffered_cond, file.path(RES_PATH, "chr17_buffered_conditional.txt"))

# ---- 12.1 GO enrichment across BP / MF / CC -------------------------------
run_go_chr17 <- function(buffered_set, label) {
  lapply(c("BP", "MF", "CC"), function(ont) {
    res <- enrichGO(
      gene          = buffered_set,
      universe      = chr17_all,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res@result %>%
      dplyr::filter(p.adjust <= 0.05) %>%
      dplyr::mutate(ontology = ont, buffer_def = label)
  }) %>% dplyr::bind_rows()
}

go_marg <- run_go_chr17(chr17_buffered_marg, "Marginal")
go_cond <- run_go_chr17(chr17_buffered_cond, "Conditional")
go_by_ontology <- dplyr::bind_rows(go_marg, go_cond) %>%
  dplyr::arrange(buffer_def, ontology, p.adjust)

write_csv(go_by_ontology, file.path(RES_PATH, "GO_chr17_buffered_all_ontologies.csv"))
saveRDS(go_by_ontology,   file.path(RES_PATH, "GO_chr17_buffered_all_ontologies.rds"))

if (nrow(go_by_ontology) > 0) {
  top_terms <- go_by_ontology %>%
    dplyr::group_by(buffer_def, ontology) %>%
    dplyr::slice_min(p.adjust, n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Description = factor(Description, levels = unique(Description)),
      log_padj    = -log10(p.adjust)
    )
  
  p_go_combined <- ggplot(top_terms,
                          aes(x = log_padj, y = reorder(Description, log_padj),
                              colour = ontology, size = Count)) +
    geom_point() +
    facet_wrap(~buffer_def, scales = "free_y") +
    theme_bw() +
    labs(x = "-log10(adj. p)", y = NULL,
         colour = "Ontology", size = "Gene count",
         title = "GO enrichment: chr 17 buffered genes (vs chr 17 background)")
  
  ggsave(file.path(IMG_PATH, "GO_chr17_buffered_combined.pdf"),
         p_go_combined, width = 11, height = 7)
  saveRDS(p_go_combined, file.path(IMG_PATH, "GO_chr17_buffered_combined.rds"))
}

# ---- 12.2 MitoCarta 3.0 enrichment ----------------------------------------
mitocarta_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/Human.MitoCarta3.0.xls"

if (file.exists(mitocarta_path)) {
  mitocarta <- read_excel(mitocarta_path, sheet = "A Human MitoCarta3.0") %>%
    dplyr::as_tibble() %>%
    dplyr::pull(Symbol) %>%
    unique()
  
  run_mito_fisher <- function(scope_filter, buffer_col, scope_label) {
    df_classified %>%
      dplyr::left_join(gene_chrom, by = "name") %>%
      { if (scope_filter == "chr17") dplyr::filter(., chr == "17") else . } %>%
      dplyr::mutate(in_mito = name %in% mitocarta) %>%
      dplyr::group_by(karyotype) %>%
      { if (scope_filter == "chr17") dplyr::filter(., dplyr::n() >= 10) else . } %>%
      dplyr::summarise(
        a = sum( .data[[buffer_col]] &  in_mito),
        b = sum( .data[[buffer_col]] & !in_mito),
        c = sum(!.data[[buffer_col]] &  in_mito),
        d = sum(!.data[[buffer_col]] & !in_mito),
        .groups = "drop"
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        OR = fisher.test(matrix(c(a, b, c, d), nrow = 2))$estimate,
        p  = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(padj = p.adjust(p, method = "BH"),
                    scope = scope_label,
                    buffer_def = buffer_col)
  }
  
  mito_tests <- dplyr::bind_rows(
    run_mito_fisher("chr17",  "comp_Buffer",      "chr 17 — marginal"),
    run_mito_fisher("chr17",  "comp_Buffer_cond", "chr 17 — conditional"),
    run_mito_fisher("genome", "comp_Buffer",      "Genome-wide — marginal"),
    run_mito_fisher("genome", "comp_Buffer_cond", "Genome-wide — conditional")
  )
  saveRDS(mito_tests, file.path(RES_PATH, "mitocarta_enrichment_all_scopes.rds"))
  
  mito_combined <- mito_tests %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
      log2_OR       = log2(pmax(OR, 0.1)),
      star          = sig_label(padj)
    )
  
  p_mito <- ggplot(mito_combined,
                   aes(x = karyotype_lab, y = log2_OR, fill = scope)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = star,
                  y = log2_OR + sign(log2_OR) * 0.15),
              position = position_dodge(width = 0.8), size = 3.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw() +
    labs(x = "Karyotype",
         y = "log2(OR) MitoCarta enrichment in buffered set",
         fill = "Scope × buffer definition")
  
  ggsave(file.path(IMG_PATH, "mitocarta_buffer_enrichment.pdf"),
         p_mito, width = 9, height = 5)
  saveRDS(p_mito, file.path(IMG_PATH, "mitocarta_buffer_enrichment.rds"))
}

# ---- 12.3 Per-chromosome mito GO heatmap ---------------------------------
# Use the conditional buffer set (more powerful, cleaner signal expected).
buffered_by_chrom <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(comp_Buffer_cond, !is.na(chr)) %>%
  dplyr::distinct(name, chr) %>%
  split(.$chr) %>%
  lapply(function(x) unique(x$name))

bg_by_chrom <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(!is.na(chr)) %>%
  dplyr::distinct(name, chr) %>%
  split(.$chr) %>%
  lapply(function(x) unique(x$name))

keep_chroms <- names(buffered_by_chrom)[sapply(buffered_by_chrom, length) >= 15]

go_per_chrom <- lapply(keep_chroms, function(ch) {
  res <- enrichGO(
    gene          = buffered_by_chrom[[ch]],
    universe      = bg_by_chrom[[ch]],
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  res@result %>%
    dplyr::filter(p.adjust <= 0.05) %>%
    dplyr::mutate(chr = ch)
}) %>% dplyr::bind_rows()

saveRDS(go_per_chrom, file.path(RES_PATH, "GO_CC_buffered_per_chrom.rds"))

if (nrow(go_per_chrom) > 0) {
  mito_terms <- go_per_chrom %>%
    dplyr::filter(stringr::str_detect(Description,
                                      regex("mitochond|respiratory|oxidative phosphor|mitoribosom",
                                            ignore_case = TRUE))) %>%
    dplyr::mutate(
      chr      = factor(chr, levels = chrom_levels),
      log_padj = -log10(p.adjust)
    )
  
  if (nrow(mito_terms) > 0) {
    p_mito_per_chrom <- ggplot(mito_terms,
                               aes(x = chr, y = Description, fill = log_padj)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = Count), size = 2.8) +
      scale_fill_gradient(low = "white", high = "firebrick",
                          name = "-log10(adj.p)") +
      theme_bw() +
      labs(x = "Chromosome", y = NULL,
           title = "Mitochondrial GO CC enrichment among buffered (conditional) genes, per chromosome",
           caption = "Cells labeled with gene count.")
    
    ggsave(file.path(IMG_PATH, "mito_GO_per_chrom.pdf"),
           p_mito_per_chrom, width = 10, height = 5)
    saveRDS(p_mito_per_chrom, file.path(IMG_PATH, "mito_GO_per_chrom.rds"))
  }
}

# ---- 12.4 LFC illustration of chr 17 mito genes ---------------------------
mito_chr17_buffered_genes <- unique(unlist(strsplit(
  go_by_ontology %>%
    dplyr::filter(buffer_def == "Conditional",
                  stringr::str_detect(Description,
                                      regex("mitochond", ignore_case = TRUE))) %>%
    dplyr::pull(geneID),
  "/"
)))

if (length(mito_chr17_buffered_genes) > 0) {
  
  lfc_illustration <- df_classified %>%
    dplyr::filter(name %in% mito_chr17_buffered_genes) %>%
    dplyr::filter(karyotype %in% c("2:1", "2:2")) %>%
    dplyr::select(name, karyotype, DNA_lfc, RNA_lfc, Protein_lfc) %>%
    tidyr::pivot_longer(c(DNA_lfc, RNA_lfc, Protein_lfc),
                        names_to = "level", values_to = "lfc") %>%
    dplyr::mutate(
      level = factor(sub("_lfc", "", level), levels = c("DNA", "RNA", "Protein")),
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
    )
  
  p_lfc_mito <- ggplot(lfc_illustration,
                       aes(x = level, y = lfc, group = name)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(alpha = 0.25, colour = "grey40") +
    geom_point(alpha = 0.4, size = 1) +
    stat_summary(aes(group = 1), fun = median,
                 geom = "line", colour = "firebrick", linewidth = 1.2) +
    stat_summary(aes(group = 1), fun = median,
                 geom = "point", colour = "firebrick", size = 3) +
    facet_wrap(~karyotype_lab) +
    theme_bw() +
    labs(x = NULL, y = "log2 fold-change vs diploid",
         title = "Buffering of chr 17 mitochondrial genes (conditional set)",
         caption = "Each line = one gene. Red = median across genes.")
  
  ggsave(file.path(IMG_PATH, "chr17_mito_LFC_buffering.pdf"),
         p_lfc_mito, width = 8, height = 4)
  saveRDS(p_lfc_mito, file.path(IMG_PATH, "chr17_mito_LFC_buffering.rds"))
}

# ===========================================================================
# 13. Summary table
# ===========================================================================

summary_tbl <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class_label, values_from = n, values_fill = 0)

write_csv(summary_tbl, file.path(RES_PATH, "compensation_class_summary.csv"))

# Also a per-karyotype summary of how many genes pass each gate
gate_summary <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    n_total              = dplyr::n(),
    n_joint_sig          = sum(joint_sig),
    n_comp_RNA           = sum(comp_RNA),
    n_comp_Protein       = sum(comp_Protein),
    n_comp_Buffer        = sum(comp_Buffer),
    n_comp_Buffer_cond   = sum(comp_Buffer_cond),
    n_gated_compensated  = sum(class_label != "Full dosage"),
    .groups = "drop"
  )

write_csv(gate_summary, file.path(RES_PATH, "gate_summary_per_karyotype.csv"))
saveRDS(gate_summary,    file.path(RES_PATH, "gate_summary_per_karyotype.rds"))

cat("Done.\n")