# ===========================================================================
# analysis_v5.R  —  Angular classification in CS space
#
# What changed vs v4:
#
#   CLASSIFICATION — angular approach replacing marginal z-scores
#   ---------------------------------------------------------------
#   The core insight: once we have CS_RNA and CS_Protein (both sign-corrected
#   so positive = compensation), classifying a gene means asking TWO things:
#
#     1. Is the gene outside the diploid noise cloud at all?
#        → Mahalanobis distance in CS space, chi-squared gate at ALPHA.
#
#     2. If yes, in which DIRECTION did it escape?
#        → Angle θ of the (CS_RNA, CS_Protein) vector in WHITENED space
#          (each axis divided by its diploid σ before computing the angle,
#           so the noise ellipse becomes a circle and θ is meaningful).
#
#   This replaces the two independent marginal z-score tests with a single
#   geometric test. Key properties:
#     • The diploid RNA-Protein correlation (ρ = 0.32) is naturally accounted
#       for by computing Mahalanobis distance rather than Euclidean distance.
#     • Buffering falls out automatically from the angle: θ > 45° means
#       protein is more attenuated than RNA in the whitened space, i.e.
#       post-transcriptional buffering. No separate CS_Buffer test needed.
#     • The classification is in CS space (not LFC space), so the diploid
#       ellipse is correctly centered at (0,0) for all karyotypes.
#
#   ANGLE CONVENTION (in whitened CS space)
#   ----------------------------------------
#     θ = atan2(CS_Protein / sd_P,  CS_RNA / sd_R)   in degrees
#
#     -45° ≤ θ <  45°  →  "RNA compensated"          (RNA >> Protein)
#      45° ≤ θ < 135°  →  "Transcriptional + buffered" (both, Protein ≥ RNA)
#     135° ≤ θ < 180°  →  "Protein only"              (Protein, RNA near 0
#      or θ < -135°          or slightly negative)
#     θ in lower half  →  "Over-dosage / reversed"    (negative CS scores;
#       (< -45°)              gene amplified past DNA)
#
#   MAHALANOBIS GATE
#   ----------------
#   Mahalanobis distance is computed on (CS_RNA, CS_Protein) using the
#   DIPLOID covariance Σ_CS derived from Σ_LFC:
#
#     Σ_CS = A * Σ_LFC * A^T,   A = diag(sign(DNA_lfc))
#
#   Because sign(DNA_lfc) is ±1, Σ_CS = Σ_LFC (the sign flip does not
#   change the covariance structure). So we reuse Sinv directly.
#
#   A gene is classified as "Full dosage" if its Mahalanobis distance does
#   not exceed the chi-squared threshold at ALPHA (df = 2). Only genes
#   outside the cloud are assigned an angular class.
#
#   COVARIANCE SANITY CHECK (section 3a)
#   --------------------------------------
#   We check whether the CS correlation structure differs across karyotypes
#   from the diploid ρ. If it does, a karyotype-specific Mahalanobis would
#   be more appropriate.
#
#   RETAINED FROM v4
#   -----------------
#   • CS score computation (section 2)
#   • Up/Down direction classification using noise-model z-score (section 3b)
#   • All downstream analyses (chr analysis, CORUM, GO, MitoCarta, chr 17)
#     updated to use the new class_label vocabulary.
#   • Bootstrap CIs on compensation fractions (section 5)
#   • Bivariate diagnostic scatter now in CS space (section 3e) — this is
#     the plot that was broken in v4 because it used LFC space.
#
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

cat(paste0("Running with sf_method: ", sf_method,
           " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))

source("utils.R")
source("scripts/constants.R")
source("scripts/getters.R")

IMG_PATH        <- paste0("img/sf_",     sf_method, "_stable_", use_stable, "_v5")
RES_PATH        <- paste0("results/sf_", sf_method, "_stable_", use_stable, "_v5")
DF_DNA_PATH     <- "results/DNA_lfc.rds"
DF_CS_PATH      <- "results/multiOmic/CS_scores_prot_and_rna.rds"
NOISE_MODEL_PATH <- "results/multiOmic/diploid_noise_bivariate.RDS"

dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05

# Chi-squared threshold for the Mahalanobis gate (df = 2)
MAH_THR <- qchisq(1 - ALPHA, df = 2)   # ~5.99
cat(sprintf("Mahalanobis gate: chi2(df=2) threshold at alpha=%.2f → %.4f\n",
            ALPHA, MAH_THR))

# ===========================================================================
# 1. Load data and bivariate noise model
# ===========================================================================

df_dna <- readRDS(DF_DNA_PATH) %>%
  dplyr::rename(DNA_lfc = lfc)

df_raw <- readRDS(DF_CS_PATH)

# Bivariate null: list(mu, Sigma, Sigma_inv, rho, cond_PgivenR, n)
# Built on diploid (RNA_lfc, Protein_lfc); mu ≈ (0, 0)
noise_biv <- readRDS(NOISE_MODEL_PATH)

mu_R  <- noise_biv$mu[1];  sd_R <- sqrt(noise_biv$Sigma[1, 1])
mu_P  <- noise_biv$mu[2];  sd_P <- sqrt(noise_biv$Sigma[2, 2])
rho   <- noise_biv$rho
Sig   <- noise_biv$Sigma       # 2x2 covariance of (RNA_lfc, Protein_lfc)
Sinv  <- noise_biv$Sigma_inv   # its inverse

cat(sprintf("Bivariate null: mu=(%.4f, %.4f), sigma=(%.3f, %.3f), rho=%.3f\n",
            mu_R, mu_P, sd_R, sd_P, rho))

# The sign-flip A = diag(sign(DNA_lfc)) transforms (RNA_lfc, Protein_lfc) →
# (CS_RNA, CS_Protein). Because A is ±1 diagonal, Cov(CS) = A * Sig * A = Sig.
# So Sinv is also the inverse covariance of (CS_RNA, CS_Protein). No
# recomputation needed.

# ===========================================================================
# 2. Build per-(gene, karyotype) table with CS scores
# ===========================================================================

df_wide <- df_raw %>%
  dplyr::select(name, karyotype, omic, lfc, pval) %>%
  tidyr::pivot_wider(names_from = omic, values_from = c(lfc, pval)) %>%
  dplyr::rename(RNA_lfc      = lfc_RNA,      Protein_lfc  = lfc_Protein,
                RNA_pval     = pval_RNA,     Protein_pval = pval_Protein) %>%
  dplyr::left_join(df_dna %>% dplyr::select(name, karyotype, DNA_lfc),
                   by = c("name", "karyotype")) %>%
  dplyr::filter(!is.na(RNA_lfc), !is.na(Protein_lfc), !is.na(DNA_lfc)) %>%
  dplyr::mutate(
    # sign(DNA_lfc) correction: positive CS = compensation for both gains & losses
    CS_RNA     = (DNA_lfc - RNA_lfc)     * sign(DNA_lfc),
    CS_Protein = (DNA_lfc - Protein_lfc) * sign(DNA_lfc),
    CS_Buffer  = (RNA_lfc - Protein_lfc) * sign(DNA_lfc)   # kept for diagnostics
  )

saveRDS(df_wide, file.path(RES_PATH, "df_wide_CS.rds"))

# ===========================================================================
# 3. Angular classification in CS space
#
#   Step 1 — Mahalanobis gate: is (CS_RNA, CS_Protein) outside the diploid
#            cloud? Uses chi-squared(df=2) at ALPHA.
#
#   Step 2 — Angle in whitened CS space: determines compensation mode.
#            Whitening: divide each CS by its diploid sigma so the ellipse
#            becomes a circle and the angle is comparable across axes.
#
#            theta = atan2(CS_Protein/sd_P, CS_RNA/sd_R)  [degrees]
#
#   Angle → class mapping:
#
#     [ -45,  45)  "RNA compensated"
#                  RNA attenuation dominates; protein may also be attenuated
#                  but less so. Classic transcriptional dosage compensation.
#
#     [  45, 135)  "Transcriptional + buffered"
#                  Both RNA and protein are attenuated, but protein is MORE
#                  attenuated than RNA in the whitened space. This is the
#                  post-transcriptional buffering signal.
#
#     [ 135, 180)  "Protein only compensated"
#     [-180,-135)  CS_RNA ≈ 0 or slightly negative; protein compensated alone.
#                  Could reflect post-transcriptional-only regulation.
#
#     [  -45,-135) LOWER half of circle — CS scores are negative, meaning
#     (everything   the omic EXCEEDS the DNA expectation (over-dosage).
#      with theta   Biologically: gene is amplified past dosage prediction.
#      < -45°)
#
# ===========================================================================

# --- Mahalanobis on (CS_RNA, CS_Protein) using diploid Sinv ----------------
cs_mat <- as.matrix(df_wide[, c("CS_RNA", "CS_Protein")])
# Centre at (0,0) — mu ≈ 0 but subtract anyway for correctness
cs_ctr <- sweep(cs_mat, 2, c(mu_R, mu_P), "-")
mah2   <- rowSums((cs_ctr %*% Sinv) * cs_ctr)   # D^2 ~ chi2(df=2) under null

df_classified <- df_wide %>%
  dplyr::mutate(

    # Mahalanobis distance squared and gate
    mahal2       = mah2,
    outside_cloud = mahal2 > MAH_THR,

    # Whitened CS scores for angle computation
    cs_rna_w  = CS_RNA     / sd_R,
    cs_prot_w = CS_Protein / sd_P,

    # Angle in degrees: atan2(y, x) with y = whitened Protein, x = whitened RNA
    theta = atan2(cs_prot_w, cs_rna_w) * 180 / pi,

    # Angular compensation class — only assigned if outside the diploid cloud
    class_label = dplyr::case_when(

      # Inside cloud: indistinguishable from diploid noise
      !outside_cloud                        ~ "Full dosage",

      # Upper-right quadrant variants (positive CS for both — true compensation)
      theta >= -45  & theta <  45           ~ "RNA compensated",
      theta >=  45  & theta < 135           ~ "Transcriptional + buffered",

      # Left of vertical: protein compensated, RNA near zero or reversed
      theta >= 135  | theta < -135          ~ "Protein only compensated",

      # Lower half: negative CS, gene exceeds DNA dosage
      theta >= -135 & theta <  -45          ~ "Over-dosage"
    )

  )

# Convenience boolean flags derived from the angular class
# (useful for downstream Fisher tests and enrichment analyses)
df_classified <- df_classified %>%
  dplyr::mutate(
    comp_RNA     = class_label %in% c("RNA compensated",
                                      "Transcriptional + buffered"),
    comp_Protein = class_label %in% c("RNA compensated",
                                      "Transcriptional + buffered",
                                      "Protein only compensated"),
    comp_Buffer  = class_label == "Transcriptional + buffered"
  )

# Sanity checks
stopifnot(all(!is.na(df_classified$class_label)))
stopifnot(all(df_classified$outside_cloud | df_classified$class_label == "Full dosage"))

cat("Classification counts:\n")
print(dplyr::count(df_classified, class_label))

saveRDS(df_classified, file.path(RES_PATH, "df_classified_v5.rds"))

# ---------------------------------------------------------------------------
# 3a. Covariance sanity check: does ρ(CS_RNA, CS_Protein) vary by karyotype?
#     If yes, a karyotype-specific Mahalanobis would be more appropriate.
# ---------------------------------------------------------------------------

rho_by_karyotype <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    rho_cs    = cor(CS_RNA, CS_Protein, method = "spearman"),
    rho_lfc   = cor(RNA_lfc, Protein_lfc, method = "spearman"),
    n         = dplyr::n(),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

cat(sprintf("\nDiploid null rho = %.3f. Per-karyotype CS rho:\n", rho))
print(rho_by_karyotype)
saveRDS(rho_by_karyotype, file.path(RES_PATH, "rho_by_karyotype.rds"))

p_rho_check <- ggplot(rho_by_karyotype,
                      aes(x = karyotype_lab, y = rho_cs)) +
  geom_col(fill = "steelblue3", width = 0.6) +
  geom_hline(yintercept = rho, linetype = "dashed", colour = "firebrick3") +
  annotate("text", x = Inf, y = rho + 0.02,
           label = sprintf("Diploid rho = %.3f", rho),
           hjust = 1.1, colour = "firebrick3", size = 3.5) +
  theme_bw() +
  labs(x = "Karyotype",
       y = "Spearman rho(CS_RNA, CS_Protein)",
       title = "CS correlation per karyotype vs diploid null",
       caption = "Large deviations from dashed line suggest karyotype-specific covariance.")

ggsave(file.path(IMG_PATH, "rho_by_karyotype.pdf"), p_rho_check, width = 6, height = 4)
saveRDS(p_rho_check, file.path(IMG_PATH, "rho_by_karyotype.rds"))

# ---------------------------------------------------------------------------
# 3b. Up / Down / Not-differential classification per omic
#     Uses only the noise-model z-score (two-sided), no DE p-value.
# ---------------------------------------------------------------------------

classify_direction <- function(lfc, mu, sigma, alpha) {
  z       <- (lfc - mu) / sigma
  p_noise <- 2 * pnorm(-abs(z))
  dplyr::case_when(
    p_noise > alpha ~ "Not differential",
    lfc > mu        ~ "Up",
    lfc < mu        ~ "Down",
    TRUE            ~ "Not differential"
  )
}

df_classified <- df_classified %>%
  dplyr::mutate(
    dir_RNA     = classify_direction(RNA_lfc,     mu_R, sd_R, ALPHA),
    dir_Protein = classify_direction(Protein_lfc, mu_P, sd_P, ALPHA)
  )

# ---------------------------------------------------------------------------
# 3c. Up / Down summary per karyotype x omic
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
# 3d. Up/Down plots
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

ggsave(file.path(IMG_PATH, "up_down_fraction_per_karyotype.pdf"),
       p_updown_frac, width = 8, height = 4)
saveRDS(p_updown_frac, file.path(IMG_PATH, "up_down_fraction_per_karyotype.rds"))

p_asymmetry <- dir_asymmetry %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    omic          = factor(omic, levels = c("RNA", "Protein")),
    sig_label     = dplyr::case_when(
      is.na(binom_padj)   ~ "",
      binom_padj <= 0.001 ~ "***",
      binom_padj <= 0.01  ~ "**",
      binom_padj <= 0.05  ~ "*",
      TRUE                ~ "ns"
    )
  ) %>%
  ggplot(aes(x = karyotype_lab, y = asymmetry, fill = omic)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sig_label,
                y = asymmetry + sign(asymmetry) * 0.04),
            position = position_dodge(width = 0.8), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "Karyotype", y = "(Up - Down) / (Up + Down)", fill = "Omic",
       caption = "Binomial test vs 50/50, BH-adjusted.")

ggsave(file.path(IMG_PATH, "up_down_asymmetry.pdf"), p_asymmetry, width = 7, height = 4)
saveRDS(p_asymmetry, file.path(IMG_PATH, "up_down_asymmetry.rds"))

# ---------------------------------------------------------------------------
# 3e. Primary diagnostic scatter: CS space with Mahalanobis ellipse
#
#     This is the corrected version of the broken v4 plot. The ellipse is
#     now in CS space where it is correctly centered at (0,0) for all
#     karyotypes, and the angular class boundaries are overlaid.
# ---------------------------------------------------------------------------

# Build diploid ellipse in CS space
make_ellipse_cs <- function(mu_cs, Sigma_cs, k, npts = 300) {
  theta  <- seq(0, 2 * pi, length.out = npts)
  circle <- rbind(cos(theta), sin(theta))
  L      <- chol(Sigma_cs)
  pts    <- t(mu_cs + k * t(L) %*% circle)
  dplyr::tibble(x = pts[, 1], y = pts[, 2], level = paste0(k, "sigma"))
}

# Sigma_CS = Sigma_LFC (sign flip doesn't change covariance)
mu_cs  <- c(mu_R, mu_P)
ell_cs <- dplyr::bind_rows(
  make_ellipse_cs(mu_cs, Sig, sqrt(MAH_THR)),  # 95% boundary = sqrt(chi2 threshold)
  make_ellipse_cs(mu_cs, Sig, 1),
  make_ellipse_cs(mu_cs, Sig, 2)
)

# Angular boundary lines in CS space
# theta = 45 in whitened space → slope = sd_P / sd_R in CS space
slope_45  <-  sd_P / sd_R    # boundary between RNA-only and Buffered
slope_m45 <- -sd_P / sd_R    # boundary between RNA-only and Over-dosage
# theta = 135 → slope same as -45 but in second quadrant → vertical boundary
# at CS_RNA = 0 approximately; draw as vertical line at x = 0

class_colors <- c(
  "Full dosage"                  = "grey70",
  "RNA compensated"              = "steelblue3",
  "Transcriptional + buffered"   = "firebrick3",
  "Protein only compensated"     = "darkorange2",
  "Over-dosage"                  = "purple3"
)

p_cs_scatter <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    class_label   = factor(class_label, levels = names(class_colors))
  ) %>%
  ggplot(aes(x = CS_RNA, y = CS_Protein)) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  # Angular boundary: theta = 45 in whitened space
  geom_abline(slope = slope_45,  intercept = 0,
              linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_abline(slope = slope_m45, intercept = 0,
              linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  # Points coloured by angular class
  geom_point(aes(colour = class_label), alpha = 0.4, size = 0.5) +
  # Diploid ellipses
  geom_path(data = ell_cs,
            aes(x, y, group = level, linetype = level),
            colour = "black", inherit.aes = FALSE, linewidth = 0.5) +
  scale_colour_manual(values = class_colors, name = "Class") +
  scale_linetype_manual(
    values = c("1sigma" = "dotted", "2sigma" = "dashed",
               paste0(round(sqrt(MAH_THR), 2), "sigma") = "solid"),
    name = "Diploid null"
  ) +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(
    x        = "CS_RNA  = (DNA - RNA) × sign(DNA)",
    y        = "CS_Protein  = (DNA - Protein) × sign(DNA)",
    subtitle = sprintf(
      "Mahalanobis gate: chi2(df=2) p < %.2f  |  rho = %.2f  |  dashed lines = angular boundaries",
      ALPHA, rho),
    caption  = paste0(
      "Solid ellipse = 95% Mahalanobis boundary. ",
      "Dashed diagonals = theta = ±45° in whitened space (sd_P/sd_R = ",
      round(slope_45, 2), ").")
  )

ggsave(file.path(IMG_PATH, "cs_scatter_angular_class.pdf"),
       p_cs_scatter, width = 10, height = 8)
saveRDS(p_cs_scatter, file.path(IMG_PATH, "cs_scatter_angular_class.rds"))

# Also a polar-coordinate rose plot of theta, for genes outside the cloud
p_theta_rose <- df_classified %>%
  dplyr::filter(outside_cloud) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    theta_bin     = cut(theta, breaks = seq(-180, 180, by = 15),
                        include.lowest = TRUE)
  ) %>%
  dplyr::count(karyotype_lab, theta_bin) %>%
  ggplot(aes(x = theta_bin, y = n, fill = theta_bin)) +
  geom_col(show.legend = FALSE) +
  coord_polar(start = -pi / 2) +
  facet_wrap(~karyotype_lab) +
  theme_bw() +
  labs(x = "Angle theta (degrees)", y = "Number of genes",
       title = "Angular distribution of compensation direction",
       subtitle = "Genes outside the diploid cloud only",
       caption = "0° = pure RNA compensation. 90° = pure protein compensation (buffered).")

ggsave(file.path(IMG_PATH, "theta_rose_per_karyotype.pdf"),
       p_theta_rose, width = 9, height = 7)
saveRDS(p_theta_rose, file.path(IMG_PATH, "theta_rose_per_karyotype.rds"))

# ===========================================================================
# 4. Sanity check: LFC density
# ===========================================================================

df_long_lfc <- df_classified %>%
  dplyr::select(name, karyotype, RNA_lfc, Protein_lfc) %>%
  tidyr::pivot_longer(c(RNA_lfc, Protein_lfc),
                      names_to = "omic", values_to = "lfc") %>%
  dplyr::mutate(
    omic          = sub("_lfc", "", omic),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

p_lfc_density <- ggplot(df_long_lfc, aes(x = lfc, fill = omic, colour = omic)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~karyotype_lab, scales = "free_y") +
  scale_fill_manual(values = omic_colors) +
  scale_colour_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "log2 fold-change vs diploid", y = "Density")

ggsave(file.path(IMG_PATH, "lfc_density_vs_dosage.pdf"),
       p_lfc_density, width = 9, height = 6)
saveRDS(p_lfc_density, file.path(IMG_PATH, "lfc_density_vs_dosage.rds"))

# Also: CS score densities to verify they're centred at 0
df_long_cs <- df_classified %>%
  dplyr::select(name, karyotype, CS_RNA, CS_Protein) %>%
  tidyr::pivot_longer(c(CS_RNA, CS_Protein),
                      names_to = "omic", values_to = "cs") %>%
  dplyr::mutate(
    omic          = sub("CS_", "", omic),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

p_cs_density <- ggplot(df_long_cs, aes(x = cs, fill = omic, colour = omic)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~karyotype_lab, scales = "free_y") +
  scale_fill_manual(values = omic_colors) +
  scale_colour_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "Compensation Score (CS)", y = "Density",
       caption = "CS > 0 means attenuation relative to DNA. Should shift right in aneuploid karyotypes.")

ggsave(file.path(IMG_PATH, "cs_density_per_karyotype.pdf"),
       p_cs_density, width = 9, height = 6)
saveRDS(p_cs_density, file.path(IMG_PATH, "cs_density_per_karyotype.rds"))

# ===========================================================================
# 5. Genome-wide fraction per class per karyotype, with bootstrap CIs
# ===========================================================================

boot_fraction <- function(x, n_boot = 1000) {
  if (length(x) == 0) return(c(mean = NA, low = NA, high = NA))
  reps <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(reps),
    low  = unname(quantile(reps, 0.025)),
    high = unname(quantile(reps, 0.975)))
}

frac_by_karyotype <- df_classified %>%
  dplyr::select(karyotype, comp_RNA, comp_Protein, comp_Buffer) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein, comp_Buffer),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(level = sub("comp_", "", level)) %>%
  dplyr::group_by(karyotype, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_comp)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    level         = factor(level, levels = c("RNA", "Protein", "Buffer"))
  )

p_frac_compensated <- ggplot(frac_by_karyotype,
                             aes(x = karyotype_lab, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes",
       fill = "Compensation type",
       caption = paste0(
         "Angular classification. Gate: Mahalanobis chi2(df=2) p < ", ALPHA, ". ",
         "Buffer = Transcriptional + buffered class."))

ggsave(file.path(IMG_PATH, "fraction_compensated_per_karyotype.pdf"),
       p_frac_compensated, width = 7, height = 4)
saveRDS(p_frac_compensated,
        file.path(IMG_PATH, "fraction_compensated_per_karyotype.rds"))

# Class distribution stacked bar
class_dist <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    class_label   = factor(class_label, levels = names(class_colors))
  )

p_class_stack <- ggplot(class_dist,
                        aes(x = karyotype_lab, y = frac, fill = class_label)) +
  geom_col() +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes",
       fill = "Compensation class") +
  theme(legend.position = "right")

ggsave(file.path(IMG_PATH, "class_distribution_per_karyotype.pdf"),
       p_class_stack, width = 8, height = 5)
saveRDS(p_class_stack, file.path(IMG_PATH, "class_distribution_per_karyotype.rds"))

# ===========================================================================
# 6. Gain vs loss stratification
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
  dplyr::select(cn_direction, comp_RNA, comp_Protein, comp_Buffer) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein, comp_Buffer),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(level = sub("comp_", "", level)) %>%
  dplyr::group_by(cn_direction, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_comp)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(level = factor(level, levels = c("RNA", "Protein", "Buffer")))

p_frac_direction <- ggplot(frac_by_direction,
                           aes(x = cn_direction, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "Copy-number change", y = "Fraction compensated", fill = "Level")

ggsave(file.path(IMG_PATH, "fraction_compensated_gain_vs_loss.pdf"),
       p_frac_direction, width = 6, height = 4)
saveRDS(p_frac_direction, file.path(IMG_PATH, "fraction_compensated_gain_vs_loss.rds"))

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

if (nrow(gene_chrom) == 0) stop("No genes mapped to chromosomes.")

chrom_levels <- c(as.character(1:22), "X", "Y")

chrom_coverage <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(!is.na(chr)) %>%
  dplyr::group_by(chr, karyotype) %>%
  dplyr::summarise(
    n_genes      = dplyr::n(),
    n_rna        = sum(comp_RNA),
    n_prot       = sum(comp_Protein),
    n_buffer     = sum(comp_Buffer),
    frac_RNA     = mean(comp_RNA),
    frac_Protein = mean(comp_Protein),
    frac_Buffer  = mean(comp_Buffer),
    .groups = "drop"
  ) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_levels))

saveRDS(chrom_coverage, file.path(RES_PATH, "chrom_coverage.rds"))

# ---- 7.1 Coverage heatmap -------------------------------------------------
p_chrom_coverage <- chrom_coverage %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  ) %>%
  ggplot(aes(x = chr, y = karyotype_lab, fill = n_genes)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_genes), size = 2.5) +
  scale_fill_gradient(low = "white", high = "darkgreen", name = "n genes") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype", title = "Coverage per chromosome")

ggsave(file.path(IMG_PATH, "chrom_coverage_heatmap.pdf"),
       p_chrom_coverage, width = 11, height = 4)

# ---- 7.2 Fraction heatmap -------------------------------------------------
frac_by_chrom <- chrom_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::select(karyotype, chr, frac_RNA, frac_Protein, frac_Buffer) %>%
  tidyr::pivot_longer(c(frac_RNA, frac_Protein, frac_Buffer),
                      names_to = "level", values_to = "frac") %>%
  dplyr::mutate(
    level         = sub("frac_", "", level),
    level         = factor(level, levels = c("RNA", "Protein", "Buffer")),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

# ---- 7.3 Enrichment vs genome-wide baseline per karyotype -----------------
genome_wide <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    genome_RNA     = mean(comp_RNA),
    genome_Protein = mean(comp_Protein),
    genome_Buffer  = mean(comp_Buffer),
    .groups = "drop"
  )

enrichment_by_chrom <- chrom_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::left_join(genome_wide, by = "karyotype") %>%
  dplyr::mutate(
    enr_RNA     = log2((frac_RNA     + 0.01) / (genome_RNA     + 0.01)),
    enr_Protein = log2((frac_Protein + 0.01) / (genome_Protein + 0.01)),
    enr_Buffer  = log2((frac_Buffer  + 0.01) / (genome_Buffer  + 0.01))
  ) %>%
  dplyr::select(karyotype, chr, enr_RNA, enr_Protein, enr_Buffer) %>%
  tidyr::pivot_longer(c(enr_RNA, enr_Protein, enr_Buffer),
                      names_to = "level", values_to = "log2_enr") %>%
  dplyr::mutate(
    level         = sub("enr_", "", level),
    level         = factor(level, levels = c("RNA", "Protein", "Buffer")),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

saveRDS(enrichment_by_chrom, file.path(RES_PATH, "chrom_enrichment.rds"))

# ---- 7.4 Fisher tests per (karyotype, chromosome) -------------------------
sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ ""
  )
}

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
    dplyr::summarise(a = sum(is_comp), b = sum(!is_comp), .groups = "drop") %>%
    dplyr::group_by(karyotype) %>%
    dplyr::mutate(a_other = sum(a) - a, b_other = sum(b) - b) %>%
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
  fisher_per_chrom("comp_Buffer")
) %>%
  dplyr::mutate(
    level         = factor(level, levels = c("RNA", "Protein", "Buffer")),
    chr           = factor(chr,   levels = chrom_levels),
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
       title = "Raw fraction compensated per chromosome")

ggsave(file.path(IMG_PATH, "chrom_compensation_heatmap.pdf"),
       p_chrom_heatmap_raw, width = 11, height = 6)

p_chrom_heatmap_enr <- ggplot(enrichment_by_chrom,
                              aes(x = chr, y = karyotype_lab, fill = log2_enr)) +
  geom_tile(colour = "white") +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "log2 enr.") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Per-chromosome enrichment vs karyotype baseline")

ggsave(file.path(IMG_PATH, "chrom_enrichment_heatmap.pdf"),
       p_chrom_heatmap_enr, width = 11, height = 6)

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
                       midpoint = 0, name = "log2(OR)") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Fisher enrichment per chromosome")

ggsave(file.path(IMG_PATH, "chrom_fisher_heatmap.pdf"),
       p_chrom_heatmap_fisher, width = 11, height = 6)

# ===========================================================================
# 8. CORUM complex enrichment
# ===========================================================================

corum_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/CORUM_gene_list.txt"

if (file.exists(corum_path)) {
  corum_genes <- read.delim(corum_path) %>% dplyr::pull(GeneSym)

  df_classified <- df_classified %>%
    dplyr::mutate(in_complex = name %in% corum_genes)

  complex_enrichment <- df_classified %>%
    dplyr::group_by(karyotype, in_complex) %>%
    dplyr::summarise(
      frac_buffered = mean(comp_Buffer),
      frac_rna_comp = mean(comp_RNA),
      n = dplyr::n(),
      .groups = "drop"
    )

  fisher_buffered <- df_classified %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(
      p  = fisher.test(table(in_complex, comp_Buffer))$p.value,
      OR = fisher.test(table(in_complex, comp_Buffer))$estimate,
      .groups = "drop"
    ) %>%
    dplyr::mutate(padj = p.adjust(p, method = "BH"))

  saveRDS(complex_enrichment, file.path(RES_PATH, "complex_enrichment.rds"))
  saveRDS(fisher_buffered,    file.path(RES_PATH, "complex_fisher_buffered.rds"))

  complex_plot_df <- complex_enrichment %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
    )

  pvals_complex <- fisher_buffered %>%
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
      data        = pvals_complex,
      aes(xmin = xmin, xmax = xmax,
          annotations = annotation,
          y_position  = y_position,
          group       = karyotype_lab),
      manual      = TRUE, textsize = 3.2,
      tip_length  = 0.01, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "steelblue"),
                      labels = c("Non-complex", "CORUM complex")) +
    theme_bw() +
    labs(x = "Karyotype",
         y = "Fraction buffered (Transcriptional + buffered class)",
         fill = "",
         caption = "Fisher's exact test, BH-adjusted.")

  ggsave(file.path(IMG_PATH, "complex_buffering.pdf"), p_complex, width = 7, height = 4)
  saveRDS(p_complex, file.path(IMG_PATH, "complex_buffering.rds"))
}

# ===========================================================================
# 9. GO enrichment per compensation class
# ===========================================================================

genes_by_class <- df_classified %>%
  dplyr::filter(class_label %in% c("RNA compensated",
                                   "Transcriptional + buffered",
                                   "Protein only compensated",
                                   "Full dosage")) %>%
  dplyr::distinct(name, class_label) %>%
  split(.$class_label) %>%
  lapply(function(x) unique(x$name))

if (length(genes_by_class) > 1) {
  enrich_res <- compareCluster(
    geneClusters  = genes_by_class,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH"
  )

  dir.create(file.path(RES_PATH, "enrichment"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(enrich_res, file.path(RES_PATH, "enrichment", "compensation_class_GO.rds"))

  p_enrich <- dotplot(enrich_res, showCategory = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(IMG_PATH, "GO_by_class.pdf"), p_enrich, width = 10, height = 8)
  saveRDS(p_enrich, file.path(IMG_PATH, "GO_by_class.rds"))
}

# ===========================================================================
# 10. Chr 17 deep dive
# ===========================================================================

df_chr17 <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::mutate(is_chr17 = !is.na(chr) & chr == "17")

chr17_buffered <- df_chr17 %>%
  dplyr::filter(chr == "17", comp_Buffer) %>%
  dplyr::pull(name) %>% unique()

chr17_all <- df_chr17 %>%
  dplyr::filter(chr == "17") %>%
  dplyr::pull(name) %>% unique()

message(sprintf("chr 17 buffered: %d | chr 17 all: %d",
                length(chr17_buffered), length(chr17_all)))

# ---- 10.1 GO enrichment ---------------------------------------------------
go_by_ontology <- lapply(c("BP", "MF", "CC"), function(ont) {
  res <- enrichGO(
    gene          = chr17_buffered,
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
    dplyr::mutate(ontology = ont)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(ontology, p.adjust)

write_csv(go_by_ontology, file.path(RES_PATH, "GO_chr17_buffered_all_ontologies.csv"))
saveRDS(go_by_ontology,   file.path(RES_PATH, "GO_chr17_buffered_all_ontologies.rds"))

if (nrow(go_by_ontology) > 0) {
  top_terms <- go_by_ontology %>%
    dplyr::group_by(ontology) %>%
    dplyr::slice_min(p.adjust, n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log_padj = -log10(p.adjust))

  p_go_combined <- ggplot(top_terms,
                          aes(x = log_padj,
                              y = reorder(Description, log_padj),
                              colour = ontology, size = Count)) +
    geom_point() +
    theme_bw() +
    labs(x = "-log10(adj. p)", y = NULL,
         colour = "Ontology", size = "Gene count",
         title = "GO enrichment: chr 17 buffered genes vs chr 17 background")

  ggsave(file.path(IMG_PATH, "GO_chr17_buffered_combined.pdf"),
         p_go_combined, width = 9, height = 6)
  saveRDS(p_go_combined, file.path(IMG_PATH, "GO_chr17_buffered_combined.rds"))
}

# ---- 10.2 MitoCarta -------------------------------------------------------
mitocarta_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/Human.MitoCarta3.0.xls"

if (file.exists(mitocarta_path)) {
  mitocarta <- read_excel(mitocarta_path, sheet = "A Human MitoCarta3.0") %>%
    dplyr::as_tibble() %>%
    dplyr::pull(Symbol) %>%
    unique()

  run_mito_fisher <- function(scope_filter, scope_label) {
    df_classified %>%
      dplyr::left_join(gene_chrom, by = "name") %>%
      { if (scope_filter == "chr17") dplyr::filter(., chr == "17") else . } %>%
      dplyr::mutate(in_mito = name %in% mitocarta) %>%
      dplyr::group_by(karyotype) %>%
      dplyr::filter(dplyr::n() >= 10) %>%
      dplyr::summarise(
        a = sum( comp_Buffer &  in_mito),
        b = sum( comp_Buffer & !in_mito),
        c = sum(!comp_Buffer &  in_mito),
        d = sum(!comp_Buffer & !in_mito),
        .groups = "drop"
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        OR = fisher.test(matrix(c(a, b, c, d), nrow = 2))$estimate,
        p  = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(padj = p.adjust(p, method = "BH"), scope = scope_label)
  }

  mito_tests <- dplyr::bind_rows(
    run_mito_fisher("chr17",  "chr 17"),
    run_mito_fisher("genome", "Genome-wide")
  )
  saveRDS(mito_tests, file.path(RES_PATH, "mitocarta_enrichment.rds"))

  p_mito <- mito_tests %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
      log2_OR       = log2(pmax(OR, 0.1)),
      star          = sig_label(padj)
    ) %>%
    ggplot(aes(x = karyotype_lab, y = log2_OR, fill = scope)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = star, y = log2_OR + sign(log2_OR) * 0.15),
              position = position_dodge(width = 0.8), size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw() +
    labs(x = "Karyotype",
         y = "log2(OR) MitoCarta enrichment in buffered set",
         fill = "Scope")

  ggsave(file.path(IMG_PATH, "mitocarta_buffer_enrichment.pdf"),
         p_mito, width = 7, height = 4)
  saveRDS(p_mito, file.path(IMG_PATH, "mitocarta_buffer_enrichment.rds"))
}

# ===========================================================================
# 11. Summary tables
# ===========================================================================

summary_tbl <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class_label, values_from = n, values_fill = 0)

write_csv(summary_tbl, file.path(RES_PATH, "compensation_class_summary.csv"))

gate_summary <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    n_total              = dplyr::n(),
    n_outside_cloud      = sum(outside_cloud),
    n_full_dosage        = sum(class_label == "Full dosage"),
    n_rna_compensated    = sum(class_label == "RNA compensated"),
    n_trans_buffered     = sum(class_label == "Transcriptional + buffered"),
    n_prot_only          = sum(class_label == "Protein only compensated"),
    n_over_dosage        = sum(class_label == "Over-dosage"),
    frac_outside_cloud   = mean(outside_cloud),
    frac_comp_RNA        = mean(comp_RNA),
    frac_comp_Protein    = mean(comp_Protein),
    frac_comp_Buffer     = mean(comp_Buffer),
    median_theta_outside = median(theta[outside_cloud], na.rm = TRUE),
    .groups = "drop"
  )

write_csv(gate_summary, file.path(RES_PATH, "gate_summary_per_karyotype.csv"))
saveRDS(gate_summary,    file.path(RES_PATH, "gate_summary_per_karyotype.rds"))

cat("Done.\n")
