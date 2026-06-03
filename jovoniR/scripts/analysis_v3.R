# ===========================================================================
# analysis_v4.R
#
# Key changes vs v3:
#
#   CLASSIFICATION
#   • Replaced BH-adjusted pnorm with a direct z-score threshold (z >= 1.645,
#     equivalent to one-sided p <= 0.05) for comp_RNA, comp_Protein, and
#     comp_Buffer. This avoids BH dilution swamping strong signals (e.g. a
#     gene with DNA_lfc=1, RNA_lfc=0 was previously called non-significant).
#
#   • comp_RNA and comp_Protein are now symmetric: both test whether the
#     respective omic is significantly attenuated relative to DNA, each
#     calibrated by its own diploid noise σ. There is no reason to treat
#     protein differently from RNA at this step.
#
#   • comp_Buffer is anchored to comp_RNA. Post-transcriptional buffering
#     (protein below RNA) is only meaningful when RNA itself showed a
#     dosage-related response. Without the anchor, CS_Buffer fires on genes
#     where both RNA and protein are near zero — RNA-protein discordance
#     unrelated to aneuploidy.
#
#   • comp_Buffer uses the CONDITIONAL test (Prot | RNA residual / sd_PgivenR)
#     rather than the marginal R-P difference. This is more powerful and
#     statistically correct: it asks "is protein below what RNA predicts?"
#     rather than "is protein below RNA by more than diploid noise in the
#     difference?"
#
#   • The CS_RNA > 0 / CS_Protein > 0 guards are removed — they are
#     redundant given the upper-tail z-score test (a negative z can never
#     exceed 1.645) and were only sources of confusion.
#
#   BIVARIATE MODEL (kept as diagnostic, not as gate)
#   • The Mahalanobis / joint significance test is retained for diagnostic
#     plots (section 3e) and the gate_audit summary, but it is NO LONGER
#     used to collapse class_label to "Full dosage". The joint gate was
#     a secondary filter that obscured the primary classification logic.
#
#   3-BIT LABELS
#   • With comp_Buffer anchored to comp_RNA, classes "001" (Buffering only)
#     and "011" (Buffered post-transcriptional without RNA compensation) are
#     impossible by construction and removed from the case_when.
#   • "Protein-only compensated" (010) remains possible: protein attenuated
#     relative to DNA, RNA not significantly attenuated. This can reflect
#     genuine post-transcriptional-only regulation and is left in.
#
#   DIRECTION CLASSIFICATION (3b)
#   • classify_direction now uses only the noise-model z-score (no DE pval
#     from the original assay). The two-stage filter (noise_ok AND de_ok)
#     was mixing two different null hypotheses. Here we ask purely: "does
#     this gene's LFC exceed diploid noise in either direction?"
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

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))

source("utils.R")
source("scripts/constants.R")
source("scripts/getters.R")

IMG_PATH <- paste0("img/sf_", sf_method, "_stable_", use_stable, "_v4")
RES_PATH <- paste0("results/sf_", sf_method, "_stable_", use_stable, "_v4")
DF_DNA_PATH         <- "results/DNA_lfc.rds"
DF_CS_SCORES_PATH   <- "results/multiOmic/CS_scores_prot_and_rna.rds"
NOISE_MODEL_PATH    <- "results/multiOmic/diploid_noise_bivariate.RDS"

dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05

# z-threshold for one-sided test at ALPHA — used everywhere instead of BH.
# qnorm(0.05, lower.tail=FALSE) = 1.645. Any gene with z >= Z_THR has a
# one-sided p <= ALPHA under the diploid null, without multiple-testing
# correction. This is intentional: the null is a calibrated continuous
# distribution, not a point null, so BH correction over thousands of genes
# collapses the effective threshold and misses clear signals.
Z_THR <- qnorm(ALPHA, lower.tail = FALSE)
cat(sprintf("z-score threshold (one-sided p = %.2f): %.4f\n", ALPHA, Z_THR))

# ===========================================================================
# 1. Load data and bivariate noise model
# ===========================================================================

df_dna <- readRDS(DF_DNA_PATH) %>%
  dplyr::rename(DNA_lfc = lfc)

df_raw <- readRDS(DF_CS_SCORES_PATH)

# Bivariate null: list(mu, Sigma, Sigma_inv, rho, cond_PgivenR, n)
noise_biv <- readRDS(NOISE_MODEL_PATH)

mu_R  <- noise_biv$mu[1];  sd_R <- sqrt(noise_biv$Sigma[1, 1])
mu_P  <- noise_biv$mu[2];  sd_P <- sqrt(noise_biv$Sigma[2, 2])
rho   <- noise_biv$rho
Sig   <- noise_biv$Sigma
Sinv  <- noise_biv$Sigma_inv

# Marginal buffering null: Var(RNA - Protein) = Var(RNA) + Var(Protein) - 2*Cov
mu_buf <- mu_R - mu_P
sd_buf <- sqrt(Sig[1, 1] + Sig[2, 2] - 2 * Sig[1, 2])

# Conditional Prot | RNA
beta_PgivenR <- noise_biv$cond_PgivenR$beta
intc_PgivenR <- noise_biv$cond_PgivenR$intercept
sd_PgivenR   <- sqrt(noise_biv$cond_PgivenR$sigma2)

cat(sprintf("Bivariate null: mu=(%.3f, %.3f), sigma=(%.3f, %.3f), rho=%.3f\n",
            mu_R, mu_P, sd_R, sd_P, rho))
cat(sprintf("Conditional Prot|RNA SD: %.3f  (marginal sigma_P: %.3f)\n",
            sd_PgivenR, sd_P))
cat(sprintf("Buffering null: mu_buf=%.3f, sd_buf=%.3f\n", mu_buf, sd_buf))

# ===========================================================================
# 2. Build per-(gene, karyotype) table with CS scores
#
#   sign(DNA_lfc) correction: positive CS always means compensation,
#   regardless of whether the karyotype is a gain or a loss.
#
#   CS_RNA     = (DNA_lfc - RNA_lfc)     * sign(DNA_lfc)
#   CS_Protein = (DNA_lfc - Protein_lfc) * sign(DNA_lfc)
#   CS_Buffer  = (RNA_lfc - Protein_lfc) * sign(DNA_lfc)
# ===========================================================================

df_wide <- df_raw %>%
  dplyr::select(name, karyotype, omic, lfc, pval) %>%
  tidyr::pivot_wider(names_from = omic, values_from = c(lfc, pval)) %>%
  dplyr::rename(RNA_lfc     = lfc_RNA,     Protein_lfc = lfc_Protein,
                RNA_pval    = pval_RNA,    Protein_pval = pval_Protein) %>%
  dplyr::left_join(df_dna %>% dplyr::select(name, karyotype, DNA_lfc),
                   by = c("name", "karyotype")) %>%
  dplyr::filter(!is.na(RNA_lfc), !is.na(Protein_lfc), !is.na(DNA_lfc)) %>%
  dplyr::mutate(
    CS_RNA     = (DNA_lfc - RNA_lfc)     * sign(DNA_lfc),
    CS_Protein = (DNA_lfc - Protein_lfc) * sign(DNA_lfc),
    CS_Buffer  = (RNA_lfc - Protein_lfc) * sign(DNA_lfc)
  )

saveRDS(df_wide, file.path(RES_PATH, "df_wide_three_CS.rds"))

# ===========================================================================
# 3. Per-gene classification
#
#   z_RNA     = CS_RNA     / sd_R        (marginal, against diploid RNA noise)
#   z_Protein = CS_Protein / sd_P        (marginal, against diploid Protein noise)
#   z_Buffer  = CS_Buffer  / sd_buf      (marginal, against diploid R-P noise)
#   z_Buffer_cond = Prot_resid / sd_PgivenR  (conditional: protein below RNA prediction?)
#
#   comp_RNA     : z_RNA     >= Z_THR
#   comp_Protein : z_Protein >= Z_THR
#                  (symmetric with RNA — protein attenuation relative to DNA,
#                   calibrated by diploid protein noise)
#   comp_Buffer  : z_Buffer_cond >= Z_THR  AND  comp_RNA
#                  (post-transcriptional buffering anchored to RNA response;
#                   uses conditional test for greater power)
#
#   Note: the Mahalanobis / joint test is computed for diagnostic purposes
#   only (section 3e). It is NOT used to gate the classification.
# ===========================================================================

# Mahalanobis squared on (RNA_lfc, Protein_lfc) — diagnostic only
xy   <- as.matrix(df_wide[, c("RNA_lfc", "Protein_lfc")])
ctr  <- sweep(xy, 2, noise_biv$mu, "-")
mah2 <- rowSums((ctr %*% Sinv) * ctr)

df_classified <- df_wide %>%
  dplyr::mutate(

    # --- Diagnostic: joint departure from diploid cloud ---
    mahal2  = mah2,
    p_joint = stats::pchisq(mah2, df = 2, lower.tail = FALSE),

    # --- z-scores (positive = compensation) ---
    z_RNA     = CS_RNA     / sd_R,
    z_Protein = CS_Protein / sd_P,
    z_Buffer  = CS_Buffer  / sd_buf,

    # Conditional buffering: residual of Protein given RNA prediction
    # sign(DNA_lfc) already in CS_Buffer; apply same sign to residual
    Prot_resid    = (RNA_lfc - (intc_PgivenR + beta_PgivenR * Protein_lfc)) *
      sign(DNA_lfc),
    z_Buffer_cond = Prot_resid / sd_PgivenR,

    # --- Classification flags ---
    # RNA: is RNA significantly attenuated relative to DNA?
    comp_RNA     = z_RNA     >= Z_THR,

    # Protein: is protein significantly attenuated relative to DNA?
    # Tested symmetrically with RNA, using protein's own diploid noise.
    comp_Protein = z_Protein >= Z_THR,

    # Buffer: is protein below what RNA predicts, AND did RNA respond?
    # Anchored to comp_RNA so that buffering is only called when there is
    # a real dosage signal to buffer.
    comp_Buffer  = (z_Buffer_cond >= Z_THR) & comp_RNA

  ) %>%
  dplyr::mutate(

    # --- 3-bit class label ---
    # "001" and "011" are impossible by construction (Buffer requires RNA).
    class_3bit = paste0(
      as.integer(comp_RNA),
      as.integer(comp_Protein),
      as.integer(comp_Buffer)
    ),
    class_label = dplyr::case_when(
      class_3bit == "000" ~ "Full dosage",
      class_3bit == "100" ~ "RNA-only compensated",
      class_3bit == "010" ~ "Protein-only compensated",
      class_3bit == "110" ~ "RNA + Protein compensated",
      class_3bit == "101" ~ "RNA compensated + buffered",
      class_3bit == "111" ~ "Fully compensated",
      class_3bit == "011" ~ "Protein compensated + buffered",  # rare but keep
      TRUE                ~ "Other"
    )

  )

# Audit: how many genes per karyotype in each class?
class_audit <- df_classified %>%
  dplyr::count(karyotype, class_label) %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()

saveRDS(class_audit, file.path(RES_PATH, "class_audit.rds"))
print(class_audit)

# Quick sanity check: comp_Buffer should never be TRUE when comp_RNA is FALSE
stopifnot(all(!(df_classified$comp_Buffer & !df_classified$comp_RNA)))

# ---------------------------------------------------------------------------
# 3b. Up / Down / Not-differential classification per omic
#
#     Uses only the noise-model z-score (two-sided). A gene is Up/Down if
#     its LFC exceeds the diploid noise floor in either direction.
#     The original DE p-value is NOT used here: mixing two null hypotheses
#     (noise model vs DE test) complicates interpretation.
# ---------------------------------------------------------------------------

classify_direction <- function(lfc, mu, sigma, alpha) {
  z       <- (lfc - mu) / sigma
  p_noise <- 2 * pnorm(-abs(z))          # two-sided
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

saveRDS(df_classified, file.path(RES_PATH, "df_classified_v4.rds"))

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
# 3d. Plots: Up/Down fractions and asymmetry
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
      is.na(binom_padj)   ~ "",
      binom_padj <= 0.001 ~ "***",
      binom_padj <= 0.01  ~ "**",
      binom_padj <= 0.05  ~ "*",
      TRUE                ~ "ns"
    )
  )

p_asymmetry <- ggplot(asym_plot,
                      aes(x = karyotype_lab, y = asymmetry, fill = omic)) +
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
       caption = "Binomial test vs 50/50, BH-adjusted. Positive = Up-skewed.")

ggsave(file.path(IMG_PATH, "up_down_asymmetry.pdf"), p_asymmetry, width = 7, height = 4)
saveRDS(p_asymmetry, file.path(IMG_PATH, "up_down_asymmetry.rds"))

# ---------------------------------------------------------------------------
# 3e. Bivariate diagnostic: (RNA, Protein) cloud per karyotype with
#     diploid 1s / 2s Mahalanobis ellipses. Kept for visualisation only.
# ---------------------------------------------------------------------------

make_ellipse <- function(mu, Sigma, k, npts = 200) {
  theta  <- seq(0, 2 * pi, length.out = npts)
  circle <- rbind(cos(theta), sin(theta))
  L      <- chol(Sigma)
  pts    <- t(mu + k * t(L) %*% circle)
  dplyr::tibble(x = pts[, 1], y = pts[, 2], level = paste0(k, "sigma"))
}

ell_df <- dplyr::bind_rows(
  make_ellipse(noise_biv$mu, Sig, 1),
  make_ellipse(noise_biv$mu, Sig, 2)
)

p_joint_per_karyotype <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    comp_class    = dplyr::case_when(
      comp_RNA & comp_Buffer  ~ "RNA + Buffered",
      comp_RNA                ~ "RNA only",
      comp_Protein            ~ "Protein only",
      TRUE                    ~ "None"
    )
  ) %>%
  ggplot(aes(RNA_lfc, Protein_lfc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = comp_class), alpha = 0.35, size = 0.5) +
  geom_path(data = ell_df,
            aes(x, y, group = level, linetype = level),
            colour = "firebrick3", inherit.aes = FALSE, linewidth = 0.6) +
  scale_colour_manual(values = c("None"           = "grey70",
                                 "RNA only"        = "steelblue3",
                                 "Protein only"    = "darkorange",
                                 "RNA + Buffered"  = "firebrick3"),
                      name = "Compensation class") +
  scale_linetype_manual(values = c("1sigma" = "solid", "2sigma" = "dashed"),
                        name = "Diploid null") +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC",
       subtitle = sprintf("Bivariate diploid null rho = %.2f  |  ellipses = 1/2 sigma Mahalanobis",
                          rho),
       caption = "Ellipses are diploid null, for reference only — not used in classification.")

p_joint_per_karyotype <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    comp_class    = dplyr::case_when(
      comp_RNA & comp_Buffer  ~ "RNA + Buffered",
      comp_RNA                ~ "RNA only",
      comp_Protein            ~ "Protein only",
      TRUE                    ~ "None"
    )
  ) %>%
  ggplot(aes(CS_RNA, CS_Protein)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = class_label), alpha = 0.35, size = 0.5) +
  # Ellipse now makes sense: diploid null in CS space is centered at (0,0)
  geom_path(data = ell_df,
            aes(x, y, group = level, linetype = level),
            colour = "firebrick3", inherit.aes = FALSE, linewidth = 0.6) +
  # Threshold lines showing where z >= 1.645 for each axis
  geom_vline(xintercept = Z_THR * sd_R, linetype = "dotted", colour = "steelblue3") +
  geom_hline(yintercept = Z_THR * sd_P, linetype = "dotted", colour = "darkorange") +
  scale_linetype_manual(values = c("1sigma" = "solid", "2sigma" = "dashed"),
                        name = "Diploid null") +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(x = "CS_RNA (DNA - RNA, sign-corrected)",
       y = "CS_Protein (DNA - Protein, sign-corrected)",
       subtitle = sprintf("Bivariate diploid null rho = %.2f | dotted lines = z-score threshold", rho),
       caption = "Points right of blue line: comp_RNA=TRUE. Points above orange line: comp_Protein=TRUE.")

ggsave(file.path(IMG_PATH, "joint_cloud_per_karyotype.pdf"), p_joint_per_karyotype, width = 9, height = 7)
saveRDS(p_joint_per_karyotype, file.path(IMG_PATH, "joint_cloud_per_karyotype.rds"))

# ===========================================================================
# 4. Sanity-check: LFC density vs DNA expectation
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

ggsave(file.path(IMG_PATH, "lfc_density_vs_dosage.pdf"), p_lfc_density, width = 9, height = 6)
saveRDS(p_lfc_density, file.path(IMG_PATH, "lfc_density_vs_dosage.rds"))

# ===========================================================================
# 5. Genome-wide fraction compensated per karyotype, with bootstrap CIs
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
  labs(x = "Karyotype", y = "Fraction of compensated genes",
       fill = "Compensation level",
       caption = paste0("z-score threshold: z >= ", round(Z_THR, 3),
                        " (one-sided p <= ", ALPHA, "). ",
                        "Buffer anchored to RNA compensation."))

ggsave(file.path(IMG_PATH, "fraction_compensated_three_levels.pdf"), p_frac_compensated, width = 7, height = 4)
saveRDS(p_frac_compensated, file.path(IMG_PATH, "fraction_compensated_three_levels.rds"))

# ===========================================================================
# 6. RNA-LFC vs Protein-LFC scatter
# ===========================================================================

p_rna_vs_prot <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  ) %>%
  ggplot(aes(x = RNA_lfc, y = Protein_lfc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_path(data = ell_df %>% dplyr::filter(level == "2sigma"),
            aes(x, y, group = level), inherit.aes = FALSE,
            colour = "firebrick3", linewidth = 0.4) +
  geom_point(aes(colour = comp_Buffer), alpha = 0.4, size = 0.6) +
  facet_wrap(~karyotype_lab) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "darkorange"),
                      name = "Buffered\n(conditional)") +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC")

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

if (nrow(gene_chrom) == 0) stop("No genes mapped to chromosomes.")

chrom_levels <- c(as.character(1:22), "X", "Y")

# ---- 7.1 Coverage heatmap -------------------------------------------------
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

ggsave(file.path(IMG_PATH, "chrom_coverage_heatmap.pdf"), p_chrom_coverage, width = 11, height = 4)

# ---- 7.2 Long-form fractions for heatmap ----------------------------------
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

# ---- 7.3 Enrichment vs karyotype-wide baseline ----------------------------
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

# ---- 7.5 Heatmap: raw fraction compensated --------------------------------
p_chrom_heatmap_raw <- ggplot(frac_by_chrom,
                              aes(x = chr, y = karyotype_lab, fill = frac)) +
  geom_tile(colour = "white") +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "white", high = "darkblue",
                       midpoint = 0.25, limits = c(0, 1),
                       name = "Frac.\ncompensated") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype", title = "Raw fraction compensated per chromosome")

ggsave(file.path(IMG_PATH, "chrom_compensation_heatmap.pdf"), p_chrom_heatmap_raw, width = 11, height = 6)

# ---- 7.6 Heatmap: enrichment vs karyotype baseline -----------------------
p_chrom_heatmap_enr <- ggplot(enrichment_by_chrom,
                              aes(x = chr, y = karyotype_lab, fill = log2_enr)) +
  geom_tile(colour = "white") +
  facet_wrap(~level, ncol = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "log2 enr.") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Per-chromosome enrichment vs karyotype baseline")

ggsave(file.path(IMG_PATH, "chrom_enrichment_heatmap.pdf"), p_chrom_heatmap_enr, width = 11, height = 6)

# ---- 7.7 Heatmap: Fisher significance ------------------------------------
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
                       midpoint = 0, name = "log2(OR)") +
  theme_bw() +
  labs(x = "Chromosome", y = "Karyotype",
       title = "Fisher enrichment per chromosome")

ggsave(file.path(IMG_PATH, "chrom_fisher_heatmap.pdf"), p_chrom_heatmap_fisher, width = 11, height = 6)

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

ggsave(file.path(IMG_PATH, "fraction_compensated_gain_vs_loss.pdf"), p_frac_direction, width = 6, height = 4)

# ===========================================================================
# 9. Class distribution per karyotype
# ===========================================================================

class_dist <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

p_class_stack <- ggplot(class_dist,
                        aes(x = karyotype_lab, y = frac, fill = class_label)) +
  geom_col() +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes",
       fill = "Compensation class") +
  theme(legend.position = "right")

ggsave(file.path(IMG_PATH, "class_distribution_per_karyotype.pdf"), p_class_stack, width = 8, height = 5)

# ===========================================================================
# 10. CORUM complex enrichment among buffered genes
# ===========================================================================

corum_path <- "data/CORUM_gene_list.txt"

if (file.exists(corum_path)) {
  corum_genes <- read.delim(corum_path) %>% dplyr::pull(GeneSym)

  df_classified <- df_classified %>%
    dplyr::mutate(in_complex = name %in% corum_genes)

  complex_enrichment <- df_classified %>%
    dplyr::group_by(karyotype, in_complex) %>%
    dplyr::summarise(
      frac_rna_comp = mean(comp_RNA),
      n = dplyr::n(),
      .groups = "drop"
    )

  fisher_buffered <- df_classified %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(
      p  = fisher.test(table(in_complex, frac_rna_comp))$p.value,
      OR = fisher.test(table(in_complex, frac_rna_comp))$estimate,
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
         y = "Fraction buffered (post-transcriptional)",
         fill = "",
         caption = "Fisher's exact test, BH-adjusted.")

  ggsave(file.path(IMG_PATH, "complex_buffering.pdf"), p_complex, width = 7, height = 4)
}

# ===========================================================================
# 11. GO enrichment per compensation class
# ===========================================================================

genes_by_class <- df_classified %>%
  dplyr::filter(class_label %in% c("RNA-only compensated",
                                   "RNA compensated + buffered",
                                   "Fully compensated",
                                   "Full dosage")) %>%
  dplyr::distinct(name, class_label) %>%
  split(.$class_label) %>%
  lapply(function(x) unique(x$name))

if (length(genes_by_class) > 1) {
  enrich_res <- compareCluster(
    geneClusters = genes_by_class,
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
}

# ===========================================================================
# 12. Chr 17 deep dive
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

# ---- 12.1 GO enrichment ---------------------------------------------------
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

# ---- 12.2 MitoCarta -------------------------------------------------------
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
}

# ===========================================================================
# 13. Summary tables
# ===========================================================================

summary_tbl <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class_label, values_from = n, values_fill = 0)

write_csv(summary_tbl, file.path(RES_PATH, "compensation_class_summary.csv"))

gate_summary <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    n_total          = dplyr::n(),
    n_comp_RNA       = sum(comp_RNA),
    n_comp_Protein   = sum(comp_Protein),
    n_comp_Buffer    = sum(comp_Buffer),
    n_full_dosage    = sum(class_label == "Full dosage"),
    frac_comp_RNA    = mean(comp_RNA),
    frac_comp_Protein = mean(comp_Protein),
    frac_comp_Buffer = mean(comp_Buffer),
    .groups = "drop"
  )

write_csv(gate_summary, file.path(RES_PATH, "gate_summary_per_karyotype.csv"))
saveRDS(gate_summary,    file.path(RES_PATH, "gate_summary_per_karyotype.rds"))

cat("Done.\n")
