
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(readxl)

source("utils.R")
source("scripts/constants.R")
source("scripts/getters.R")

# ── Annotation for significance brackets ─────────────────────────────────────
sig_label <- function(p) dplyr::case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""
)

IMG_PATH <- paste0("img/analysis_sf_", SF_METHOD, "_stable_", USE_STABLE)
RES_PATH <- paste0("results/analysis_sf_", SF_METHOD, "_stable_", USE_STABLE)
DF_DNA_PATH         <- "results/DNA_lfc.rds"
DF_CS_SCORES_PATH   <- "results/multiOmic/CS_scores_prot_and_rna.rds"
NOISE_MODEL_PATH    <- "results/multiOmic/diploid_noise.RDS"
NOISE_MODEL_BIV_PATH    <- "results/multiOmic/diploid_noise_bivariate.RDS"
dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load data and bivariate noise model
# ===========================================================================

df_dna <- readRDS(DF_DNA_PATH) %>%
  dplyr::rename(DNA_lfc = lfc)

df_raw <- readRDS(DF_CS_SCORES_PATH)

# Bivariate null: list(mu, Sigma, Sigma_inv, rho, cond_PgivenR, n)
noise     <- readRDS(NOISE_MODEL_PATH)      # diploid_noise tibble: mu, sigma, nu per omic
noise_biv <- readRDS(NOISE_MODEL_BIV_PATH) # diploid_noise_biv: Sigma, Sigma_inv, nu, rho

# Marginal parameters — always from noise (the t-fitted marginals)
mu_R  <- noise$mu[noise$quantity == "CS_RNA"]
sd_R  <- noise$sigma[noise$quantity == "CS_RNA"]
nu_R  <- noise$nu[noise$quantity == "CS_RNA"]

mu_P  <- noise$mu[noise$quantity == "CS_Protein"]
sd_P  <- noise$sigma[noise$quantity == "CS_Protein"]
nu_P  <- noise$nu[noise$quantity == "CS_Protein"]

# Thresholds — from the correct per-omic df
Z_THR_RNA  <- qt(ALPHA, df = nu_R, lower.tail = FALSE)
Z_THR_PROT <- qt(ALPHA, df = nu_P, lower.tail = FALSE)

# Joint parameters — only from noise_biv, only for Mahalanobis
Sinv     <- noise_biv$Sigma_inv
mu_biv   <- noise_biv$mu
nu_joint <- noise_biv$nu

cat(sprintf("z-threshold RNA  (df=%.1f, p=%.2f): %.4f\n", nu_R, ALPHA, Z_THR_RNA))
cat(sprintf("z-threshold Prot (df=%.1f, p=%.2f): %.4f\n", nu_P, ALPHA, Z_THR_PROT))
cat(sprintf("Marginal null: mu=(%.3f, %.3f), sigma=(%.3f, %.3f), rho=%.3f\n",
            mu_R, mu_P, sd_R, sd_P, noise_biv$rho))

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
  dplyr::rename(RNA_lfc      = lfc_RNA,      Protein_lfc  = lfc_Protein,
                RNA_pval     = pval_RNA,     Protein_pval = pval_Protein) %>%
  dplyr::left_join(df_dna %>% dplyr::select(name, karyotype, DNA_lfc),
                   by = c("name", "karyotype")) %>%
  dplyr::filter(!is.na(RNA_lfc), !is.na(Protein_lfc), !is.na(DNA_lfc)) %>%
  dplyr::mutate(
    CS_RNA     = (DNA_lfc - RNA_lfc)     * sign(DNA_lfc),
    CS_Protein = (DNA_lfc - Protein_lfc) * sign(DNA_lfc)
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
ctr  <- sweep(xy, 2, mu_biv, "-")
mah2 <- rowSums((ctr %*% Sinv) * ctr)

df_classified <- df_wide %>%
  dplyr::mutate(
    
    # Diagnostic: joint departure — p-value from F(2, nu_joint), correct for bivariate-t
    mahal2  = mah2,
    p_joint = pf(mah2 / 2, df1 = 2, df2 = nu_joint, lower.tail = FALSE),
    
    # z-scores: centre and scale by marginal t parameters
    z_RNA     = (CS_RNA     - mu_R) / sd_R,
    z_Protein = (CS_Protein - mu_P) / sd_P,
    
    # one-sided p-values under the fitted marginal t
    p_RNA     = pt(z_RNA,     df = nu_R, lower.tail = FALSE),
    p_Protein = pt(z_Protein, df = nu_P, lower.tail = FALSE),
    
    # Classification flags — omic-specific thresholds
    comp_RNA     = z_RNA     >= Z_THR_RNA,
    comp_Protein = z_Protein >= Z_THR_PROT
    
  ) %>%
  dplyr::mutate(
    
    # 2-bit class label
    class_2bit = paste0(as.integer(comp_RNA), as.integer(comp_Protein)),
    class_label = dplyr::case_when(
      class_2bit == "00" ~ "Full dosage sensitive",
      class_2bit == "10" ~ "RNA compensated",
      class_2bit == "01" ~ "Protein compensated",
      class_2bit == "11" ~ "Fully compensated",
      TRUE               ~ "Other"
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

# ---------------------------------------------------------------------------
# 3b. Up / Down / Not-differential classification per omic
#
#     Uses only the noise-model z-score (two-sided). A gene is Up/Down if
#     its LFC exceeds the diploid noise floor in either direction.
#     The original DE p-value is NOT used here: mixing two null hypotheses
#     (noise model vs DE test) complicates interpretation.
# ---------------------------------------------------------------------------

classify_direction_t <- function(lfc, mu, sigma, nu, alpha) {
  z       <- (lfc - mu) / sigma
  p_noise <- 2 * pt(-abs(z), df = nu)
  dplyr::case_when(
    p_noise > alpha ~ "Not differential",
    lfc > mu        ~ "Up",
    lfc < mu        ~ "Down",
    TRUE            ~ "Not differential"
  )
}

df_classified <- df_classified %>%
  dplyr::mutate(
    dir_RNA     = classify_direction_t(RNA_lfc,     mu_R, sd_R, nu_R, ALPHA),
    dir_Protein = classify_direction_t(Protein_lfc, mu_P, sd_P, nu_P, ALPHA)
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

make_t_contour <- function(mu, Sigma, nu, prob, npts = 300) {
  r     <- sqrt(2 * qf(prob, df1 = 2, df2 = nu))
  theta <- seq(0, 2 * pi, length.out = npts)
  L     <- chol(Sigma)
  pts   <- t(mu + r * t(L) %*% rbind(cos(theta), sin(theta)))
  dplyr::tibble(x = pts[, 1], y = pts[, 2],
                level = paste0(round(prob * 100), "% density"))
}

contour_df <- dplyr::bind_rows(
  make_t_contour(noise_biv$mu, noise_biv$Sigma, nu_joint, 0.50),
  make_t_contour(noise_biv$mu, noise_biv$Sigma, nu_joint, 0.90),
  make_t_contour(noise_biv$mu, noise_biv$Sigma, nu_joint, 0.95)
)

p_lfc_space <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping),
    comp_class    = dplyr::case_when(
      comp_RNA & comp_Protein ~ "RNA + Protein",
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
  geom_path(data = contour_df,
            aes(x, y, group = level, linetype = level),
            colour = "firebrick3", inherit.aes = FALSE, linewidth = 0.6) +
  scale_colour_manual(
    values = c("None"          = "grey70",
               "RNA only"      = "steelblue3",
               "Protein only"  = "darkorange",
               "RNA + Protein" = "firebrick3"),
    name = "Compensation class"
  ) +
  scale_linetype_manual(
    values = c("50% density" = "solid",
               "90% density" = "dashed",
               "95% density" = "dotted"),
    name = "Bivariate-t null"
  ) +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC")

contour_df_cs <- dplyr::bind_rows(
  make_t_contour(c(mu_R, mu_P), noise_biv$Sigma, nu_joint, 0.50),
  make_t_contour(c(mu_R, mu_P), noise_biv$Sigma, nu_joint, 0.90),
  make_t_contour(c(mu_R, mu_P), noise_biv$Sigma, nu_joint, 0.95)
)

p_cs_space <- df_classified %>%
  dplyr::mutate(
    karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping)
  ) %>%
  ggplot(aes(CS_RNA, CS_Protein)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = class_label), alpha = 0.35, size = 0.5) +
  geom_path(data = contour_df_cs,
            aes(x, y, group = level, linetype = level),
            colour = "firebrick3", inherit.aes = FALSE, linewidth = 0.6) +
  # Threshold lines: Z_THR * sigma + mu (t-threshold, not Gaussian)
  geom_vline(xintercept = Z_THR_RNA  * sd_R + mu_R,
             linetype = "dotted", colour = "steelblue3") +
  geom_hline(yintercept = Z_THR_PROT * sd_P + mu_P,
             linetype = "dotted", colour = "darkorange") +
  scale_colour_manual(
    values = c("Full dosage sensitive" = "grey70",
               "RNA compensated"       = "steelblue3",
               "Protein compensated"   = "darkorange",
               "Fully compensated"     = "firebrick3"),
    name = "Class"
  ) +
  scale_linetype_manual(
    values = c("50% density" = "solid",
               "90% density" = "dashed",
               "95% density" = "dotted"),
    name = "Bivariate-t null"
  ) +
  facet_wrap(~karyotype_lab) +
  coord_fixed() +
  theme_bw() +
  labs(x = "CS_RNA (DNA - RNA, sign-corrected)",
       y = "CS_Protein (DNA - Protein, sign-corrected)")

ggsave(file.path(IMG_PATH, "joint_cloud_per_karyotype_LFC.pdf"), p_lfc_space, width = 9, height = 7)
ggsave(file.path(IMG_PATH, "joint_cloud_per_karyotype_CS.pdf"), p_cs_space, width = 9, height = 7)
saveRDS(p_lfc_space, file.path(IMG_PATH, "joint_cloud_per_karyotype_LFC.rds"))
saveRDS(p_cs_space, file.path(IMG_PATH, "joint_cloud_per_karyotype_CS.rds"))

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
  dplyr::select(karyotype, comp_RNA, comp_Protein) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein),
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

frac_dosage_sensitive <- df_classified %>%
  dplyr::select(karyotype, comp_RNA, comp_Protein) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(
    level    = sub("comp_", "", level),
    is_dosage = !is_comp          # <── invert HERE, before bootstrapping
  ) %>%
  dplyr::group_by(karyotype, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_dosage)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    level         = factor(level, levels = c("RNA", "Protein", "Buffer"))
  )

pval_comp <- df_classified %>%
  dplyr::select(karyotype, comp_RNA, comp_Protein) %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    p_value = {
      n  <- dplyr::n()
      x1 <- sum(comp_RNA,     na.rm = TRUE)
      x2 <- sum(comp_Protein, na.rm = TRUE)
      prop.test(c(x1, x2), c(n, n), correct = FALSE)$p.value
    },
    y_pos = max(
      mean(comp_RNA,     na.rm = TRUE),
      mean(comp_Protein, na.rm = TRUE)
    ) + 0.06,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

pval_ds <- df_classified %>%
  dplyr::select(karyotype, comp_RNA, comp_Protein) %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    p_value = {
      n  <- dplyr::n()
      x1 <- sum(!comp_RNA,     na.rm = TRUE)   # <── NOT compensated
      x2 <- sum(!comp_Protein, na.rm = TRUE)
      prop.test(c(x1, x2), c(n, n), correct = FALSE)$p.value
    },
    y_pos = max(
      mean(!comp_RNA,     na.rm = TRUE),
      mean(!comp_Protein, na.rm = TRUE)
    ) + 0.06,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )


p_frac_compensated <- ggplot(frac_by_karyotype,
                             aes(x = karyotype_lab, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_text(data = pval_comp,
            aes(x = karyotype_lab, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = omic_colors) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of compensated genes",
       fill = "Compensation level")

p_frac_dosage_sensitive <- ggplot(frac_dosage_sensitive,
                                  aes(x = karyotype_lab, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_text(data = pval_ds,
            aes(x = karyotype_lab, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = omic_colors) +
  # scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of dosage-sensitive genes",
       fill = "Compensation level")

ggsave(file.path(IMG_PATH, "fraction_compensated_three_levels.pdf"),
       p_frac_compensated, width = 7, height = 4)
saveRDS(p_frac_compensated,
        file.path(IMG_PATH, "fraction_compensated_three_levels.rds"))

ggsave(file.path(IMG_PATH, "fraction_dosage_sensitive_three_levels.pdf"),
       p_frac_dosage_sensitive, width = 7, height = 4)
saveRDS(p_frac_dosage_sensitive,
        file.path(IMG_PATH, "fraction_dosage_sensitive_three_levels.rds"))

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
  geom_point(alpha = 0.4, size = 0.6, colour = "grey50") +
  geom_path(data = contour_df,
            aes(x, y, group = level), inherit.aes = FALSE,
            colour = "firebrick3", linewidth = 0.4) +
  facet_wrap(~karyotype_lab) +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC")

ggsave(file.path(IMG_PATH, "rna_vs_protein_lfc.pdf"), p_rna_vs_prot, width = 9, height = 6)
saveRDS(p_rna_vs_prot, file.path(IMG_PATH, "rna_vs_protein_lfc.rds"))

# ===========================================================================
# 7. Per-chromosome-arm compensation analysis
# ===========================================================================

# ---- 7.0 Fetch chromosome arm assignments via biomaRt ---------------------
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_chrom <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "chromosome_name",
                 "band", "start_position"),
  filters    = "hgnc_symbol",
  values     = unique(df_classified$name),
  mart       = mart
) %>%
  dplyr::as_tibble() %>%
  dplyr::rename(name = hgnc_symbol, chr = chromosome_name) %>%
  dplyr::filter(chr %in% c(as.character(1:22), "X", "Y")) %>%
  dplyr::mutate(
    # Cytogenetic band starts with p or q → chromosome arm
    arm     = dplyr::case_when(
      stringr::str_starts(band, "p") ~ paste0(chr, "p"),
      stringr::str_starts(band, "q") ~ paste0(chr, "q"),
      TRUE                           ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(arm)) %>%
  dplyr::distinct(name, .keep_all = TRUE)

n_total  <- length(unique(df_classified$name))
n_mapped <- nrow(gene_chrom)
message(sprintf("Chromosome arm mapping: %d / %d genes mapped (%.1f%%)",
                n_mapped, n_total, 100 * n_mapped / n_total))

if (nrow(gene_chrom) == 0) stop("No genes mapped to chromosome arms.")

# Ordered arm levels: 1p, 1q, 2p, 2q, ..., 22p, 22q, Xp, Xq
arm_levels <- purrr::map(
  c(as.character(1:22), "X", "Y"),
  ~ paste0(.x, c("p", "q"))
) %>%
  unlist() %>%
  intersect(unique(gene_chrom$arm))  # keep only arms present in data

saveRDS(gene_chrom, file.path(RES_PATH, "gene_chrom_arm.rds"))

# ---- 7.1 Coverage per arm -------------------------------------------------
arm_coverage <- df_classified %>%
  dplyr::left_join(gene_chrom %>%
                     dplyr::select(name, chr, arm),
                   by = "name") %>%
  dplyr::filter(!is.na(arm)) %>%
  dplyr::group_by(arm, chr, karyotype) %>%
  dplyr::summarise(
    n_genes        = dplyr::n(),
    # Compensated counts
    n_rna_comp     = sum(comp_RNA),
    n_prot_comp    = sum(comp_Protein),
    # Dosage sensitive = NOT compensated at that level
    n_rna_ds       = sum(!comp_RNA),
    n_prot_ds      = sum(!comp_Protein),
    # Fractions
    frac_RNA_comp  = mean(comp_RNA),
    frac_Prot_comp = mean(comp_Protein),
    frac_RNA_ds    = mean(!comp_RNA),
    frac_Prot_ds   = mean(!comp_Protein),
    .groups        = "drop"
  ) %>%
  dplyr::mutate(
    arm           = factor(arm, levels = arm_levels),
    karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping)
  )

saveRDS(arm_coverage, file.path(RES_PATH, "arm_coverage.rds"))

# ---- 7.2 Coverage heatmap (n genes per arm) --------------------------------
p_arm_coverage <- arm_coverage %>%
  ggplot(aes(x = arm, y = karyotype_lab, fill = n_genes)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_genes), size = 2) +
  scale_fill_gradient(low = "white", high = "darkgreen", name = "n genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
  labs(x = "Chromosome arm", y = "Karyotype",
       title = "Coverage per chromosome arm")

ggsave(file.path(IMG_PATH, "arm_coverage_heatmap.pdf"), p_arm_coverage, width = 16, height = 4)

# ---- 7.3 Fisher tests per (karyotype, arm) ---------------------------------
# Tests: is compensation / dosage sensitivity enriched on this arm
# vs the rest of the genome for this karyotype?

fisher_one <- function(a, b, a_other, b_other) {
  ft <- suppressWarnings(
    fisher.test(matrix(c(a, b, a_other, b_other), nrow = 2))
  )
  dplyr::tibble(p = ft$p.value, OR = unname(ft$estimate))
}

fisher_per_arm <- function(count_col, total_col = "n_genes") {
  arm_coverage %>%
    dplyr::filter(n_genes >= 10) %>%
    dplyr::rename(focal   = !!count_col,
                  n_total = !!total_col) %>%
    dplyr::mutate(b = n_total - focal) %>%
    dplyr::group_by(karyotype) %>%
    dplyr::mutate(
      a_other = sum(focal) - focal,
      b_other = sum(b)     - b
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(fisher_one(focal, b, a_other, b_other)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(karyotype) %>%
    dplyr::mutate(fisher_padj = p.adjust(p, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(metric = count_col)
}

arm_fisher <- dplyr::bind_rows(
  fisher_per_arm("n_rna_comp"),
  fisher_per_arm("n_prot_comp"),
  fisher_per_arm("n_rna_ds"),
  fisher_per_arm("n_prot_ds")
) %>%
  dplyr::mutate(
    omic = dplyr::case_when(
      stringr::str_detect(metric, "rna")  ~ "RNA",
      stringr::str_detect(metric, "prot") ~ "Protein"
    ),
    type = dplyr::case_when(
      stringr::str_detect(metric, "comp") ~ "Compensated",
      stringr::str_detect(metric, "ds")   ~ "Dosage sensitive"
    ),
    omic = factor(omic, levels = c("RNA", "Protein")),
    type = factor(type, levels = c("Compensated", "Dosage sensitive")),
    arm  = factor(arm,  levels = arm_levels),
    karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping)
  )

saveRDS(arm_fisher, file.path(RES_PATH, "arm_fisher.rds"))

# ---- 7.4 Long-form fractions for raw heatmap ------------------------------
frac_by_arm <- arm_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::select(arm, karyotype, karyotype_lab,
                frac_RNA_comp, frac_Prot_comp,
                frac_RNA_ds,   frac_Prot_ds) %>%
  tidyr::pivot_longer(
    cols      = c(frac_RNA_comp, frac_Prot_comp, frac_RNA_ds, frac_Prot_ds),
    names_to  = "metric",
    values_to = "frac"
  ) %>%
  dplyr::mutate(
    omic = dplyr::case_when(
      stringr::str_detect(metric, "RNA")  ~ "RNA",
      stringr::str_detect(metric, "Prot") ~ "Protein"
    ),
    type = dplyr::case_when(
      stringr::str_detect(metric, "comp") ~ "Compensated",
      stringr::str_detect(metric, "ds")   ~ "Dosage sensitive"
    ),
    omic = factor(omic, levels = c("RNA", "Protein")),
    type = factor(type, levels = c("Compensated", "Dosage sensitive"))
  )

# ---- 7.5 Heatmap: raw fractions (2x2 facet: omic x type) ------------------
p_arm_heatmap_raw <- ggplot(frac_by_arm,
                            aes(x = arm, y = karyotype_lab, fill = frac)) +
  geom_tile(colour = "white") +
  facet_grid(type ~ omic) +
  scale_fill_gradient2(
    low      = "white",
    high     = "darkblue",
    midpoint = 0.3,
    limits   = c(0, 1),
    name     = "Fraction"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        strip.text  = element_text(face = "bold")) +
  labs(x = "Chromosome arm", y = "Karyotype",
       title = "Raw fraction compensated / dosage sensitive per arm")

ggsave(file.path(IMG_PATH, "arm_fraction_heatmap.pdf"),
       p_arm_heatmap_raw, width = 18, height = 8)

# ---- 7.6 Heatmap: Fisher OR + significance (main figure) ------------------
p_arm_heatmap_fisher <- arm_fisher %>%
  dplyr::mutate(
    log2_OR = log2(pmin(pmax(OR, 0.1), 10)),
    star    = sig_label(fisher_padj)
  ) %>%
  ggplot(aes(x = arm, y = karyotype_lab, fill = log2_OR)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 2.5, vjust = 0.75) +
  facet_grid(type ~ omic) +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0,
    limits   = c(-log2(10), log2(10)),
    name     = "log2(OR)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    strip.text  = element_text(face = "bold")
  ) +
  labs(
    x       = "Chromosome arm",
    y       = "Karyotype",
    title   = "Per-arm enrichment: compensated vs dosage sensitive",
    caption = "Fisher exact test vs rest of genome. BH-adjusted. * p<0.05, ** p<0.01, *** p<0.001."
  )

ggsave(file.path(IMG_PATH, "arm_fisher_heatmap.pdf"),
       p_arm_heatmap_fisher, width = 18, height = 10)
saveRDS(p_arm_heatmap_fisher,
        file.path(IMG_PATH, "arm_fisher_heatmap.rds"))

# ---- 7.7 Arm-level summary: top enriched/depleted arms per karyotype ------
top_arms <- arm_fisher %>%
  dplyr::filter(fisher_padj < 0.05) %>%
  dplyr::mutate(log2_OR = log2(pmin(pmax(OR, 0.1), 10))) %>%
  dplyr::arrange(karyotype, type, omic, dplyr::desc(abs(log2_OR))) %>%
  dplyr::group_by(karyotype, type, omic) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

saveRDS(top_arms, file.path(RES_PATH, "top_arms_fisher.rds"))
write.csv(top_arms %>%
            dplyr::select(karyotype, arm, omic, type,
                          OR, p, fisher_padj, log2_OR),
          file.path(RES_PATH, "top_arms_fisher.csv"),
          row.names = FALSE)

cat("\nTop enriched/depleted arms (padj < 0.05):\n")
print(top_arms %>%
        dplyr::select(karyotype_lab, arm, omic, type,
                      OR, fisher_padj) %>%
        dplyr::arrange(karyotype_lab, type, omic))

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
  dplyr::select(cn_direction, comp_RNA, comp_Protein) %>%
  tidyr::pivot_longer(c(comp_RNA, comp_Protein),
                      names_to = "level", values_to = "is_comp") %>%
  dplyr::mutate(level = sub("comp_", "", level)) %>%
  dplyr::group_by(cn_direction, level) %>%
  dplyr::summarise(stats = list(boot_fraction(is_comp)), .groups = "drop") %>%
  tidyr::unnest_wider(stats) %>%
  dplyr::mutate(level = factor(level, levels = c("RNA", "Protein")))

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
  corum_genes <- read.delim(corum_path) %>%
    dplyr::pull(GeneSym)
  
  df_classified <- df_classified %>%
    dplyr::mutate(in_complex = name %in% corum_genes)
  
  # ── Compute gap ───────────────────────────────────────────────────────────────
  df_gap <- df_classified %>%
    dplyr::mutate(
      gap           = RNA_lfc - Protein_lfc,
      karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping),
      complex_lab   = ifelse(in_complex, "CORUM complex", "Non-complex")
    )
  
  # ── Test 1: CORUM vs non-CORUM within each karyotype ─────────────────────────
  # Welch t-test on gap ~ in_complex per karyotype
  tests_complex <- df_gap %>%
    dplyr::group_by(karyotype, karyotype_lab) %>%
    dplyr::group_modify(~ {
      tt <- t.test(gap ~ in_complex, data = .x, var.equal = FALSE)
      dplyr::tibble(
        p          = tt$p.value,
        mean_non   = tt$estimate[1],   # in_complex = FALSE
        mean_corum = tt$estimate[2],   # in_complex = TRUE
        delta      = tt$estimate[2] - tt$estimate[1],
        ci_lo      = tt$conf.int[1],
        ci_hi      = tt$conf.int[2]
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(padj = p.adjust(p, method = "BH"))
  
  # ── Test 2: pairwise karyotype comparisons, stratified by complex membership ──
  # Are gains different from losses in how much gap they show?
  tests_karyotype <- df_gap %>%
    dplyr::group_by(in_complex, complex_lab) %>%
    dplyr::group_modify(~ {
      # Kruskal-Wallis first: is there any difference across karyotypes?
      kw <- kruskal.test(gap ~ karyotype, data = .x)
      # Pairwise Wilcoxon for post-hoc (robust to heavy tails)
      pw <- pairwise.wilcox.test(
        .x$gap, .x$karyotype,
        p.adjust.method = "BH"
      )
      dplyr::tibble(
        kw_p      = kw$p.value,
        pw_result = list(pw$p.value)   # matrix of pairwise p-values
      )
    }) %>%
    dplyr::ungroup()
  
  # ── Summary stats for plotting ────────────────────────────────────────────────
  gap_summary <- df_gap %>%
    dplyr::group_by(karyotype_lab, complex_lab) %>%
    dplyr::summarise(
      mean_gap = mean(gap,  na.rm = TRUE),
      sd_gap   = sd(gap,    na.rm = TRUE),
      n        = dplyr::n(),
      se_gap   = sd_gap / sqrt(n),
      ci_lo    = mean_gap - qt(0.975, df = n - 1) * se_gap,
      ci_hi    = mean_gap + qt(0.975, df = n - 1) * se_gap,
      .groups  = "drop"
    )
  
  saveRDS(df_gap,          file.path(RES_PATH, "df_gap.rds"))
  saveRDS(tests_complex,   file.path(RES_PATH, "gap_tests_complex.rds"))
  saveRDS(tests_karyotype, file.path(RES_PATH, "gap_tests_karyotype.rds"))
  
  annot_df <- tests_complex %>%
    dplyr::left_join(
      gap_summary %>%
        dplyr::group_by(karyotype_lab) %>%
        dplyr::summarise(y_top = max(ci_hi, na.rm = TRUE), .groups = "drop"),
      by = "karyotype_lab"
    ) %>%
    dplyr::mutate(
      y_position = y_top + 0.05,
      annotation = sig_label(padj)
    )
  
  library(ggpubr)
  
  # ── Plot 1: mean gap ± CI with ggpubr brackets ───────────────────────────────
  p_gap <- ggplot(gap_summary,
                  aes(x = karyotype_lab, y = mean_gap,
                      colour = complex_lab, group = complex_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    # p-value labels from tests_complex, placed just above the top CI
    geom_text(
      data = tests_complex %>%
        dplyr::left_join(
          gap_summary %>%
            dplyr::group_by(karyotype_lab) %>%
            dplyr::summarise(y_top = max(ci_hi, na.rm = TRUE), .groups = "drop"),
          by = "karyotype_lab"
        ) %>%
        dplyr::mutate(
          y_position = y_top + 0.03,
          label      = dplyr::case_when(
            padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ "ns"
          )
        ),
      aes(x = karyotype_lab, y = y_position, label = label),
      inherit.aes = FALSE,
      colour      = "black",
      size        = 4
    ) +
    scale_colour_manual(values = c("CORUM complex" = "steelblue3",
                                   "Non-complex"   = "grey50")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    labs(x       = "Karyotype",
         y       = "Mean RNA LFC − Protein LFC",
         colour  = "")
  
  # ── Plot 2: violin + boxplot with ggpubr stat_compare_means ──────────────────
  # Comparisons for the CORUM vs non-CORUM bracket within each karyotype facet
  complex_comparison <- list(c("CORUM complex", "Non-complex"))
  
  p_gap_violin <- ggplot(df_gap,
                         aes(x = complex_lab, y = gap, fill = complex_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_violin(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.3,
                 fill = "white") +
    stat_compare_means(
      comparisons  = complex_comparison,
      method       = "t.test",        # Welch t-test, consistent with tests_complex
      label        = "p.signif",      # **/*** labels; use "p.format" for numeric
      tip.length   = 0.01,
      bracket.size = 0.4
    ) +
    stat_compare_means(
      method   = "kruskal.test",      # overall across all x groups — not meaningful
      label.y  = -Inf,                # here x only has 2 levels so skip, see note
      hide.ns  = TRUE
    ) +
    scale_fill_manual(values = c("CORUM complex" = "steelblue3",
                                 "Non-complex"   = "grey70")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    facet_wrap(~karyotype_lab, nrow = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x       = "Complex membership",
         y       = "RNA LFC − Protein LFC",
         caption = "Welch t-test per karyotype, * p<0.05, ** p<0.01, *** p<0.001.")
  
  ggsave(file.path(IMG_PATH, "gap_complex_means.pdf"),  p_gap,           width = 8, height = 4)
  ggsave(file.path(IMG_PATH, "gap_complex_violin.pdf"), p_gap_violin,    width = 10, height = 4)
}

# ===========================================================================
# 11. GO enrichment per compensation class
# ===========================================================================

enrich_by_direction <- purrr::map(
  c("1:0", "2:1", "2:2"),
  function(dir) {
    
    genes <- df_classified %>%
      dplyr::filter(
        karyotype   == dir,
        class_label %in% c("RNA compensated", "Protein compensated",
                           "Fully compensated", "Full dosage sensitive")
      ) %>%
      dplyr::distinct(name, class_label) %>%
      split(.$class_label) %>%
      lapply(function(x) unique(x$name))
    
    if (length(genes) < 2) return(NULL)
    
    compareCluster(
      geneClusters = genes,
      fun           = "enrichGO",
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      pool          = FALSE
    )
  }
) %>%
  purrr::set_names(c("1:0", "2:1", "2:2")) %>%
  purrr::compact()

plots_direction <- purrr::imap(enrich_by_direction, ~ {
  dotplot(.x, showCategory = 8) +
    ggtitle(.y) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
})

p_direction <- patchwork::wrap_plots(plots_direction, ncol = 3) +
  patchwork::plot_annotation(title = "GO enrichment: Gains vs Losses")

ggsave(file.path(IMG_PATH, "GO_by_class_gain_loss.pdf"),
       p_direction, width = 18, height = 9)


# ---- 12.2 MitoCarta -------------------------------------------------------
mitocarta_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/Human.MitoCarta3.0.xls"

if (file.exists(mitocarta_path)) {
  mitocarta <- read_excel(mitocarta_path, sheet = "A Human MitoCarta3.0") %>%
    dplyr::as_tibble() %>%
    dplyr::pull(Symbol) %>%
    unique()
  
  # ── Gap-based MitoCarta analysis ───────────────────────────────────────────
  # Instead of the binary comp_Buffer flag, we test whether mitochondrial genes
  # show a larger RNA - Protein gap, using a Welch t-test on the continuous gap.
  # This is consistent with the rest of the gap-based analysis.
  
  df_mito <- df_classified %>%
    dplyr::mutate(
      gap    = RNA_lfc - Protein_lfc,
      in_mito = name %in% mitocarta
    )
  
  # ── Test: is the gap larger in MitoCarta genes per karyotype? ────────────────
  run_mito_ttest <- function(df, scope_label) {
    df %>%
      dplyr::group_by(karyotype) %>%
      dplyr::filter(dplyr::n() >= 10) %>%
      dplyr::group_modify(~ {
        tt <- t.test(gap ~ in_mito, data = .x, var.equal = FALSE)
        dplyr::tibble(
          mean_non_mito = tt$estimate[1],
          mean_mito     = tt$estimate[2],
          delta         = tt$estimate[2] - tt$estimate[1],  # mito - non-mito
          ci_lo         = tt$conf.int[1],
          ci_hi         = tt$conf.int[2],
          p             = tt$p.value
        )
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        padj  = p.adjust(p, method = "BH"),
        scope = scope_label
      )
  }
  
  mito_tests <- run_mito_ttest(df_mito, "Genome-wide")
  
  saveRDS(mito_tests, file.path(RES_PATH, "mitocarta_gap.rds"))
  
  # ── Summary for plotting: mean gap ± CI per karyotype × mito membership ──────
  mito_summary <- df_mito %>%
    dplyr::group_by(karyotype, in_mito) %>%
    dplyr::filter(dplyr::n() >= 10) %>%
    dplyr::summarise(
      mean_gap = mean(gap,  na.rm = TRUE),
      sd_gap   = sd(gap,    na.rm = TRUE),
      n        = dplyr::n(),
      se_gap   = sd_gap / sqrt(n),
      ci_lo    = mean_gap - qt(0.975, df = n - 1) * se_gap,
      ci_hi    = mean_gap + qt(0.975, df = n - 1) * se_gap,
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping),
      mito_lab      = ifelse(in_mito, "MitoCarta", "Non-mitochondrial")
    )
  
  # ── Significance annotations ──────────────────────────────────────────────────
  annot_mito <- mito_tests %>%
    dplyr::left_join(
      mito_summary %>%
        dplyr::group_by(karyotype) %>%
        dplyr::summarise(y_top = max(ci_hi, na.rm = TRUE), .groups = "drop"),
      by = "karyotype"
    ) %>%
    dplyr::mutate(
      karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping),
      y_position    = y_top + 0.03,
      label         = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ "ns"
      )
    )
  
  # ── Plot ──────────────────────────────────────────────────────────────────────
  p_mito <- ggplot(mito_summary,
                   aes(x = karyotype_lab, y = mean_gap,
                       colour = mito_lab, group = mito_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    geom_text(
      data        = annot_mito,
      aes(x = karyotype_lab, y = y_position, label = label),
      inherit.aes = FALSE,
      colour      = "black",
      size        = 4,
      fontface    = "bold"
    ) +
    scale_colour_manual(
      values = c("MitoCarta"          = "firebrick3",
                 "Non-mitochondrial"  = "grey50")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_bw() +
    labs(x       = "Karyotype",
         y       = "Mean RNA LFC − Protein LFC (gap)",
         colour  = "")
  
  # ── Violin for full distribution ──────────────────────────────────────────────
  p_mito_violin <- ggplot(
    df_mito %>%
      dplyr::mutate(
        karyotype_lab = factor(karyotype_mapping[karyotype], levels = karyotype_mapping),
        mito_lab      = ifelse(in_mito, "MitoCarta", "Non-mitochondrial")
      ),
    aes(x = mito_lab, y = gap, fill = mito_lab)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_violin(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(width = 0.2, outlier.size = 0.3, outlier.alpha = 0.3,
                 fill = "white") +
    stat_compare_means(
      method     = "t.test",
      label      = "p.signif",
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = c("MitoCarta"         = "firebrick3",
                 "Non-mitochondrial" = "grey70")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    facet_wrap(~karyotype_lab, nrow = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x       = "Gene set",
         y       = "RNA LFC − Protein LFC (gap)")
  
  ggsave(file.path(IMG_PATH, "mitocarta_gap_means.pdf"),
         p_mito,        width = 7,  height = 4)
  ggsave(file.path(IMG_PATH, "mitocarta_gap_violin.pdf"),
         p_mito_violin, width = 10, height = 4)
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
    n_full_dosage    = sum(class_label == "Full dosage"),
    frac_comp_RNA    = mean(comp_RNA),
    frac_comp_Protein = mean(comp_Protein),
    .groups = "drop"
  )

write_csv(gate_summary, file.path(RES_PATH, "gate_summary_per_karyotype.csv"))
saveRDS(gate_summary,    file.path(RES_PATH, "gate_summary_per_karyotype.rds"))

cat("Done.\n")
