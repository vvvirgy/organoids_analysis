rm(list = ls())
library(tidyverse)

# ── Paths ──────────────────────────────────────────────────────────────────
PROT_BOOT_PATH <- "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/noise_model/protein/dep_diploid_genes_results.rds"
RNA_BOOT_PATH  <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/diploid_bootstrap.rds"
OUT_PATH       <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise.RDS"
OUT_PATH_BIV   <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise_bivariate.RDS"
IMG_PATH       <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/img/"

# ── Load bootstrap data ────────────────────────────────────────────────────
rna_df <- readRDS(RNA_BOOT_PATH) %>%
  dplyr::filter(min_mean_expr > 0.1, non_zero_percent > 1, n_samples > 12) %>%
  dplyr::select(gene, iteration, lfc) %>%
  dplyr::rename(lfc_rna = lfc)

prot_df <- readRDS(PROT_BOOT_PATH) %>%
  dplyr::select(gene, iteration, lfc) %>%
  dplyr::rename(lfc_prot = lfc) %>%
  na.omit()

# Matched draws: same gene × iteration in both omics.
# This is the ONLY data source for the bivariate fit, so the marginals are
# now also estimated from this matched set (was: separate filters per omic).
df_matched <- dplyr::inner_join(rna_df, prot_df, by = c("gene", "iteration")) %>%
  dplyr::mutate(lfc_buff = lfc_rna - lfc_prot)

cat("Matched draws used for bivariate fit:", nrow(df_matched), "\n")
cat("Unique genes in matched set:        ", dplyr::n_distinct(df_matched$gene), "\n")

# ── Fit bivariate Gaussian on (RNA, Protein) ───────────────────────────────
# Σ is estimated from matched draws → ρ carries the actual coupling.
xy   <- as.matrix(df_matched[, c("lfc_rna", "lfc_prot")])
mu2  <- colMeans(xy)                     # length-2 vector: (μ_R, μ_P)
Sig  <- stats::cov(xy)                   # 2×2 covariance
rho  <- Sig[1, 2] / sqrt(Sig[1, 1] * Sig[2, 2])

cat("Bivariate fit:\n")
cat("  μ_RNA  =", mu2[1], "  σ_RNA  =", sqrt(Sig[1, 1]), "\n")
cat("  μ_Prot =", mu2[2], "  σ_Prot =", sqrt(Sig[2, 2]), "\n")
cat("  ρ      =", rho, "\n")

# ── Derive the three marginal nulls analytically from Σ ────────────────────
# RNA marginal:    N(μ_R, σ_R²)
# Protein marg.:   N(μ_P, σ_P²)
# Buffer = R − P:  N(μ_R − μ_P, σ_R² + σ_P² − 2ρσ_Rσ_P)
sd_R    <- sqrt(Sig[1, 1])
sd_P    <- sqrt(Sig[2, 2])
mu_buf  <- mu2[1] - mu2[2]
var_buf <- Sig[1, 1] + Sig[2, 2] - 2 * Sig[1, 2]
sd_buf  <- sqrt(var_buf)

# Sanity check: var_buf derived from Σ should equal the empirical var of lfc_buff
cat("Var(R−P) from Σ:        ", var_buf, "\n")
cat("Var(R−P) empirical:     ", var(df_matched$lfc_buff), "\n")

diploid_noise <- dplyr::tibble(
  quantity = c("CS_RNA", "CS_Protein", "CS_Buffering"),
  mu       = c(mu2[1],   mu2[2],       mu_buf),
  sigma    = c(sd_R,     sd_P,         sd_buf),
  n        = c(nrow(df_matched), nrow(df_matched), nrow(df_matched))
)

# Bivariate object for the joint Mahalanobis test downstream
diploid_noise_biv <- list(
  mu          = mu2,           # named: lfc_rna, lfc_prot
  Sigma       = Sig,
  Sigma_inv   = solve(Sig),
  rho         = rho,
  # Conditional model: Prot | RNA = r  ~  N(μ_P + β(r − μ_R), σ²_{P|R})
  # Used for the buffering-as-conditional-test interpretation.
  cond_PgivenR = list(
    beta      = Sig[1, 2] / Sig[1, 1],        # ρ · σ_P / σ_R
    intercept = mu2[2] - (Sig[1, 2] / Sig[1, 1]) * mu2[1],
    sigma2    = Sig[2, 2] - Sig[1, 2]^2 / Sig[1, 1]   # σ_P² (1 − ρ²)
  ),
  n           = nrow(df_matched)
)

print(diploid_noise)
print(diploid_noise_biv)

saveRDS(diploid_noise,     OUT_PATH)        # backward-compatible: same shape as before
saveRDS(diploid_noise_biv, OUT_PATH_BIV)    # NEW: full bivariate object

# ── Diagnostic plots ───────────────────────────────────────────────────────
# (a) Three marginal densities with ±1σ lines, as before
df_density <- dplyr::bind_rows(
  dplyr::tibble(lfc = df_matched$lfc_rna,  quantity = "CS_RNA"),
  dplyr::tibble(lfc = df_matched$lfc_prot, quantity = "CS_Protein"),
  dplyr::tibble(lfc = df_matched$lfc_buff, quantity = "CS_Buffering")
) %>%
  dplyr::mutate(quantity = factor(quantity,
                                  levels = c("CS_RNA", "CS_Protein", "CS_Buffering")))

sigma_lines <- diploid_noise %>%
  dplyr::mutate(quantity = factor(quantity,
                                  levels = c("CS_RNA", "CS_Protein", "CS_Buffering")))

p_marg <- ggplot(df_density, aes(x = lfc, fill = quantity)) +
  geom_density(alpha = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_vline(data = sigma_lines, aes(xintercept =  sigma),
             linetype = "dotted", colour = "firebrick3") +
  geom_vline(data = sigma_lines, aes(xintercept = -sigma),
             linetype = "dotted", colour = "firebrick3") +
  facet_wrap(~quantity, nrow = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "Diploid LFC (null)", y = "Density",
       subtitle = "Marginals derived from bivariate fit. Red dotted lines = ±1σ") +
  theme(legend.position = "none")

ggsave(file.path(IMG_PATH, "noise_model_null_densities_biv.pdf"),
       p_marg, width = 9, height = 3.5)
ggsave(file.path(IMG_PATH, "noise_model_null_densities_biv.png"),
       p_marg, width = 9, height = 3.5, dpi = 450)

# (b) Joint scatter with 1σ / 2σ Mahalanobis ellipses
# Ellipse: points x with (x−μ)ᵀ Σ⁻¹ (x−μ) = k². For k=1 → 1-σ, k=2 → 2-σ.
make_ellipse <- function(mu, Sigma, k, npts = 200) {
  theta <- seq(0, 2 * pi, length.out = npts)
  circle <- rbind(cos(theta), sin(theta))
  L <- chol(Sigma)                       # L'L = Σ
  pts <- t(mu + k * t(L) %*% circle)
  dplyr::tibble(x = pts[, 1], y = pts[, 2], level = paste0(k, "σ"))
}
ell_df <- dplyr::bind_rows(
  make_ellipse(mu2, Sig, 1),
  make_ellipse(mu2, Sig, 2)
)

# Downsample points for visualisation if very large
plot_pts <- df_matched
if (nrow(plot_pts) > 20000) plot_pts <- dplyr::slice_sample(plot_pts, n = 20000)

p_joint <- ggplot(plot_pts, aes(lfc_rna, lfc_prot)) +
  geom_point(alpha = 0.08, size = 0.4) +
  geom_path(data = ell_df, aes(x, y, group = level, colour = level),
            linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  scale_colour_manual(values = c(`1σ` = "firebrick3", `2σ` = "steelblue3")) +
  coord_fixed() +
  theme_bw() +
  labs(x = "Diploid RNA LFC", y = "Diploid Protein LFC",
       colour = "Mahalanobis",
       subtitle = sprintf("Bivariate null: ρ = %.3f, n = %d matched draws",
                          rho, nrow(df_matched)))

ggsave(file.path(IMG_PATH, "noise_model_joint_ellipses.pdf"),
       p_joint, width = 5.5, height = 5)
ggsave(file.path(IMG_PATH, "noise_model_joint_ellipses.png"),
       p_joint, width = 5.5, height = 5, dpi = 450)

library(moments)   # skewness, kurtosis, jarque.test
library(MASS)      # for mvn alternative if needed

# ── Per-marginal normality diagnostics ────────────────────────────────────────
marginal_diag <- dplyr::tibble(
  quantity = c("RNA", "Protein", "Buffer"),
  data     = list(df_matched$lfc_rna, df_matched$lfc_prot, df_matched$lfc_buff)
) %>%
  dplyr::mutate(
    skewness  = purrr::map_dbl(data, moments::skewness),
    ex_kurt   = purrr::map_dbl(data, moments::kurtosis) - 3,  # excess kurtosis
    sw_p      = purrr::map_dbl(data, ~ {
      # Shapiro-Wilk is unreliable for n > 5000; subsample
      x <- .x
      if (length(x) > 5000) x <- sample(x, 5000)
      shapiro.test(x)$p.value
    }),
    jb_p      = purrr::map_dbl(data, ~ moments::jarque.test(.x)$p.value)
  ) %>%
  dplyr::select(-data)

print(marginal_diag)
# Rule of thumb: |skewness| > 0.5 is meaningful; excess kurtosis > 1 → heavy tails

# ── Q-Q plots for the three marginals ─────────────────────────────────────────
qq_df <- dplyr::bind_rows(
  dplyr::tibble(lfc = df_matched$lfc_rna,  quantity = "RNA"),
  dplyr::tibble(lfc = df_matched$lfc_prot, quantity = "Protein"),
  dplyr::tibble(lfc = df_matched$lfc_buff, quantity = "Buffer")
) %>%
  dplyr::mutate(quantity = factor(quantity, levels = c("RNA", "Protein", "Buffer"))) %>%
  dplyr::group_by(quantity) %>%
  dplyr::mutate(
    theoretical = qnorm(ppoints(dplyr::n()))[rank(lfc)]
  ) %>%
  dplyr::ungroup()

p_qq <- ggplot(qq_df, aes(x = theoretical, y = lfc)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_abline(slope  = 1, intercept = 0,           # wrong: should use per-group mu/sigma
              linetype = "dashed", colour = "grey40") +
  # Better: draw the fitted Gaussian line using per-group mu/sigma
  geom_abline(
    data = diploid_noise %>%
      dplyr::mutate(quantity = sub("CS_", "", quantity),
                    quantity = factor(quantity, levels = c("RNA", "Protein", "Buffer"))),
    aes(slope = sigma, intercept = mu),
    colour = "firebrick3", linewidth = 0.7
  ) +
  facet_wrap(~quantity, scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical quantiles (N(0,1))", y = "Sample LFC",
       subtitle = "Red line = fitted Gaussian. Deviations in tails → non-normality")

ggsave(file.path(IMG_PATH, "noise_model_qqplots.pdf"), p_qq, width = 9, height = 3.5)




library(MASS)

t_fit_rna  <- MASS::fitdistr(df_matched$lfc_rna,  "t")
t_fit_prot <- MASS::fitdistr(df_matched$lfc_prot, "t")

cat("RNA  — df:", t_fit_rna$estimate["df"],
    " mu:", t_fit_rna$estimate["m"],
    " sigma:", t_fit_rna$estimate["s"], "\n")
cat("Prot — df:", t_fit_prot$estimate["df"],
    " mu:", t_fit_prot$estimate["m"],
    " sigma:", t_fit_prot$estimate["s"], "\n")



# Compare thresholds at your typical alpha levels
alphas <- c(0.05, 0.025, 0.01, 0.005)

threshold_comparison <- dplyr::tibble(
  alpha         = alphas,
  gaussian      = qnorm(1 - alphas),
  t_RNA         = qt(1 - alphas, df = t_fit_rna$estimate["df"]),
  t_Protein     = qt(1 - alphas, df = t_fit_prot$estimate["df"])
) %>%
  dplyr::mutate(
    inflation_RNA  = t_RNA  / gaussian - 1,
    inflation_Prot = t_Protein / gaussian - 1
  )

print(threshold_comparison)
