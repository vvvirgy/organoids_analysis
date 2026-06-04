
rm(list = ls())
library(tidyverse)
library(MASS)
library(mvtnorm)

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

df_matched <- dplyr::inner_join(rna_df, prot_df, by = c("gene", "iteration")) %>%
  dplyr::mutate(lfc_buff = lfc_rna - lfc_prot)

cat("Matched draws used for fit:", nrow(df_matched), "\n")
cat("Unique genes in matched set:", dplyr::n_distinct(df_matched$gene), "\n")

# ── Fit marginal t distributions ───────────────────────────────────────────
# MASS::fitdistr estimates (m, s, df) via MLE — strictly better than
# moment matching for heavy-tailed data.
t_fit_rna  <- MASS::fitdistr(df_matched$lfc_rna,  "t")
t_fit_prot <- MASS::fitdistr(df_matched$lfc_prot, "t")

nu_rna  <- t_fit_rna$estimate["df"]
nu_prot <- t_fit_prot$estimate["df"]

cat("\nMarginal t fits:\n")
cat("  RNA    — mu:", t_fit_rna$estimate["m"],
    " sigma:", t_fit_rna$estimate["s"],
    " df:", nu_rna, "\n")
cat("  Protein — mu:", t_fit_prot$estimate["m"],
    " sigma:", t_fit_prot$estimate["s"],
    " df:", nu_prot, "\n")

# ── Fit bivariate-t on (RNA, Protein) ─────────────────────────────────────
# The bivariate-t requires a single shared df. We use the harmonic mean of
# the two marginal dfs — less sensitive to the smaller value than the
# arithmetic mean, less extreme than min().
nu_joint <- 2 / (1/nu_rna + 1/nu_prot)
cat("\nJoint df (harmonic mean of RNA and Protein):", nu_joint, "\n")

# Under a multivariate-t with df=ν, Cov(X) = Σ · ν/(ν−2).
# The sample covariance estimates Cov(X), NOT Σ directly → correct it.
xy      <- as.matrix(df_matched[, c("lfc_rna", "lfc_prot")])
mu_biv  <- c(t_fit_rna$estimate["m"], t_fit_prot$estimate["m"])
Sig_raw <- stats::cov(xy)
Sig_t   <- Sig_raw * (nu_joint - 2) / nu_joint   # scale matrix Σ (not covariance)
rho_t   <- Sig_t[1,2] / sqrt(Sig_t[1,1] * Sig_t[2,2])

cat("\nBivariate-t fit:\n")
cat("  mu_RNA  =", mu_biv[1], "  sigma_RNA  =", sqrt(Sig_t[1,1]), "\n")
cat("  mu_Prot =", mu_biv[2], "  sigma_Prot =", sqrt(Sig_t[2,2]), "\n")
cat("  rho     =", rho_t, "\n")
cat("  nu      =", nu_joint, "\n")

# ── Noise tibble (marginal t parameters) ──────────────────────────────────
diploid_noise <- dplyr::tibble(
  quantity = c("CS_RNA", "CS_Protein"),
  mu       = c(t_fit_rna$estimate["m"],
               t_fit_prot$estimate["m"]),
  sigma    = c(t_fit_rna$estimate["s"],
               t_fit_prot$estimate["s"]),
  nu       = c(nu_rna, nu_prot)
  # Thresholds: use qt(1 - ALPHA, df = nu) per row, NOT qnorm()
)

# ── Bivariate-t object ─────────────────────────────────────────────────────
# Conditional model: Prot | RNA = r  ~  t_ν(μ_{P|R}, σ²_{P|R})
# where μ_{P|R} = μ_P + β(r − μ_R),  β = Σ_PR / Σ_RR
# and   σ²_{P|R} = Σ_PP − Σ_PR²/Σ_RR  (scaled by (ν+1)/(ν + d_M) in strict form,
#       but the location/scale below is sufficient for threshold derivation)
beta_cond      <- Sig_t[1,2] / Sig_t[1,1]
intercept_cond <- mu_biv[2] - beta_cond * mu_biv[1]
sigma2_cond    <- Sig_t[2,2] - Sig_t[1,2]^2 / Sig_t[1,1]

diploid_noise_biv <- list(
  mu          = mu_biv,
  Sigma       = Sig_t,          # scale matrix (NOT covariance)
  Sigma_cov   = Sig_raw,        # sample covariance, kept for reference
  Sigma_inv   = solve(Sig_t),
  rho         = rho_t,
  nu          = nu_joint,
  cond_PgivenR = list(
    beta      = beta_cond,
    intercept = intercept_cond,
    sigma2    = sigma2_cond,
    nu        = nu_joint + 1    # conditional df increases by 1 for bivariate-t
  ),
  n           = nrow(df_matched)
)

print(diploid_noise)
saveRDS(diploid_noise,     OUT_PATH)
saveRDS(diploid_noise_biv, OUT_PATH_BIV)

# ── Diagnostic plots ───────────────────────────────────────────────────────

# (a) Marginal densities: empirical vs fitted t vs fitted Gaussian ──────────
plot_grid <- dplyr::bind_rows(
  dplyr::tibble(lfc = df_matched$lfc_rna,  quantity = "CS_RNA"),
  dplyr::tibble(lfc = df_matched$lfc_prot, quantity = "CS_Protein")
) %>%
  dplyr::mutate(quantity = factor(quantity,
                                  levels = c("CS_RNA", "CS_Protein")))

# Build fitted-density curves for both t and Gaussian per panel
fit_curves <- purrr::pmap_dfr(
  list(
    qty   = c("CS_RNA",    "CS_Protein"),
    mu_t  = c(t_fit_rna$estimate["m"],  t_fit_prot$estimate["m"]),
    s_t   = c(t_fit_rna$estimate["s"],  t_fit_prot$estimate["s"]),
    df_t  = c(nu_rna,      nu_prot),
    mu_g  = c(t_fit_rna$estimate["m"],  t_fit_prot$estimate["m"]),
    s_g   = c(sqrt(Sig_raw[1,1]),       sqrt(Sig_raw[2,2]))
  ),
  function(qty, mu_t, s_t, df_t, mu_g, s_g) {
    x <- seq(
      quantile(df_matched[[switch(qty,
                                  CS_RNA = "lfc_rna", CS_Protein = "lfc_prot", CS_Buffering = "lfc_buff"
      )]], 0.001),
      quantile(df_matched[[switch(qty,
                                  CS_RNA = "lfc_rna", CS_Protein = "lfc_prot", CS_Buffering = "lfc_buff"
      )]], 0.999),
      length.out = 500
    )
    dplyr::bind_rows(
      dplyr::tibble(x, y = dt((x - mu_t)/s_t, df = df_t)/s_t,
                    model = "t (fitted)",   quantity = qty),
      dplyr::tibble(x, y = dnorm(x, mu_g, s_g),
                    model = "Gaussian",     quantity = qty)
    )
  }
) %>%
  dplyr::mutate(quantity = factor(quantity,
                                  levels = c("CS_RNA", "CS_Protein", "CS_Buffering")))

p_marg <- ggplot(plot_grid, aes(x = lfc)) +
  geom_density(aes(fill = quantity), alpha = 0.4) +
  geom_line(data = fit_curves, aes(x = x, y = y, colour = model,
                                   linetype = model), linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("t (fitted)" = "firebrick3", "Gaussian" = "steelblue3")) +
  scale_linetype_manual(values = c("t (fitted)" = "solid",   "Gaussian" = "dashed")) +
  facet_wrap(~quantity, nrow = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "Diploid LFC (null)", y = "Density",
       colour = "Model", linetype = "Model",
       subtitle = "Empirical density (fill) vs fitted t (red) and Gaussian (blue)") +
  theme(legend.position = "bottom")

ggsave(file.path(IMG_PATH, "noise_model_null_densities_t.pdf"), p_marg, width = 10, height = 4)

# (b) Q-Q plots against fitted t (not Gaussian) ────────────────────────────
qq_df <- dplyr::bind_rows(
  dplyr::tibble(lfc = df_matched$lfc_rna,  quantity = "RNA"),
  dplyr::tibble(lfc = df_matched$lfc_prot, quantity = "Protein")
) %>%
  dplyr::mutate(quantity = factor(quantity, levels = c("RNA", "Protein"))) %>%
  dplyr::group_by(quantity) %>%
  dplyr::mutate(
    theoretical_gauss = qnorm(ppoints(dplyr::n()))[rank(lfc, ties.method = "first")],
    theoretical_t     = qt(ppoints(dplyr::n()), df = diploid_noise$nu[
      match(dplyr::cur_group()$quantity,
            sub("CS_", "", diploid_noise$quantity))
    ])[rank(lfc, ties.method = "first")]
  ) %>%
  dplyr::ungroup()

abline_df <- diploid_noise %>%
  dplyr::mutate(
    quantity = factor(sub("CS_", "", quantity), levels = c("RNA", "Protein"))
  )

# Gaussian Q-Q
p_qq_gauss <- ggplot(qq_df, aes(x = theoretical_gauss, y = lfc)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_abline(data = abline_df, aes(slope = sigma, intercept = mu),
              colour = "steelblue3", linewidth = 0.8) +
  facet_wrap(~quantity, scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical quantiles — Gaussian", y = "Sample LFC",
       subtitle = "Q-Q vs Gaussian fit (blue). S-shape = heavy tails")

# t Q-Q
p_qq_t <- ggplot(qq_df, aes(x = theoretical_t, y = lfc)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_abline(data = abline_df, aes(slope = sigma, intercept = mu),
              colour = "firebrick3", linewidth = 0.8) +
  facet_wrap(~quantity, scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical quantiles — t (fitted df)", y = "Sample LFC",
       subtitle = "Q-Q vs fitted t (red). Points should lie on line if t fits well")

ggsave(file.path(IMG_PATH, "noise_model_qq_gaussian.pdf"), p_qq_gauss, width = 9, height = 3.5)
ggsave(file.path(IMG_PATH, "noise_model_qq_t.pdf"),        p_qq_t,     width = 9, height = 3.5)

# (c) Joint scatter with bivariate-t contours ──────────────────────────────
# Contours are iso-density curves of the fitted bivariate-t,
# analogous to the Mahalanobis ellipses for the Gaussian case.
make_t_contour <- function(mu, Sigma, nu, prob, npts = 300) {
  # The p-th iso-density contour of a bivariate-t with df=nu is an ellipse
  # scaled by sqrt(nu * (p^(-2/nu) - 1) / (nu - 2)) — i.e. the quantile of
  # the F distribution maps back to a radius.
  # Simpler: use the fact that Mahalanobis² ~ F(2, nu)*2 for bivariate-t.
  r     <- sqrt(2 * qf(prob, df1 = 2, df2 = nu))
  theta <- seq(0, 2*pi, length.out = npts)
  L     <- chol(Sigma)
  pts   <- t(mu + r * t(L) %*% rbind(cos(theta), sin(theta)))
  dplyr::tibble(x = pts[,1], y = pts[,2],
                level = paste0(round(prob*100), "% density"))
}

contour_df <- dplyr::bind_rows(
  make_t_contour(mu_biv, Sig_t, nu_joint, 0.50),
  make_t_contour(mu_biv, Sig_t, nu_joint, 0.90),
  make_t_contour(mu_biv, Sig_t, nu_joint, 0.95)
)

plot_pts <- df_matched
if (nrow(plot_pts) > 20000) plot_pts <- dplyr::slice_sample(plot_pts, n = 20000)

p_joint <- ggplot(plot_pts, aes(lfc_rna, lfc_prot)) +
  geom_point(alpha = 0.08, size = 0.4) +
  geom_path(data = contour_df,
            aes(x, y, group = level, colour = level), linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  scale_colour_manual(values = c(
    "50% density" = "firebrick3",
    "90% density" = "steelblue3",
    "95% density" = "forestgreen"
  )) +
  coord_fixed() +
  theme_bw() +
  labs(x = "Diploid RNA LFC", y = "Diploid Protein LFC",
       colour = "Bivariate-t contour",
       subtitle = sprintf("Bivariate-t null: rho = %.3f, nu = %.1f, n = %d",
                          rho_t, nu_joint, nrow(df_matched)))

ggsave(file.path(IMG_PATH, "noise_model_joint_t_contours.pdf"), p_joint, width = 5.5, height = 5)
