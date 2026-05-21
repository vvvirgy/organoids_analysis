# ---- 1. Classify genes as compensated (CS > 0 and beyond noise envelope) ----

noise_model <- get_noise_model()

df_classified <- df %>%
  dplyr::select(!c(mu, sigma)) %>% 
  dplyr::left_join(
    noise_model$df %>% dplyr::select(omic, mu = mean, sigma),
    by = "omic"
  ) %>%
  dplyr::mutate(
    # upper-tail p-value: P(X >= CS) under the per-omic noise model
    p_cs = pnorm(CS, mean = mu, sd = sigma, lower.tail = FALSE)
  ) %>%
  dplyr::group_by(omic, karyotype) %>%
  dplyr::mutate(
    padj_cs  = p.adjust(p_cs, method = "BH"),
    CS_class = ifelse(CS > 0 & padj_cs <= 0.05, "CS>0", "CS<=0")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-mu, -sigma, -p_cs, -padj_cs)


# ---- 2. Bootstrap: compare fraction of compensated genes RNA vs Protein ----

compare_omics_boot <- function(sub_df, n_genes = 500, n_boot = 1000) {
  
  rna_labels  <- sub_df$CS_class[sub_df$omic == "RNA"]
  prot_labels <- sub_df$CS_class[sub_df$omic == "Protein"]
  
  if (length(rna_labels) == 0 || length(prot_labels) == 0) return(NULL)
  
  boot_results <- replicate(n_boot, {
    s_rna  <- sample(rna_labels,  size = n_genes, replace = TRUE)
    s_prot <- sample(prot_labels, size = n_genes, replace = TRUE)
    
    f_rna  <- mean(s_rna  == "CS>0")
    f_prot <- mean(s_prot == "CS>0")
    
    c(f_rna = f_rna, f_prot = f_prot, diff = f_prot - f_rna)
  }) %>% t() %>% as.data.frame()
  
  p_val <- mean(boot_results$diff <= 0)
  
  data.frame(
    omic   = c("RNA", "Protein"),
    f_mean = c(mean(boot_results$f_rna),            mean(boot_results$f_prot)),
    low    = c(quantile(boot_results$f_rna, 0.025), quantile(boot_results$f_prot, 0.025)),
    high   = c(quantile(boot_results$f_rna, 0.975), quantile(boot_results$f_prot, 0.975)),
    p_prot_gt_rna = p_val
  )
}


# ---- 3. Run the bootstrap per karyotype ----

plot_data <- df_classified %>%
  dplyr::group_by(karyotype) %>%
  tidyr::nest() %>%
  dplyr::mutate(stats = purrr::map(data, ~ compare_omics_boot(.x, n_genes = 500))) %>%
  tidyr::unnest(stats) %>%
  dplyr::select(-data) %>%
  dplyr::mutate(sig_50 = ifelse(low > 0.5, "*", "")) %>%
  dplyr::mutate(
    p_prot_gt_rna = ifelse(
      p_prot_gt_rna <= 0.001, "**",
      ifelse(p_prot_gt_rna <= 0.05, "*", "ns")
    )
  ) %>%
  dplyr::mutate(karyotype = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype = factor(karyotype, levels = karyotype_mapping))


# ---- 4. Bracket annotation dataset for geom_signif ----

pvals_df <- plot_data %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    y_position = max(high) + 0.075,
    p_label    = unique(p_prot_gt_rna),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    xmin       = as.numeric(karyotype) - 0.2,
    xmax       = as.numeric(karyotype) + 0.2,
    annotation = p_label
  )


# ---- 5. Plot ----

p_compensated_genes_fraction <- ggplot(plot_data, aes(x = karyotype, y = f_mean, fill = omic)) +
  
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # Star if lower CI > 50%
  geom_text(
    aes(label = sig_50, y = high + 0.02),
    position = position_dodge(width = 0.8),
    size = 6
  ) +
  
  # Brackets with p-values comparing Protein vs RNA
  geom_signif(
    data = pvals_df,
    aes(
      xmin        = xmin,
      xmax        = xmax,
      annotations = annotation,
      y_position  = y_position,
      group       = karyotype
    ),
    manual       = TRUE,
    textsize     = 3.5,
    tip_length   = 0.01,
    inherit.aes  = FALSE
  ) +
  
  theme_bw() +
  scale_fill_manual(values = omic_colors) +
  labs(
    y    = "Fraction of compensated genes",
    fill = "",
    x    = "Karyotype"
  )

p_compensated_genes_fraction
