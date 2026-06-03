
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(boot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))

# source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/utils.R")
# source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/constants.R")
# source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/getters.R")
# IMG_PATH <- paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/img/sf_", sf_method, "_stable_", use_stable)
# RES_PATH <- paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/sf_", sf_method, "_stable_", use_stable)
# DF_DNA_PATH = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/DNA_lfc.rds"
# DF_CS_SCORES_PATH = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/CS_scores_prot_and_rna.rds"
# NOISE_MODEL_PATH = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/multiOmic/diploid_noise.RDS"

source("utils.R")
source("scripts/constants.R")
source("scripts/getters.R")
IMG_PATH <- paste0("img/sf_", sf_method, "_stable_", use_stable)
RES_PATH <- paste0("results/sf_", sf_method, "_stable_", use_stable)
DF_DNA_PATH = "results/DNA_lfc.rds"
DF_CS_SCORES_PATH = "results/multiOmic/CS_scores_prot_and_rna.rds"
NOISE_MODEL_PATH = "results/multiOmic/diploid_noise.RDS"

dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.01
MIN_N_PER_GROUP <- 4   # per-karyotype minimum sample size filter

# ---------------------------------------------------------------------------
# 1. Load data and noise model
# ---------------------------------------------------------------------------

df_dna <- readRDS(DF_DNA_PATH) %>%
  dplyr::rename(DNA_lfc = lfc)

df_raw <- readRDS(DF_CS_SCORES_PATH)

# Noise model now contains three quantities: RNA, Protein, Buffering

noise_model <- get_noise_model() %>%
  dplyr::mutate(omic = str_replace(quantity, "CS_", "")) %>%
  dplyr::select(mu, sigma, omic)

#noise_model$sigma = 0.3
# Expecting columns: mu, sigma, omic with values in c("RNA","Protein","Buffering")

# ---------------------------------------------------------------------------
# 2. Build the unified per-(gene, karyotype) table with all three CS scores
#
#    Conventions (positive = compensation):
#      CS_RNA      = DNA_lfc - RNA_lfc
#      CS_Protein  = DNA_lfc - Protein_lfc
#      CS_Buffer   = RNA_lfc - Protein_lfc
# ---------------------------------------------------------------------------
df_wide <- df_raw %>%
  dplyr::select(name, karyotype, omic, lfc, pval) %>%
  tidyr::pivot_wider(names_from = omic, values_from = c(lfc, pval)) %>%
  dplyr::rename(RNA_lfc = lfc_RNA, Protein_lfc = lfc_Protein,
                RNA_pval = pval_RNA, Protein_pval = pval_Protein) %>%
  dplyr::left_join(df_dna %>% dplyr::select(name, karyotype, DNA_lfc),
                   by = c("name", "karyotype")) %>%
  dplyr::filter(!is.na(RNA_lfc), !is.na(Protein_lfc), !is.na(DNA_lfc)) %>%
  dplyr::mutate(
    CS_RNA     = (DNA_lfc - RNA_lfc)     * sign(DNA_lfc),
    CS_Protein = (DNA_lfc - Protein_lfc) * sign(DNA_lfc),
    CS_Buffer  = (RNA_lfc - Protein_lfc) * sign(DNA_lfc)   # sign-correct this too
  )
saveRDS(df_wide, file.path(RES_PATH, "df_wide_three_CS.rds"))

# ---------------------------------------------------------------------------
# 3. Per-gene classification: compensated at each level?
#    Upper-tail one-sided test against the appropriate noise model.
# ---------------------------------------------------------------------------
nm <- noise_model %>% dplyr::select(omic, mu, sigma)
nm_rna    <- nm %>% dplyr::filter(omic == "RNA")
nm_prot   <- nm %>% dplyr::filter(omic == "Protein")
nm_buffer <- nm %>% dplyr::filter(omic == "Buffering")


z_threshold <- qnorm(ALPHA, lower.tail = FALSE)  # = 1.645
df_classified <- df_wide %>%
  # dplyr::mutate(
  #   p_RNA     = pnorm(CS_RNA,     mean = nm_rna$mu,    sd = nm_rna$sigma,    lower.tail = FALSE),
  #   p_Protein = pnorm(CS_Protein, mean = nm_prot$mu,   sd = nm_prot$sigma,   lower.tail = FALSE),
  #   p_Buffer  = pnorm(CS_Buffer,  mean = nm_buffer$mu, sd = nm_buffer$sigma, lower.tail = FALSE)
  # ) %>%
  # #dplyr::group_by(karyotype) %>%
  # dplyr::mutate(
  #   padj_RNA     = p.adjust(p_RNA,     method = "BH"),
  #   padj_Protein = p.adjust(p_Protein, method = "BH"),
  #   padj_Buffer  = p.adjust(p_Buffer,  method = "BH")
  # ) %>%
  # dplyr::ungroup() %>%
  # dplyr::mutate(
  #   comp_RNA     = (CS_RNA     > 0) & (padj_RNA     <= ALPHA),
  #   comp_Protein = (CS_Protein > 0) & (padj_Protein <= ALPHA),
  #   comp_Buffer  = (CS_Buffer  > 0) & (padj_Buffer  <= ALPHA)
  # ) %>%
  dplyr::mutate(
    z_RNA     = CS_RNA     / nm_rna$sigma,
    z_Protein = CS_Protein / nm_prot$sigma,
    z_Buffer  = CS_Buffer  / nm_buffer$sigma,

    comp_RNA     = z_RNA     >= z_threshold,
    comp_Protein = z_Protein >= z_threshold,
    comp_Buffer  = z_Buffer  >= z_threshold
  ) %>%
  dplyr::mutate(
    class_3bit = paste0(
      as.integer(comp_RNA), as.integer(comp_Protein), as.integer(comp_Buffer)
    ),
    class_label = dplyr::case_when(
      class_3bit == "000" ~ "Full dosage",
      class_3bit == "110" ~ "RNA-compensated",
      class_3bit == "011" ~ "Buffered (post-transcriptional)",
      class_3bit == "111" ~ "Compensated at both levels",
      class_3bit == "010" ~ "Protein-only compensated",
      class_3bit == "001" ~ "Buffering only",
      class_3bit == "100" ~ "RNA-only (unstable)",
      class_3bit == "101" ~ "RNA-comp + amplified protein",
      TRUE                 ~ "Other"
    )
  )

# ---------------------------------------------------------------------------
# 3b. Up / Down / Not-differential classification per omic
#
#     A gene is called Up (Down) if:
#       - its LFC is significantly different from zero under the noise model
#         (two-sided test against N(mu_omic, sigma_omic))
#       - AND its DE p-value (BH-adjusted within karyotype x omic) <= ALPHA
#       - AND the sign of the LFC matches the direction
#
#     Two-sided here because we want to detect deviations from zero in either
#     direction; the sign filter handles direction assignment.
# ---------------------------------------------------------------------------

classify_direction <- function(lfc, pval_adj, mu, sigma, alpha) {
  z          <- (lfc - mu) / sigma
  p_noise    <- 2 * pnorm(-abs(z))                          # two-sided noise-model p
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
    dir_RNA = classify_direction(RNA_lfc,     RNA_padj, nm_rna$mu,   nm_rna$sigma,  ALPHA),
    dir_Protein = classify_direction(Protein_lfc, Protein_padj, nm_prot$mu,  nm_prot$sigma, ALPHA)
  )

saveRDS(df_classified, file.path(RES_PATH, "df_classified_three_CS.rds"))

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

# Asymmetry score: (Up - Down) / (Up + Down). Excludes Not-differential.
dir_asymmetry <- dir_counts %>%
  dplyr::filter(direction %in% c("Up", "Down")) %>%
  dplyr::select(karyotype, omic, direction, n) %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  dplyr::mutate(
    n_diff    = Up + Down,
    asymmetry = ifelse(n_diff == 0, NA_real_, (Up - Down) / n_diff),
    # binomial test: is Up vs Down deviation from 50/50 significant?
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
# 3d. Plots
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

# ---------------------------------------------------------------------------
# 4. Sanity-check plot: distribution of LFCs vs DNA expectation
# ---------------------------------------------------------------------------

df_long_lfc <- df_classified %>%
  dplyr::select(name, karyotype, DNA_lfc, RNA_lfc, Protein_lfc) %>%
  tidyr::pivot_longer(c(RNA_lfc, Protein_lfc), names_to = "omic", values_to = "lfc") %>%
  dplyr::mutate(omic = sub("_lfc", "", omic)) %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

dna_lines <- df_classified %>%
  dplyr::distinct(karyotype, DNA_lfc) %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

p_lfc_density <- ggplot(df_long_lfc, aes(x = lfc, fill = omic, colour = omic)) +
  geom_density(alpha = 0.3) +
  #geom_vline(data = dna_lines, aes(xintercept = DNA_lfc), linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~karyotype_lab, scales = "free_y") +
  scale_fill_manual(values = omic_colors) +
  scale_colour_manual(values = omic_colors) +
  theme_bw() +
  labs(x = "log2 fold-change vs diploid", y = "Density")

ggsave(file.path(IMG_PATH, "lfc_density_vs_dosage.pdf"), p_lfc_density, width = 9, height = 6)
saveRDS(p_lfc_density, file.path(IMG_PATH, "lfc_density_vs_dosage.rds"))

# ---------------------------------------------------------------------------
# 5. Genome-wide fraction compensated per karyotype, with bootstrap CIs
#    Three levels: RNA, Protein, Buffering
# ---------------------------------------------------------------------------

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
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)) %>%
  dplyr::mutate(level = factor(level, levels = c("RNA", "Protein", "Buffer")))

p_frac_compensated <- ggplot(frac_by_karyotype,
                             aes(x = karyotype_lab, y = mean, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = low, ymax = high),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of compensated genes",
       fill = "Compensation level")

ggsave(file.path(IMG_PATH, "fraction_compensated_three_levels.pdf"), p_frac_compensated, width = 7, height = 4)
saveRDS(p_frac_compensated, file.path(IMG_PATH, "fraction_compensated_three_levels.rds"))

# ---------------------------------------------------------------------------
# 6. RNA-LFC vs Protein-LFC scatter (buffering visualization)
# ---------------------------------------------------------------------------

dna_pts <- df_classified %>%
  dplyr::distinct(karyotype, DNA_lfc) %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

p_rna_vs_prot <- df_classified %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)) %>%
  ggplot(aes(x = RNA_lfc, y = Protein_lfc)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = comp_Buffer), alpha = 0.4, size = 0.6) +
  facet_wrap(~karyotype_lab) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "darkorange")) +
  theme_bw() +
  labs(x = "RNA log2FC", y = "Protein log2FC", colour = "Buffered")

ggsave(file.path(IMG_PATH, "rna_vs_protein_lfc.pdf"), p_rna_vs_prot, width = 8, height = 6)
saveRDS(p_rna_vs_prot, file.path(IMG_PATH, "rna_vs_protein_lfc.rds"))

# ---------------------------------------------------------------------------
# 7. Per-chromosome compensation analysis
#
#    Three views:
#      (a) Raw fraction of compensated genes per (karyotype, chromosome)
#      (b) Enrichment relative to karyotype-wide baseline (log2)
#      (c) Fisher test per (karyotype, chromosome) vs rest of genome
#    Plus a chromosome-specific CORUM check for the buffered set.
# ---------------------------------------------------------------------------

# ---- 7.0 Gene -> chromosome lookup via biomaRt -----------------------------
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

# Report mapping coverage so we know how many genes are dropped
n_total   <- length(unique(df_classified$name))
n_mapped  <- nrow(gene_chrom)
message(sprintf("Chromosome mapping: %d / %d genes mapped (%.1f%%)",
                n_mapped, n_total, 100 * n_mapped / n_total))

if (nrow(gene_chrom) == 0) stop("No genes mapped to chromosomes; aborting section 7.")

chrom_levels <- c(as.character(1:22), "X", "Y")

# ---- 7.1 Per-chromosome coverage (sample sizes) ----------------------------
chrom_coverage <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(!is.na(chr)) %>%
  dplyr::group_by(chr, karyotype) %>%
  dplyr::summarise(
    n_genes  = dplyr::n(),
    n_buffer = sum(comp_Buffer),
    n_rna    = sum(comp_RNA),
    n_prot   = sum(comp_Protein),
    frac_Buffer  = mean(comp_Buffer),
    frac_RNA     = mean(comp_RNA),
    frac_Protein = mean(comp_Protein),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_levels))

saveRDS(chrom_coverage, file.path(RES_PATH, "chrom_coverage.rds"))

# Diagnostic plot: how many genes per (karyotype, chromosome)?
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

# ---- 7.2 Long-form fractions for heatmap -----------------------------------
frac_by_chrom <- chrom_coverage %>%
  dplyr::filter(n_genes >= 10) %>%
  dplyr::select(karyotype, chr, frac_RNA, frac_Protein, frac_Buffer) %>%
  tidyr::pivot_longer(c(frac_RNA, frac_Protein, frac_Buffer),
                      names_to = "level", values_to = "frac") %>%
  dplyr::mutate(level = sub("frac_", "", level),
                level = factor(level, levels = c("RNA", "Protein", "Buffer")),
                karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

# ---- 7.3 Enrichment vs karyotype-wide baseline -----------------------------
# Tells us whether a chromosome is buffering more than its karyotype's average,
# not just inheriting the genome-wide rate.
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
  dplyr::mutate(level = sub("enr_", "", level),
                level = factor(level, levels = c("RNA", "Protein", "Buffer")),
                karyotype_lab = karyotype_mapping[karyotype],
                karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

saveRDS(enrichment_by_chrom, file.path(RES_PATH, "chrom_enrichment.rds"))

# ---- 7.4 Fisher tests per (karyotype, chromosome) --------------------------
# Tests: among compensated genes in this karyotype, are they enriched on
# this chromosome vs the rest of the genome?
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
  fisher_per_chrom("comp_Buffer")
) %>%
  dplyr::mutate(
    level = factor(level, levels = c("RNA", "Protein", "Buffer")),
    chr   = factor(chr,   levels = chrom_levels),
    karyotype_lab = karyotype_mapping[karyotype],
    karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping)
  )

saveRDS(chrom_fisher, file.path(RES_PATH, "chrom_fisher_per_level.rds"))

# ---- 7.5 Heatmap 1: raw fraction compensated ------------------------------
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

ggsave(file.path(IMG_PATH, "chrom_compensation_heatmap.pdf"), p_chrom_heatmap_raw, width = 11, height = 6)
saveRDS(p_chrom_heatmap_raw, file.path(IMG_PATH, "chrom_compensation_heatmap.rds"))

# ---- 7.6 Heatmap 2: enrichment vs karyotype baseline ----------------------
# Spots where a chromosome buffers more (red) or less (blue) than its
# karyotype-wide average.
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

ggsave(file.path(IMG_PATH, "chrom_enrichment_heatmap.pdf"), p_chrom_heatmap_enr, width = 11, height = 6)
saveRDS(p_chrom_heatmap_enr, file.path(IMG_PATH, "chrom_enrichment_heatmap.rds"))

# ---- 7.7 Heatmap 3: Fisher significance (with stars) ----------------------
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
    log2_OR = log2(pmin(pmax(OR, 0.1), 10)),  # clip for color stability
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

ggsave(file.path(IMG_PATH, "chrom_fisher_heatmap.pdf"), p_chrom_heatmap_fisher, width = 11, height = 6)
saveRDS(p_chrom_heatmap_fisher, file.path(IMG_PATH, "chrom_fisher_heatmap.rds"))

# ---- 7.8 Chromosome-17-specific deep dive ---------------------------------
# Given chr 17 dominates the buffering signal, check whether that signal is
# (a) driven by complex subunits (mechanistic) or (b) just a coverage artifact.

df_chr17 <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::mutate(is_chr17 = !is.na(chr) & chr == "17")

# Are chr 17 buffered genes more often CORUM members than other buffered genes?
# Requires `in_complex` from the CORUM section above; only run if present.
if ("in_complex" %in% colnames(df_classified)) {

  df_chr17 <- df_chr17 %>%
    dplyr::mutate(in_complex = name %in% corum_genes)

  chr17_complex_test <- df_chr17 %>%
    dplyr::filter(comp_Buffer) %>%
    dplyr::group_by(karyotype) %>%
    dplyr::filter(dplyr::n() >= 10,
                  dplyr::n_distinct(is_chr17) == 2) %>%
    dplyr::summarise(
      n_chr17           = sum(is_chr17),
      n_other           = sum(!is_chr17),
      frac_chr17_corum  = mean(in_complex[is_chr17]),
      frac_other_corum  = mean(in_complex[!is_chr17]),
      fisher_p          = fisher.test(table(is_chr17, in_complex))$p.value,
      fisher_OR         = fisher.test(table(is_chr17, in_complex))$estimate,
      .groups = "drop"
    ) %>%
    dplyr::mutate(fisher_padj = p.adjust(fisher_p, method = "BH"))

  saveRDS(chr17_complex_test, file.path(RES_PATH, "chr17_buffer_corum_test.rds"))

  # Also: is chr 17 itself enriched for CORUM members vs the rest of the genome,
  # independent of buffering? Tells us if chr 17 is just complex-rich by nature.
  chr17_corum_baseline <- df_chr17 %>%
    dplyr::distinct(name, is_chr17, in_complex) %>%
    dplyr::summarise(
      frac_corum_chr17  = mean(in_complex[is_chr17]),
      frac_corum_other  = mean(in_complex[!is_chr17]),
      fisher_p          = fisher.test(table(is_chr17, in_complex))$p.value,
      fisher_OR         = fisher.test(table(is_chr17, in_complex))$estimate
    )

  saveRDS(chr17_corum_baseline, file.path(RES_PATH, "chr17_corum_baseline.rds"))

  message("Chr 17 CORUM baseline (all genes): frac_chr17 = ",
          round(chr17_corum_baseline$frac_corum_chr17, 3),
          ", frac_other = ",
          round(chr17_corum_baseline$frac_corum_other, 3),
          ", OR = ",
          round(chr17_corum_baseline$fisher_OR, 2))
}

# ---------------------------------------------------------------------------
# 8. Genome-wide landscape of complex-subunit density per chromosome
# ---------------------------------------------------------------------------

# Use the union of genes in your data (so the baseline is comparable to what
# the buffering analysis tested). Alternatively, use all protein-coding genes
# from biomaRt for a true genome-wide reference.

# ---- 8.1 CORUM density per chromosome --------------------------------------
chrom_complex_landscape <- gene_chrom %>%
  dplyr::mutate(in_complex = name %in% corum_genes) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    n_genes        = dplyr::n(),
    n_complex      = sum(in_complex),
    frac_complex   = mean(in_complex),
    .groups        = "drop"
  ) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_levels))

# Genome-wide baseline for reference
genome_complex_frac <- mean(gene_chrom$name %in% corum_genes)

saveRDS(chrom_complex_landscape, file.path(RES_PATH, "chrom_complex_landscape.rds"))

# ---- 8.2 Fisher test: is each chromosome enriched for CORUM members? -------
chrom_complex_fisher <- gene_chrom %>%
  dplyr::mutate(in_complex = name %in% corum_genes) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    a       = sum(in_complex),
    b       = sum(!in_complex),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    a_other = sum(a) - a,
    b_other = sum(b) - b
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    fisher_p  = fisher.test(matrix(c(a, b, a_other, b_other), nrow = 2))$p.value,
    fisher_OR = fisher.test(matrix(c(a, b, a_other, b_other), nrow = 2))$estimate
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    fisher_padj = p.adjust(fisher_p, method = "BH"),
    chr         = factor(chr, levels = chrom_levels)
  )

saveRDS(chrom_complex_fisher, file.path(RES_PATH, "chrom_complex_fisher.rds"))

# ---- 8.3 Plot: complex density per chromosome ------------------------------
p_complex_landscape <- chrom_complex_landscape %>%
  dplyr::left_join(chrom_complex_fisher %>% dplyr::select(chr, fisher_padj, fisher_OR),
                   by = "chr") %>%
  dplyr::mutate(
    star = dplyr::case_when(
      fisher_padj <= 0.001 ~ "***",
      fisher_padj <= 0.01  ~ "**",
      fisher_padj <= 0.05  ~ "*",
      TRUE                 ~ ""
    ),
    above_baseline = frac_complex > genome_complex_frac
  ) %>%
  ggplot(aes(x = chr, y = frac_complex, fill = above_baseline)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = genome_complex_frac,
             linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = star, y = frac_complex + 0.005), size = 4) +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "steelblue"),
                    guide = "none") +
  theme_bw() +
  labs(x = "Chromosome",
       y = "Fraction of genes in CORUM complexes",
       caption = sprintf("Dashed = genome-wide baseline (%.1f%%). Stars = BH-adjusted Fisher p.",
                         100 * genome_complex_frac))

ggsave(file.path(IMG_PATH, "chrom_complex_landscape.pdf"), p_complex_landscape, width = 10, height = 4)
saveRDS(p_complex_landscape, file.path(IMG_PATH, "chrom_complex_landscape.rds"))

# ---- 8.4 Direct comparison: complex density vs buffering enrichment --------
# The key plot for the chr 17 story.

buffer_enr_per_chrom <- chrom_fisher %>%
  dplyr::filter(level == "Buffer") %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    mean_log2_OR = mean(log2(pmax(OR, 0.1)), na.rm = TRUE),
    any_signif   = any(fisher_padj <= 0.05, na.rm = TRUE),
    .groups      = "drop"
  )

complex_vs_buffer <- chrom_complex_landscape %>%
  dplyr::left_join(buffer_enr_per_chrom, by = "chr") %>%
  dplyr::filter(!is.na(mean_log2_OR))

# Correlation: do complex-rich chromosomes buffer more?
cor_test <- cor.test(complex_vs_buffer$frac_complex,
                     complex_vs_buffer$mean_log2_OR,
                     method = "spearman")

p_complex_vs_buffer <- complex_vs_buffer %>%
  ggplot(aes(x = frac_complex, y = mean_log2_OR)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", alpha = 0.2) +
  geom_point(aes(colour = any_signif, size = n_genes), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = chr), size = 3, max.overlaps = Inf) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                      name = "Any sig.\nbuffering\nenrichment") +
  scale_size_continuous(name = "n genes\nin data") +
  theme_bw() +
  labs(x = "Fraction of chromosome genes in CORUM complexes",
       y = "Mean log2(OR) for buffering enrichment\n(averaged over karyotypes)",
       caption = sprintf("Spearman rho = %.2f, p = %.2g",
                         cor_test$estimate, cor_test$p.value))

ggsave(file.path(IMG_PATH, "complex_density_vs_buffering.pdf"), p_complex_vs_buffer, width = 7, height = 5)
saveRDS(p_complex_vs_buffer, file.path(IMG_PATH, "complex_density_vs_buffering.rds"))

# ---------------------------------------------------------------------------
# 9. Functional enrichment of chr 17 buffered genes
#
#    The picture from earlier: chr 17 buffered genes are dominated by
#    nuclear-encoded mitochondrial proteins (mitoribosomal subunits, OXPHOS
#    assembly factors, SLC25 carriers). This section formalizes that finding
#    and checks whether the mechanism is chr 17-specific or universal.
# ---------------------------------------------------------------------------

# ---- 9.0 Gene sets --------------------------------------------------------
chr17_buffered <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(chr == "17", comp_Buffer) %>%
  dplyr::pull(name) %>% unique()

chr17_all <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(chr == "17") %>%
  dplyr::pull(name) %>% unique()

all_measured <- unique(df_classified$name)

message(sprintf("chr 17 buffered: %d genes | chr 17 all: %d | all measured: %d",
                length(chr17_buffered), length(chr17_all), length(all_measured)))

# ---- 9.1 GO enrichment across BP / MF / CC -------------------------------
# Using chr 17 background — asks what's special about buffered chr 17 genes
# relative to other chr 17 genes (rather than relative to the whole genome,
# which would conflate "is on chr 17" with "is buffered").

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
saveRDS(go_by_ontology, file.path(RES_PATH, "GO_chr17_buffered_all_ontologies.rds"))

# Combined dotplot across ontologies
if (nrow(go_by_ontology) > 0) {
  top_terms <- go_by_ontology %>%
    dplyr::group_by(ontology) %>%
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
    theme_bw() +
    labs(x = "-log10(adj. p)", y = NULL,
         colour = "Ontology", size = "Gene count",
         title = "GO enrichment: chr 17 buffered genes vs chr 17 background")

  ggsave(file.path(IMG_PATH, "GO_chr17_buffered_combined.pdf"), p_go_combined, width = 9, height = 6)
  saveRDS(p_go_combined, file.path(IMG_PATH, "GO_chr17_buffered_combined.rds"))
}

# ---- 9.2 MitoCarta 3.0 enrichment ----------------------------------------
# MitoCarta is the gold-standard, focused inventory of nuclear-encoded
# mitochondrial proteins. Cleaner annotation than GO for this question.
# Download from: https://www.broadinstitute.org/mitocarta/
# Expected columns: Symbol, MitoCarta3.0_List (TRUE/FALSE), MitoCarta3.0_SubMitoLocalization

mitocarta_path <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/data/Human.MitoCarta3.0.xls"


if (file.exists(mitocarta_path)) {
  mitocarta <- read_excel(mitocarta_path, sheet = "A Human MitoCarta3.0") %>%
    dplyr::as_tibble() %>%
    dplyr::pull(Symbol) %>%
    unique()

  # Fisher test per karyotype, restricted to chr 17:
  # "Among chr 17 genes, are MitoCarta members enriched in the buffered set?"
  mito_chr17_test <- df_classified %>%
    dplyr::left_join(gene_chrom, by = "name") %>%
    dplyr::filter(chr == "17") %>%
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
    dplyr::mutate(padj = p.adjust(p, method = "BH"))

  saveRDS(mito_chr17_test, file.path(RES_PATH, "mitocarta_chr17_enrichment.rds"))

  # Same test genome-wide (not restricted to chr 17):
  # "Among ALL genes, are MitoCarta members enriched in the buffered set?"
  mito_genome_test <- df_classified %>%
    dplyr::mutate(in_mito = name %in% mitocarta) %>%
    dplyr::group_by(karyotype) %>%
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
    dplyr::mutate(padj = p.adjust(p, method = "BH"))

  saveRDS(mito_genome_test, file.path(RES_PATH, "mitocarta_genome_enrichment.rds"))

  # Side-by-side plot: enrichment on chr 17 vs genome-wide, per karyotype
  mito_combined <- dplyr::bind_rows(
    mito_chr17_test  %>% dplyr::mutate(scope = "chr 17 only"),
    mito_genome_test %>% dplyr::mutate(scope = "Genome-wide")
  ) %>%
    dplyr::mutate(
      karyotype_lab = karyotype_mapping[karyotype],
      karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping),
      log2_OR       = log2(pmax(OR, 0.1)),
      star = dplyr::case_when(
        padj <= 0.001 ~ "***",
        padj <= 0.01  ~ "**",
        padj <= 0.05  ~ "*",
        TRUE          ~ ""
      )
    )

  p_mito <- ggplot(mito_combined,
                   aes(x = karyotype_lab, y = log2_OR, fill = scope)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = star,
                  y = log2_OR + sign(log2_OR) * 0.15),
              position = position_dodge(width = 0.8), size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw() +
    labs(x = "Karyotype",
         y = "log2(OR) MitoCarta enrichment in buffered set",
         fill = "Scope")

  ggsave(file.path(IMG_PATH, "mitocarta_buffer_enrichment.pdf"), p_mito, width = 7, height = 4)
  saveRDS(p_mito, file.path(IMG_PATH, "mitocarta_buffer_enrichment.rds"))
}

# ---- 9.3 Is mitochondrial buffering chr 17-specific or universal? --------
# Run GO CC per chromosome on the buffered set, looking specifically for
# mitochondrial terms. If only chr 17 shows them, chr 17 is special. If many
# chromosomes do, the mechanism is universal and chr 17 just has more mito
# genes than average.

buffered_by_chrom <- df_classified %>%
  dplyr::left_join(gene_chrom, by = "name") %>%
  dplyr::filter(comp_Buffer, !is.na(chr)) %>%
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

# Heatmap: mitochondrial-related terms across chromosomes
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
           title = "Mitochondrial GO CC enrichment among buffered genes, per chromosome",
           caption = "Cells labeled with gene count. Empty cells = term not significant.")

    ggsave(file.path(IMG_PATH, "mito_GO_per_chrom.pdf"), p_mito_per_chrom, width = 10, height = 5)
    saveRDS(p_mito_per_chrom, file.path(IMG_PATH, "mito_GO_per_chrom.rds"))
  }
}

# ---- 9.4 Concrete illustration: LFC of chr 17 mito genes ------------------
# Pull the actual mitochondrial genes from the chr 17 buffered set and show
# their RNA vs Protein LFC. The classic buffering signature: RNA tracks DNA,
# protein collapses toward zero.

mito_chr17_buffered_genes <- unique(unlist(strsplit(
  go_by_ontology %>%
    dplyr::filter(stringr::str_detect(Description,
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
         title = "Buffering of chr 17 mitochondrial genes",
         caption = "Each line = one gene. Red = median across genes.")

  ggsave(file.path(IMG_PATH, "chr17_mito_LFC_buffering.pdf"),
         p_lfc_mito, width = 8, height = 4)
  saveRDS(p_lfc_mito, file.path(IMG_PATH, "chr17_mito_LFC_buffering.rds"))
}

# ---------------------------------------------------------------------------
# 8. Gain vs loss stratification
# ---------------------------------------------------------------------------

df_classified <- df_classified %>%
  tidyr::separate(karyotype, into = c("A", "B"), sep = ":", remove = FALSE, convert = TRUE) %>%
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
  theme_bw() +
  labs(x = "Copy-number change", y = "Fraction compensated",
       fill = "Level")

ggsave(file.path(IMG_PATH, "fraction_compensated_gain_vs_loss.pdf"), p_frac_direction, width = 6, height = 4)

# ---------------------------------------------------------------------------
# 9. 3-bit class distribution per karyotype
# ---------------------------------------------------------------------------

class_dist <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(karyotype) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(karyotype_lab = karyotype_mapping[karyotype]) %>%
  dplyr::mutate(karyotype_lab = factor(karyotype_lab, levels = karyotype_mapping))

p_class_stack <- ggplot(class_dist, aes(x = karyotype_lab, y = frac, fill = class_label)) +
  geom_col() +
  theme_bw() +
  labs(x = "Karyotype", y = "Fraction of genes", fill = "Compensation class") +
  theme(legend.position = "right")

ggsave(file.path(IMG_PATH, "class_distribution_per_karyotype.pdf"), p_class_stack, width = 8, height = 5)

# ---------------------------------------------------------------------------
# 10. CORUM complex enrichment among buffered genes
# ---------------------------------------------------------------------------

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

  # Fisher test per karyotype: are complex subunits enriched among buffered genes?
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

  # Build annotation table for geom_signif: one bracket per karyotype
  sig_label <- function(p) {
    dplyr::case_when(
      is.na(p)   ~ "",
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE       ~ "ns"
    )
  }

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
      data = pvals_complex,
      aes(
        xmin        = xmin,
        xmax        = xmax,
        annotations = annotation,
        y_position  = y_position,
        group       = karyotype_lab
      ),
      manual      = TRUE,
      textsize    = 3.2,
      tip_length  = 0.01,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "steelblue"),
                      labels = c("Non-complex", "CORUM complex")) +
    theme_bw() +
    labs(x = "Karyotype", y = "Fraction buffered (post-transcriptional)",
         fill = "",
         caption = "Fisher's exact test, BH-adjusted across karyotypes. OR = odds ratio for complex membership among buffered genes.")

  ggsave(file.path(IMG_PATH, "complex_buffering.pdf"), p_complex, width = 7, height = 4)
}

# ---------------------------------------------------------------------------
# 11. GO enrichment per compensation class
# ---------------------------------------------------------------------------

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
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH"
  )

  dir.create(file.path(RES_PATH, "enrichment"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(enrich_res, file.path(RES_PATH, "enrichment", "compensation_class_GO.rds"))

  p_enrich <- dotplot(enrich_res, showCategory = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(IMG_PATH, "GO_by_class.pdf"), p_enrich, width = 10, height = 8)
}

# ---------------------------------------------------------------------------
# 12. Summary table
# ---------------------------------------------------------------------------

summary_tbl <- df_classified %>%
  dplyr::group_by(karyotype, class_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class_label, values_from = n, values_fill = 0)

write_csv(summary_tbl, file.path(RES_PATH, "compensation_class_summary.csv"))
