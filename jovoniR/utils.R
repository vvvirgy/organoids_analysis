
#SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_smad2_karyo_all_organoids.rds"
SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/processed_data/scRNA/scexp_karyo_all_organoids_filt.rds"
META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

filter_sce = function(sce) {
  meta = sce@colData
  counts = sce@assays@data$counts
  
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)
  
  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  
  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .2
  cell_outliers_filter <- mad5_filter | feat100_filter | mit_prop_filter
  
  counts = as.matrix(counts)
  counts <- counts[, !cell_outliers_filter]
  meta <- meta[!cell_outliers_filter, ]
  
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]
  counts <- counts[,colnames(counts) %in% rownames(meta)]
  
  list(counts = counts, meta = meta)
}

filter_sce_v2 = function(sce, 
                         min_counts_floor = 500,
                         min_features = 200,
                         mad_high = 10,
                         mad_low = 2,
                         mito_threshold = 0.2,
                         min_gene_mean = 0.01) {
  meta = sce@colData
  counts = sce@assays@data$counts
  
  # Per-cell QC metrics
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)
  mitochondrial_genes <- grepl("^MT-", rownames(counts))
  mito_prop <- colSums(counts[mitochondrial_genes, ]) / total_counts
  
  # Per-sample MAD-based filters on log10(counts)
  log_counts <- log10(total_counts + 1)
  sample_ids <- meta$sample_id
  
  qc_df <- tibble(
    cell = colnames(counts),
    sample_id = sample_ids,
    total_counts = total_counts,
    total_features = total_features,
    log_counts = log_counts,
    mito_prop = mito_prop
  ) %>%
    group_by(sample_id) %>%
    mutate(
      med_log = median(log_counts),
      mad_log = mad(log_counts),
      high_count = log_counts > med_log + mad_high * mad_log,
      low_count  = log_counts < med_log - mad_low  * mad_log
    ) %>%
    ungroup()
  
  # Absolute thresholds (applied globally, not per-sample)
  absolute_low_count   <- total_counts < min_counts_floor
  low_feature_filter   <- total_features < min_features
  high_mito_filter     <- mito_prop > mito_threshold
  
  # Combine all cell-level filters
  cell_outliers_filter <- qc_df$high_count | qc_df$low_count |
    absolute_low_count | low_feature_filter |
    high_mito_filter
  
  # Report what's being removed
  cat("Cells removed by filter:\n")
  cat(sprintf("  High count (per-sample MAD): %d\n", sum(qc_df$high_count)))
  cat(sprintf("  Low count  (per-sample MAD): %d\n", sum(qc_df$low_count)))
  cat(sprintf("  Below absolute count floor:  %d\n", sum(absolute_low_count)))
  cat(sprintf("  Below feature floor:         %d\n", sum(low_feature_filter)))
  cat(sprintf("  High mito proportion:        %d\n", sum(high_mito_filter)))
  cat(sprintf("  Total cells removed:         %d / %d (%.1f%%)\n",
              sum(cell_outliers_filter), length(cell_outliers_filter),
              100 * mean(cell_outliers_filter)))
  
  # Apply cell filter
  counts <- as.matrix(counts)
  counts <- counts[, !cell_outliers_filter]
  meta   <- meta[!cell_outliers_filter, ]
  
  # Gene-level filter
  non_expressed_genes <- rowMeans(counts) <= min_gene_mean
  counts <- counts[!non_expressed_genes, ]
  cat(sprintf("Genes removed (mean <= %g): %d / %d\n",
              min_gene_mean, sum(non_expressed_genes), length(non_expressed_genes)))
  
  # Ensure alignment
  counts <- counts[, colnames(counts) %in% rownames(meta)]
  
  list(counts = counts, meta = meta)
}

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
