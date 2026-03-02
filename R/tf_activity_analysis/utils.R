filter_sce = function(data) {
  meta = data@meta.data
  counts = data@assays$RNA$counts
  
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)
  
  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  feat_mad_filter <- total_features > 5 * mad(total_features)
  
  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .2
  cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter |  mit_prop_filter
  
  counts = as.matrix(counts)
  
  counts <- counts[, !cell_outliers_filter]
  meta <- meta[!cell_outliers_filter, ]
  
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]
  counts <- counts[,colnames(counts) %in% rownames(meta)]
  
  input_data = list(counts = counts, meta = meta)
  
  input_data
}
