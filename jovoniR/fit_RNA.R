
args <- commandArgs(trailingOnly = TRUE)

# Default values if not provided
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)

cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))

rm(list = setdiff(ls(), c("sf_method", "use_stable")))
source("utils.R")
library(tidyverse)
library(devil)
library(SingleCellExperiment)
dir.create("results/RNA", recursive = T)

# Load Data
sce = readRDS(SCE_PATH)

filter_sce = function(sce) {
  meta = sce@colData
  counts = sce@assays@data$counts
  
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)
  
  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  feat_mad_filter <- total_features > 5 * mad(total_features)
  
  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .2
  cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter | mit_prop_filter
  
  counts = as.matrix(counts)
  counts <- counts[, !cell_outliers_filter]
  meta <- meta[!cell_outliers_filter, ]
  
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]
  counts <- counts[,colnames(counts) %in% rownames(meta)]
  
  list(counts = counts, meta = meta)
}

input_data = filter_sce(sce)
karyotypes_df_all = readRDS(META_PATH)
rm(sce)

# --- Size Factor Calculation ---
all_samples <- unique(input_data$meta$sample_id)
custom_sf_list <- lapply(all_samples, function(s) {
  p_cells <- rownames(input_data$meta[input_data$meta$sample_id == s, ])
  p_karyo <- karyotypes_df_all %>% filter(sample == s, karyotype %in% c("1:1", "2:0"))
  p_stable_genes <- intersect(p_karyo$hgnc_symbol, rownames(input_data$counts))
  
  sf <- devil:::calculate_sf(input_data$counts[p_stable_genes, p_cells], method = "normed_sum")
  names(sf) <- p_cells
  return(sf)
})

stable_size_factors <- unlist(custom_sf_list)
stable_size_factors <- stable_size_factors[colnames(input_data$counts)]
generic_size_factors = devil:::calculate_sf(input_data$counts, method = sf_method)

# Dynamic path naming
sf_suffix <- paste0(sf_method, "_stable_", use_stable)
saveRDS(stable_size_factors, paste0("results/RNA/stable_sf_", sf_suffix, ".rds"))
saveRDS(generic_size_factors, paste0("results/RNA/generic_sf_", sf_suffix, ".rds"))

if (use_stable) {
  size_factors = stable_size_factors
} else {
  size_factors = generic_size_factors
}

# --- LFC Analysis ---
all_genes = rownames(input_data$counts)
gene = "SAMD11"
all_lfc_res = lapply(all_genes, function(gene) {
  print(gene)
  karyotypes_df = karyotypes_df_all %>% dplyr::filter(gene == hgnc_symbol)
  if (nrow(karyotypes_df) == 0) return(NULL)
  
  good_sample_ids = unique(karyotypes_df$sample)
  cell_ids = rownames(input_data$meta[input_data$meta$sample_id %in% good_sample_ids,])
  
  meta = input_data$meta[cell_ids, ] %>% 
    as_tibble() %>% 
    dplyr::left_join(karyotypes_df %>% dplyr::rename(sample_id = sample), by = "sample_id")
  
  if (!"1:1" %in% unique(meta$karyotype)) return(NULL)
  
  meta$karyotype = factor(meta$karyotype, levels = c("1:1", setdiff(meta$karyotype, "1:1")))
  if (length(levels(meta$karyotype)) == 1) return(NULL)
  if (mean(input_data$counts[gene, cell_ids]) <= 0.005) return(NULL)
  
  design_matrix = model.matrix(~karyotype, meta)
  
  fit = suppressMessages(my_fit_devil(
    input_matrix = t(as.matrix(input_data$counts[gene, cell_ids])), 
    design_matrix = design_matrix, 
    overdispersion = "MOM", 
    size_factors = size_factors[cell_ids], 
    max_iter = 500, 
    verbose = FALSE
  ))
  
  lfcs = c(fit$beta / log(2))
  names(lfcs) = str_replace_all(colnames(design_matrix), "karyotype", "")
  lfcs = lfcs[!grepl("Int", names(lfcs))]
  
  # Test DE
  sample_ids = input_data$meta$sample_id
  coeffs = colnames(design_matrix)
  
  lfc_res = dplyr::tibble()
  for (c in coeffs) {
    contrast_vec = as.numeric(coeffs == c)
    test_res = devil::test_de(devil.fit = fit, contrast = contrast_vec, clusters = sample_ids) %>% dplyr::mutate(name = gene)
    
    if (!is.null(test_res)) {
      lfc_res = dplyr::bind_rows(
        lfc_res,
        test_res %>% dplyr::filter(name == gene) %>% dplyr::mutate(karyotype = str_replace(c, "karyotype", "")) %>% dplyr::filter(karyotype != "(Intercept)")
      )
    }
  }
  
  # dplyr::tibble(lfc = lfcs, karyotype = names(lfcs), name = gene)
  lfc_res
})

all_lfc_res = do.call("bind_rows", all_lfc_res)

# Save Final Results
final_output_path <- paste0("results/RNA/lfc_res_", sf_suffix, ".rds")
saveRDS(all_lfc_res, final_output_path)
cat(paste0("Process complete. Results saved to: ", final_output_path, "\n"))