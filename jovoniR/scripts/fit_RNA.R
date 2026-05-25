
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Default values if not provided
# sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
# use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)
# 
# cat(paste0("Running with sf_method: ", sf_method, " | use_stable: ", use_stable, "\n"))
#rm(list = setdiff(ls(), c("sf_method", "use_stable")))

rm(list = ls())
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/utils.R")
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/scripts/constants.R")
library(tidyverse)
library(devil)
library(SingleCellExperiment)
dir.create("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/", recursive = T)

# Load Data
sce = readRDS(SCE_PATH)
#input_data = filter_sce(sce)
input_data = filter_sce_v2(sce, min_counts_floor = 2000)

tibble(
  sample = input_data$meta$sample_id,
  total_counts = colSums(input_data$counts)
) %>%
  ggplot(aes(x = sample, y = total_counts)) +
  geom_violin(scale = "width", fill = "lightblue") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Per-sample depth distribution (post-filter)")

# 1. How many cells per sample after filtering?

karyotypes_df_all = readRDS(META_PATH)
rm(sce)
gc()

# --- Size Factor Calculation ---
all_samples <- unique(input_data$meta$sample_id)
# custom_sf_list <- lapply(all_samples, function(s) {
#   p_cells <- rownames(input_data$meta[input_data$meta$sample_id == s, ])
#   p_karyo <- karyotypes_df_all %>% filter(sample == s, karyotype %in% c("1:1", "2:0"))
#   p_stable_genes <- intersect(p_karyo$hgnc_symbol, rownames(input_data$counts))
#   
#   sf <- devil:::calculate_sf(input_data$counts[p_stable_genes, p_cells], method = SF_METHOD)
#   names(sf) <- p_cells
#   return(sf)
# })


#stable_size_factors <- unlist(custom_sf_list)
#stable_size_factors <- stable_size_factors[colnames(input_data$counts)]
generic_size_factors = devil:::calculate_sf(input_data$counts, method = SF_METHOD)

tibble(
  sample = input_data$meta$sample_id,
  total = colSums(input_data$counts),
  sf = generic_size_factors
) %>%
  ggplot(aes(total, sf, color = sample)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() + scale_y_log10() +
  theme(legend.position = "none")

gc()

# Dynamic path naming
sf_suffix <- paste0(SF_METHOD, "_stable_", USE_STABLE)
#saveRDS(stable_size_factors, paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/stable_sf_", sf_suffix, ".rds"))
saveRDS(generic_size_factors, paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/generic_sf_", sf_suffix, ".rds"))

if (USE_STABLE) {
  size_factors = stable_size_factors
} else {
  size_factors = generic_size_factors
}

# --- LFC Analysis ---
all_genes = rownames(input_data$counts)
gene = all_genes[4]
gene = "REG4"
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
  # if(sum(input_data$counts[gene, cell_ids] > 0) <= length(cell_ids) * 0.1) return(NULL)
  
  # Compute per-sample mean expression (named vector: sample_id -> mean count)
  per_sample_means <- tibble(
    sample_id = meta$sample_id,
    karyotype = meta$karyotype,
    expr = input_data$counts[gene, cell_ids]
  ) %>%
    dplyr::group_by(sample_id, karyotype) %>%
    dplyr::summarise(mean_expr_sample = mean(expr), .groups = "drop")
  
  # Build a list-column: one named vector per karyotype level
  sample_means_by_karyo <- per_sample_means %>%
    dplyr::group_by(karyotype) %>%
    dplyr::summarise(
      sample_means = list(setNames(mean_expr_sample, sample_id)),
      .groups = "drop"
    )
  
  tibble(c = c(input_data$counts[gene, cell_ids]), s = meta$sample_id, sf = size_factors[cell_ids], k = meta$karyotype) %>%
    ggplot(mapping = aes(x = k, y = log1p(c / sf), col = s)) +
    geom_boxplot()
  
  design_matrix = model.matrix(~karyotype, meta)
  
  fit = devil::fit_devil(
    input_matrix = t(as.matrix(input_data$counts[gene, cell_ids])), 
    design_matrix = design_matrix, 
    clusters = meta$sample_id,
    overdispersion = "MOM", 
    size_factors = size_factors[cell_ids], 
    max_iter = 500, 
    verbose = FALSE
  )
  
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
  
  lfc_res %>% 
    dplyr::mutate(mean_expr = mean(input_data$counts[gene, cell_ids])) %>% 
    dplyr::mutate(non_zero_percent = sum(input_data$counts[gene, cell_ids] > 0) / length(cell_ids) * 100) %>% 
    dplyr::left_join(meta %>% 
                       dplyr::select(karyotype, sample_id) %>% 
                       dplyr::distinct() %>% 
                       dplyr::group_by(karyotype) %>% 
                       dplyr::summarise(n_samples = n())) %>% 
    dplyr::left_join(
      sample_means_by_karyo %>% dplyr::mutate(karyotype = as.character(karyotype)),
      by = "karyotype"
    )
})

all_lfc_res = do.call("bind_rows", all_lfc_res)

# Save Final Results
final_output_path <- paste0("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/lfc_res_", sf_suffix, ".rds")
saveRDS(all_lfc_res, final_output_path)
cat(paste0("Process complete. Results saved to: ", final_output_path, "\n"))