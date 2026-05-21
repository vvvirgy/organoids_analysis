
rm(list = ls())
source("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/utils.R")
library(tidyverse)
library(magrittr)
library(devil)
library(SingleCellExperiment)
dir.create("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/", recursive = T)

# Read Inputs
diploid_gene_list = readRDS("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/noise_model/filtered_design_matrix_diploids.rds")
size_factors = readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/generic_sf_psinorm_stable_FALSE.rds")
sce = readRDS(SCE_PATH)
input_data = filter_sce_v2(sce, min_counts_floor = 2000)
rm(sce)
gc()

karyotypes_df_all = readRDS(META_PATH)

names(diploid_gene_list)

good_genes_df = karyotypes_df_all %>% 
  dplyr::filter(karyotype == "1:1") %>% 
  dplyr::select(sample, hgnc_symbol, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(hgnc_symbol) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 10)


# Start analysis
#genes = names(diploid_gene_list)
genes = intersect(names(diploid_gene_list), good_genes_df$hgnc_symbol)
# genes = good_genes_df$hgnc_symbol
# gene = genes[32]

# Precompute ONCE, outside both loops:
sample_to_genes <- karyotypes_df_all %>%
  dplyr::filter(karyotype == "1:1") %>%
  dplyr::select(sample, hgnc_symbol) %>%
  dplyr::distinct()

# Per-gene sample lookup (outside iter loop):
gene_samples_map <- sample_to_genes %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::summarise(samples = list(sample)) %>%
  tibble::deframe()

df_gene_boot_all <- lapply(genes, function(gene) {
  print(gene)
  
  if ((which(genes == gene) %% 100) == 0) print(which(genes == gene) / length(genes) * 100)
  # ---- Per-gene setup (done once) ----
  
  if (!gene %in% rownames(input_data$counts)) return(NULL)
  
  iterations_df = diploid_gene_list[[gene]]
  good_sample_ids = unique(iterations_df$sample)
  
  if (length(good_sample_ids) < 6) return(NULL)
  
  cell_idx <- which(input_data$meta$sample_id %in% good_sample_ids)
  cell_ids <- rownames(input_data$meta)[cell_idx]
  sample_id_per_cell <- input_data$meta$sample_id[cell_idx]
  
  # Extract the gene's expression ONCE
  gene_counts <- as.numeric(input_data$counts[gene, cell_ids])
  sf_sub      <- size_factors[cell_ids]
  
  # Gene-level filters: if they fail, skip the gene entirely (all 50 iters)
  mean_expr <- mean(gene_counts)
  if (mean_expr <= 0.05) return(NULL)
  non_zero_percent <- sum(gene_counts > 0) / length(cell_ids) * 100
  if (non_zero_percent <= 1) return(NULL)
  
  # Per-sample means: gene-level, not iter-level
  means_per_sample <- vapply(
    split(gene_counts, sample_id_per_cell),
    mean, numeric(1)
  )
  min_mean_expr <- min(means_per_sample)
  if (min_mean_expr < 0.05) return(NULL)
  
  # Pre-shape the input matrix for devil (1 x n_cells), done once
  input_mat <- matrix(gene_counts, nrow = 1, dimnames = list(gene, cell_ids))
  n_samples <- length(good_sample_ids)
  
  # ---- Per-iteration loop (only the random group changes) ----
  unique_iters = unique(iterations_df$iteration) %>% sort()
  K = 20
  first_K_iters = unique_iters[1:min(K, length(unique_iters))]
  iters = first_20_iters
  
  cell_id_to_sample_df = dplyr::tibble(sample = input_data$meta$sample_id, cell_id = rownames(input_data$meta))
  
  lapply(iters, function(iter) {
    
    meta_iteration <- iterations_df %>% 
      dplyr::filter(iteration == iter)
    
    # Build meta in the exact order of input_mat columns
    meta_iteration_join <- tibble::tibble(cell_id = colnames(input_mat)) %>%
      dplyr::left_join(cell_id_to_sample_df, by = "cell_id") %>%
      dplyr::left_join(meta_iteration,       by = "sample")
    
    design_matrix <- model.matrix(~group, meta_iteration_join)
    
    fit <- devil::fit_devil(
      input_matrix    = input_mat,
      design_matrix   = design_matrix,
      clusters        = sample_id_per_cell,
      overdispersion  = "MOM",
      size_factors    = sf_sub,
      max_iter        = 500,
      verbose         = FALSE
    )
    
    devil::test_de(fit, contrast = c(0, 1)) %>%
      dplyr::mutate(
        gene             = gene,
        iteration        = iter,
        mean_expr        = mean_expr,
        non_zero_percent = non_zero_percent,
        means_per_sample = list(unname(means_per_sample)),
        min_mean_expr    = min_mean_expr,
        n_samples        = n_samples
      )
  }) %>% dplyr::bind_rows()
  
}) %>% dplyr::bind_rows()


# # Old
# df_gene_boot_all <- lapply(genes, function(gene) {
#   print(gene)
#   
#   if ((which(genes == gene) %% 100) == 0) print(which(genes == gene) / length(genes) * 100)
#   # ---- Per-gene setup (done once) ----
#   
#   if (!gene %in% rownames(input_data$counts)) return(NULL)
#   
#   good_sample_ids <- gene_samples_map[[gene]]   # precomputed outside, see prev message
#   if (length(good_sample_ids) < 6) return(NULL)
#   
#   cell_idx <- which(input_data$meta$sample_id %in% good_sample_ids)
#   cell_ids <- rownames(input_data$meta)[cell_idx]
#   sample_id_per_cell <- input_data$meta$sample_id[cell_idx]
#   
#   # Extract the gene's expression ONCE
#   gene_counts <- as.numeric(input_data$counts[gene, cell_ids])
#   sf_sub      <- size_factors[cell_ids]
#   
#   # Gene-level filters: if they fail, skip the gene entirely (all 50 iters)
#   mean_expr <- mean(gene_counts)
#   if (mean_expr <= 0.05) return(NULL)
#   non_zero_percent <- sum(gene_counts > 0) / length(cell_ids) * 100
#   if (non_zero_percent <= 1) return(NULL)
#   
#   # Per-sample means: gene-level, not iter-level
#   means_per_sample <- vapply(
#     split(gene_counts, sample_id_per_cell),
#     mean, numeric(1)
#   )
#   min_mean_expr <- min(means_per_sample)
#   
#   if (min_mean_expr < 0.05) return(NULL)
#   
#   # Pre-shape the input matrix for devil (1 x n_cells), done once
#   input_mat <- matrix(gene_counts, nrow = 1, dimnames = list(gene, cell_ids))
#   
#   n_samples <- length(good_sample_ids)
#   
#   # ---- Per-iteration loop (only the random group changes) ----
#   
#   lapply(iters, function(iter) {
#     
#     group_by_sample <- integer(n_samples)
#     group_by_sample[sample.int(n_samples, size = floor(n_samples / 2))] <- 1L
#     
#     names(group_by_sample) <- good_sample_ids
#     group_per_cell  <- group_by_sample[sample_id_per_cell]
#     
#     # Build design matrix directly — skip model.matrix + join
#     design_matrix <- cbind(`(Intercept)` = 1, group = group_per_cell)
#     
#     fit <- devil::fit_devil(
#       input_matrix    = input_mat,
#       design_matrix   = design_matrix,
#       clusters        = sample_id_per_cell,
#       overdispersion  = "MOM",
#       size_factors    = sf_sub,
#       max_iter        = 500,
#       verbose         = FALSE
#     )
#     
#     devil::test_de(fit, contrast = c(0, 1)) %>%
#       dplyr::mutate(
#         gene             = gene,
#         iteration        = iter,
#         mean_expr        = mean_expr,
#         non_zero_percent = non_zero_percent,
#         means_per_sample = list(unname(means_per_sample)),
#         min_mean_expr    = min_mean_expr,
#         n_samples        = n_samples
#       )
#   }) %>% dplyr::bind_rows()
#   
# }) %>% dplyr::bind_rows()

# gene = "MMP19"
# df_gene_boot_all = lapply(genes, function(gene) {
#   print(gene)
#   # gene_iterations_df = diploid_gene_list[[gene]]
# 
#   # iters = unique(gene_iterations_df$iteration)
#   # iters = iters[1:(min(10, length(unique(gene_iterations_df$iteration))))]
#   # 
#   
#   df_gene_boot = lapply(iters, function(iter) {
#     # gene_df = gene_iterations_df %>% dplyr::filter(iteration == iter)
#     
#     good_sample_ids = karyotypes_df_all %>% 
#       dplyr::filter(karyotype == "1:1", hgnc_symbol == gene) %>% 
#       dplyr::select(sample, hgnc_symbol, karyotype) %>% 
#       dplyr::distinct() %>% 
#       dplyr::pull(sample)
#     
#     #good_sample_ids = unique(gene_df$sample)
#     
#     if (length(good_sample_ids) < 6) return(NULL)
#     
#     cell_ids = rownames(input_data$meta[input_data$meta$sample_id %in% good_sample_ids,])
#     
#     sample_groups_df = dplyr::tibble(sample = good_sample_ids, group = rbinom(length(good_sample_ids), 1, prob = 0.5))
#     
#     meta = input_data$meta[cell_ids, ] %>% 
#       as_tibble() %>% 
#       #dplyr::left_join(gene_df %>% dplyr::rename(sample_id = sample), by = "sample_id")
#       dplyr::left_join(sample_groups_df %>% dplyr::rename(sample_id = sample), by = "sample_id")
#     
#     if (!gene %in% rownames(input_data$counts)) return(NULL)
#     #if (mean(input_data$counts[gene, cell_ids]) <= 0.005) return(NULL)
#     if (mean(input_data$counts[gene, cell_ids]) <= 0.05) return(NULL)
#     if (sum(input_data$counts[gene, cell_ids] > 0) / length(cell_ids) <= 0.01) return(NULL)
#     
#     design_matrix = model.matrix(~group, meta)
#     
#     tibble(
#       c = input_data$counts[gene, cell_ids],
#       sample = meta$sample_id,
#       sf = size_factors[cell_ids],
#       group = meta$group
#     ) %>%
#       ggplot(mapping = aes(x = sample, y = c / sf)) +
#       geom_boxplot()
#     
#     means_per_sample = tibble(
#       c = input_data$counts[gene, cell_ids],
#       sample = meta$sample_id,
#       sf = size_factors[cell_ids],
#       group = meta$group
#     ) %>% dplyr::group_by(sample) %>% 
#       dplyr::summarise(mean_expr = mean(c))
#     min_mean_expr = min(means_per_sample$mean_expr)
#     
#     fit = devil::fit_devil(
#       input_matrix = t(as.matrix(input_data$counts[gene, cell_ids])), 
#       design_matrix = design_matrix, 
#       clusters = meta$sample_id,
#       overdispersion = "MOM", 
#       size_factors = size_factors[cell_ids], 
#       max_iter = 500, 
#       verbose = FALSE
#     )
#     
#     devil::test_de(fit, contrast = c(0,1)) %>% dplyr::mutate(gene = gene, iteration = iter) %>% 
#       dplyr::mutate(mean_expr = mean(input_data$counts[gene, cell_ids])) %>% 
#       dplyr::mutate(non_zero_percent = sum(input_data$counts[gene, cell_ids] > 0) / length(cell_ids) * 100) %>% 
#       dplyr::mutate(means_per_sample = list(means_per_sample$mean_expr)) %>% 
#       dplyr::mutate(min_mean_expr = min_mean_expr) %>% 
#       dplyr::mutate(n_samples = length(good_sample_ids))
#   }) %>% do.call("bind_rows", .)
#   
# }) %>% do.call("bind_rows", .)

saveRDS(df_gene_boot_all, "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/diploid_bootstrap.rds")

# df_gene_boot_all = readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/diploid_bootstrap.rds")
# 
# df_gene_boot_all$n_samples = lapply(df_gene_boot_all$gene, function(g) {
#   length(unique(diploid_gene_list[[g]]$sample))
# }) %>% unlist()
# 
# df_gene_boot_all %>% 
#   dplyr::filter(min_mean_expr > 0.5) %>% 
#   dplyr::filter(non_zero_percent > 1) %>% 
#   dplyr::filter(n_samples > 5) %>% 
#   ggplot(mapping = aes(x = lfc)) +
#   geom_density()
# 
# saveRDS(df_gene_boot_all, "/orfeo/cephfs/scratch/cdslab/gsantacatterina/organoids_analysis/jovoniR/results/RNA/diploid_bootstrap.rds")

