## ============================================================================
## Differential protein expression across tetraploid-vs-diploid karyotype groups
##
## Key optimizations vs. the original DEP-based version:
##   - vsn normalization is performed ONCE on the full proteomics matrix
##     (instead of once per gene x iteration). All downstream tests reuse
##     the same normalized SummarizedExperiment.
##   - Differential testing uses a direct limma call (lmFit + eBayes) on a
##     1-row matrix instead of DEP::test_diff. test_diff internally calls
##     fdrtool to estimate a null model, which fails on a single test
##     statistic ("'xmin' not less than 'xmax'"). Multiple-testing correction
##     is handled downstream across the full grid of gene x iteration tests.
## ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(DEP)
  library(SummarizedExperiment)
  library(limma)
})

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj')

## --- External sources kept as-is (not the bottleneck) -----------------------
source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/constants.R")
source("/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/organoids_analysis/R/functions_utils/make_groups.R")

## ============================================================================
## Utility functions (inlined from proteomics_de_utils.R)
## ============================================================================

#' Preprocess proteomics data into a normalized SummarizedExperiment.
#' Builds the DEP SE, filters proteins with no valid values, then vsn-normalizes.
#' Intended to be called ONCE globally on the full sample set.
dep_preprocesing <- function(data, design) {
  data <- data %>%
    dplyr::select(PG.ProteinGroups, PG.Genes, any_of(design$label))
  
  data <- make_unique(data, names = 'PG.Genes', ids = 'PG.ProteinGroups', delim = '_')
  samples_index <- which(!colnames(data) %in% c('PG.ProteinGroups', 'PG.Genes', 'name', 'ID'))
  
  # log-transformed inside make_se
  dep <- DEP::make_se(proteins_unique = data, columns = samples_index, expdesign = design)
  
  # remove proteins with no valid value in any condition
  dep_filt <- filter_missval(dep, thr = 0)
  
  # variance-stabilizing normalization (the expensive step we want to do ONCE)
  dep_norm <- normalize_vsn(dep_filt)
  
  return(dep_norm)
}

#' Run limma-based differential expression on a 1-row SE.
#'
#' Why not DEP::test_diff? test_diff internally calls fdrtool to compute
#' adjusted p-values; fdrtool fits a null model on the distribution of
#' test statistics and fails when given only one ("'xmin' not less than
#' 'xmax'"). A direct lmFit + eBayes call works for any number of rows.
#'
#' Note: eBayes variance shrinkage is only meaningful when there are many
#' rows to shrink toward. On a 1-row fit, the shrinkage is a no-op and the
#' moderated t is effectively an ordinary t-statistic.
#'
#' Returns a named list with the contrast result (gene, lfc, pvalue).
#' Multiple-testing correction across all gene x iteration tests should be
#' done downstream.
run_limma_one_gene <- function(dep_sub, ctrl = 'A') {
  expr <- SummarizedExperiment::assay(dep_sub)   # 1 x n_samples (vsn-normalized)
  cond <- factor(dep_sub$condition)
  cond <- relevel(cond, ref = ctrl)
  
  design <- stats::model.matrix(~ cond)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  
  # coefficient 2 = "non-ctrl vs ctrl" (since ctrl is the reference level)
  list(
    gene   = SummarizedExperiment::rowData(dep_sub)$name,
    lfc    = unname(fit$coefficients[, 2]),
    pvalue = unname(fit$p.value[, 2])
  )
}

## ============================================================================
## Paths
## ============================================================================

IMG_PATH <- file.path(res_path, 'noise_model')
RES_PATH <- file.path(data_path, 'noise_model')

dir.create(IMG_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(data_path, 'noise_model/protein'), recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## Load inputs
## ============================================================================

# sample-name dictionary
dict <- readRDS(file.path(data_path, 'utilities/full_dict_dna_rna_prot.rds'))

# per-gene karyotype design matrices (list keyed by gene)
design_matrix <- readRDS(file.path(RES_PATH, 'filtered_design_matrix_diploids.rds'))

# proteomics <-> genomics sample-name mapping
samples_check <- dict %>%
  dplyr::select(proteomics_code, DNA, fixed_name) %>%
  dplyr::distinct()

# raw proteomics intensities
data <- readxl::read_excel(
  file.path(data_path, 'proteomics/Results_Organoids_NoIsoforms.xlsx'),
  sheet = 1
) %>%
  as.data.frame()

data[which(is.na(data[, 'PG.Genes'])), 'PG.Genes'] <- 'Unknown'

## ============================================================================
## Reshape proteomics data and build sample annotation
## ============================================================================

# long -> coerce intensities -> wide (same as original)
data <- data %>%
  tidyr::pivot_longer(
    cols = colnames(data)[-(1:2)],
    names_to = 'sample',
    values_to = 'Intensity'
  ) %>%
  dplyr::mutate(Intensity = as.numeric(Intensity)) %>%
  dplyr::mutate(Intensity = ifelse(is.nan(Intensity), NA, Intensity))

# annotation: link replicate label -> PDO -> fixed_name used in design_matrix
ann <- tibble(
  PDO = data$sample,
  replicate = data$sample
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(PDO = gsub('_a$|_b$', '', PDO)) %>%
  dplyr::mutate(PDO = gsub('_HSR', 'HSR', PDO)) %>%
  dplyr::full_join(samples_check, by = dplyr::join_by('PDO' == 'proteomics_code')) %>%
  dplyr::filter(!is.na(PDO), !is.na(replicate))

# wide format for DEP
data <- data %>%
  tidyr::pivot_wider(names_from = sample, values_from = Intensity)

## ============================================================================
## ONE-TIME GLOBAL NORMALIZATION
##   vsn depends only on which samples are in the matrix, not on group labels
##   or which gene is being tested. We normalize the full proteomics matrix
##   once and reuse the result for every gene x iteration below.
## ============================================================================

message("Building global normalized SummarizedExperiment...")

# All samples that exist in `data` AND have a valid annotation
all_sample_labels <- intersect(colnames(data), ann$replicate)

# Minimal "design" only used to satisfy make_se; condition is a placeholder
# and is overridden per-iteration on the subset SE before differential testing.
global_design <- ann %>%
  dplyr::filter(replicate %in% all_sample_labels) %>%
  dplyr::distinct(replicate, .keep_all = TRUE) %>%
  dplyr::transmute(
    label = replicate,
    condition = "all",
    replicate = seq_len(dplyr::n())
  )

dep_global <- dep_preprocesing(data, global_design)

message(sprintf("Global SE: %d proteins x %d samples",
                nrow(dep_global), ncol(dep_global)))

## ============================================================================
## Per-gene / per-iteration differential expression
##   Reuses dep_global. For each iteration we:
##     1. Build the iteration's group assignments
##     2. Subset dep_global to (gene row) x (iteration's samples)
##     3. Re-attach the iteration's condition labels onto colData
##     4. Run lmFit + eBayes on the 1-row matrix (essentially free)
##
##   We keep at most the first 20 unique iterations per gene (fewer if
##   that gene has fewer available). Each gene contributes a tibble with
##   columns: gene, lfc, pvalue, iteration. All gene tibbles are row-bound
##   at the end.
## ============================================================================

MAX_ITERATIONS <- 20

diploid_genes <- names(design_matrix)

per_gene_results <- parallel::mclapply(diploid_genes, function(g) {
  
  message(g)
  
  des_all <- design_matrix[[g]] %>%
    dplyr::rename(hgnc_symbol = gene)
  
  # take at most the first 20 unique iterations (in their original order)
  iterations <- unique(des_all$iteration)
  iterations <- head(iterations, MAX_ITERATIONS)
  
  # locate this gene's row in the global SE once
  gene_idx <- which(rowData(dep_global)$name == g)
  if (length(gene_idx) == 0) {
    # gene was filtered out by filter_missval globally -> nothing to test
    return(NULL)
  }
  
  rows <- lapply(iterations, function(i) {
    
    des_i <- des_all %>%
      dplyr::filter(iteration == i)
    
    des_matrix <- dplyr::full_join(
      des_i, ann,
      by = dplyr::join_by('sample' == 'fixed_name')
    ) %>%
      dplyr::filter(!is.na(group), !is.na(PDO))
    
    # build per-iteration design (label / condition / replicate index)
    design_mat <- des_matrix %>%
      dplyr::group_by(group) %>%
      dplyr::rename(replicate_name = replicate) %>%
      dplyr::filter(!is.na(replicate_name)) %>%
      dplyr::mutate(replicate = seq_len(dplyr::n())) %>%
      dplyr::rename(label = replicate_name) %>%
      dplyr::rename(condition = group) %>%
      dplyr::ungroup()
    
    if (length(unique(design_mat$condition)) < 2) return(NULL)
    
    # match samples by colData(dep_global)$label (NOT colnames)
    keep_mask <- dep_global$label %in% design_mat$label
    if (sum(keep_mask) < 2) return(NULL)
    
    dep_sub <- dep_global[gene_idx, keep_mask]
    
    # re-attach this iteration's condition labels onto colData
    condition_map <- setNames(design_mat$condition, design_mat$label)
    replicate_map <- setNames(design_mat$replicate, design_mat$label)
    dep_sub$condition <- unname(condition_map[dep_sub$label])
    dep_sub$replicate <- unname(replicate_map[dep_sub$label])
    
    if (length(unique(dep_sub$condition)) < 2) return(NULL)
    
    res <- tryCatch(
      run_limma_one_gene(dep_sub, ctrl = 'A'),
      error = function(e) {
        message(sprintf("    limma failed for %s iter %s: %s", g, i, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    
    tibble::tibble(
      gene      = res$gene,
      lfc       = res$lfc,
      pvalue    = res$pvalue,
      iteration = i
    )
  })
  
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  
  dplyr::bind_rows(rows)
  
}, mc.cores = 4)

## ============================================================================
## Combine and save
## ============================================================================

# drop genes with no testable iterations and row-bind into a single tibble
per_gene_results <- Filter(Negate(is.null), per_gene_results)
results_df <- dplyr::bind_rows(per_gene_results)

message(sprintf("Final table: %d rows across %d genes",
                nrow(results_df),
                dplyr::n_distinct(results_df$gene)))

saveRDS(
  results_df,
  file.path(data_path, 'noise_model/protein/dep_diploid_genes_results.rds')
)

message("Done.")