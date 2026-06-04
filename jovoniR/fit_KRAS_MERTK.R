
args <- commandArgs(trailingOnly = TRUE)

# Default values if not provided
sf_method  <- ifelse(length(args) >= 1, args[1], "psinorm")
use_stable <- ifelse(length(args) >= 2, as.logical(args[2]), TRUE)

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

sf_suffix <- paste0(sf_method, "_stable_", use_stable)
if (use_stable) {
  size_factors = readRDS(paste0("results/RNA/stable_sf_", sf_suffix, ".rds"))
} else {
  size_factors = readRDS(paste0("results/RNA/generic_sf_", sf_suffix, ".rds"))
}

# --- LFC Analysis ---
#all_genes = rownames(input_data$counts)
all_genes = c("KRAS", "MERTK")

# Does KRAS  amplified increases KRAS expression? ####
gene = "MERTK"

karyotypes_df = karyotypes_df_all %>% dplyr::filter(gene == hgnc_symbol)
kras_status_df = karyotypes_df_all %>% 
  dplyr::filter(hgnc_symbol == "KRAS") %>% 
  dplyr::mutate(status = mut_status) %>% 
  dplyr::select(sample, status)

karyotypes_df = karyotypes_df %>% dplyr::left_join(kras_status_df) %>% na.omit()

if (nrow(karyotypes_df) == 0) return(NULL)

good_sample_ids = intersect(unique(karyotypes_df$sample), unique(kras_status_df$sample))
cell_ids = rownames(input_data$meta[input_data$meta$sample_id %in% good_sample_ids,])

meta = input_data$meta[cell_ids, ] %>% 
  as_tibble() %>% 
  dplyr::left_join(karyotypes_df %>% dplyr::rename(sample_id = sample), by = "sample_id")

if (!"1:1" %in% unique(meta$karyotype)) return(NULL)

meta$karyotype = factor(meta$karyotype, levels = c("1:1", setdiff(meta$karyotype, "1:1")))
if (length(levels(meta$karyotype)) == 1) return(NULL)
if (mean(input_data$counts[gene, cell_ids]) <= 0.005) return(NULL)

design_matrix = model.matrix(~karyotype, meta)

dplyr::tibble(counts = input_data$counts[gene, cell_ids], sf = size_factors[cell_ids], karyotype = meta$karyotype, status = meta$status) %>% 
  ggplot(mapping = aes(x = karyotype, y = log1p(counts / sf), fill = status, col = status)) +
  geom_boxplot() +
  scale_y_continuous(transform = "log1p")

# Does MERTK amplified increases KRAS expression? ####
gene = "KRAS"

karyotypes_df = karyotypes_df_all %>% dplyr::filter(gene == hgnc_symbol)
mertk_status_df = karyotypes_df_all %>% dplyr::filter(hgnc_symbol == "MERTK") %>% 
  tidyr::separate(karyotype, into = c("A", "B"), sep = ":", convert = TRUE) %>% 
  dplyr::mutate(ploidy = A + B) %>% 
  dplyr::mutate(MERTK_status = ifelse(ploidy == 2, "diploid", ifelse(ploidy < 2, "loss", "gain"))) %>% 
  dplyr::select(sample, MERTK_status)

mertk_status_df$MERTK_status[is.na(mertk_status_df$MERTK_status)]

karyotypes_df = karyotypes_df %>% dplyr::left_join(mertk_status_df) %>% na.omit()

if (nrow(karyotypes_df) == 0) return(NULL)

length(good_sample_ids)
good_sample_ids = intersect(unique(karyotypes_df$sample), unique(mertk_status_df$sample))
cell_ids = rownames(input_data$meta[input_data$meta$sample_id %in% good_sample_ids,])

meta = input_data$meta[cell_ids, ] %>% 
  as_tibble() %>% 
  dplyr::left_join(karyotypes_df %>% dplyr::rename(sample_id = sample), by = "sample_id")

if (!"1:1" %in% unique(meta$karyotype)) return(NULL)

meta$karyotype = factor(meta$karyotype, levels = c("1:1", setdiff(meta$karyotype, "1:1")))
if (length(levels(meta$karyotype)) == 1) return(NULL)
if (mean(input_data$counts[gene, cell_ids]) <= 0.005) return(NULL)

design_matrix = model.matrix(~karyotype, meta)

dplyr::tibble(KRAS_counts = input_data$counts[gene, cell_ids], sf = size_factors[cell_ids], karyotype = meta$karyotype, MERTK_status = meta$MERTK_status) %>% 
  ggplot(mapping = aes(x = karyotype, y = KRAS_counts / sf, fill = MERTK_status)) +
  geom_boxplot()

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

dplyr::tibble(lfc = lfcs, karyotype = names(lfcs), name = gene)