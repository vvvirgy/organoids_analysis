
rm(list = ls())
source("utils.R")

library(devil)
library(tidyverse)
library(ggplot2)

df_DNA_path = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/cumulative_dp.rds"
df_DNA = readRDS(df_DNA_path)

genes = df_DNA$hgnc_symbol %>% unique()
samples = df_DNA$sample %>% unique()

# Calculate SFS
sfs_baseline = lapply(samples, function(s) {
  df_DNA %>% 
    #dplyr::filter(sample == s, karyotype == "1:1") %>%
    dplyr::filter(sample == s, karyotype %in% c("1:1", "2:0")) %>% 
    dplyr::select(segment_id, segment_DP) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(segment_DP) %>% 
    unlist() %>% 
    mean(na.rm = TRUE)
}) %>% unlist()

# If some samples lack 1:1, you may need to fallback to the global median
sfs_baseline[is.na(sfs_baseline)] <- median(sfs_baseline, na.rm = TRUE)
sfs <- sfs_baseline / exp(mean(log(sfs_baseline)))
names(sfs) = samples

DNA_res = dplyr::tibble()

for (g in genes) {
  idx = which(g == genes)
  if (idx %% 100 == 0) print(idx)
  
  # Fit one gene at the time
  df_sub = df_DNA %>% dplyr::filter(hgnc_symbol == g)
  
  karyotypes = lapply(1:nrow(df_sub), function(i) {
    rep(df_sub$karyotype[i], length(df_sub[i, ]$segment_DP[[1]]))
  }) %>% unlist() 
  karyotypes = factor(karyotypes, levels = c("1:1", setdiff(karyotypes, "1:1")))
  
  clusters =  lapply(1:nrow(df_sub), function(i) {
    rep(df_sub$sample[i], length(df_sub[i, ]$segment_DP[[1]]))
  }) %>% unlist()
  
  current_sfs = sfs[clusters]
  
  counts = lapply(1:nrow(df_sub), function(i) {
    df_sub[i, ]$segment_DP[[1]]
  }) %>% unlist()
  
  design_matrix = model.matrix(~karyotype, dplyr::tibble(karyotype = karyotypes))
  
  fit <- suppressMessages(my_fit_devil(input_matrix = t(as.matrix(counts)), design_matrix = design_matrix, 
                                overdispersion = "MOM", size_factors = current_sfs, max_iter = 500, verbose = TRUE))
  
  # Test DE
  sample_ids = clusters
  coeffs = colnames(design_matrix)
  
  lfc_res = dplyr::tibble()
  for (c in coeffs) {
    contrast_vec = as.numeric(coeffs == c)
    test_res = devil::test_de(devil.fit = fit, contrast = contrast_vec, clusters = sample_ids) %>% dplyr::mutate(name = g)
    
    if (!is.null(test_res)) {
      lfc_res = dplyr::bind_rows(
        lfc_res,
        test_res %>% dplyr::filter(name == g) %>% 
          dplyr::mutate(karyotype = str_replace(c, "karyotype", "")) %>% 
          dplyr::filter(karyotype != "(Intercept)")
      )
    }
  }

  # lfcs = c(fit$beta / log(2))
  # names(lfcs) = str_replace_all(colnames(design_matrix), "karyotype", "")
  # lfcs = lfcs[!grepl("Int", names(lfcs))]
  
  DNA_res = dplyr::bind_rows(
    DNA_res, lfc_res
    # dplyr::tibble(
    #   lfc = lfcs, 
    #   karyotype = names(lfcs), 
    #   name = g
    # )  
  )
}

saveRDS(DNA_res, "results/DNA_lfc.rds")

DNA_res %>% 
  ggplot(mapping = aes(x = karyotype, y = lfc)) +
  geom_boxplot()

DNA_res %>% 
  ggplot(mapping = aes(x = karyotype, y = lfc)) +
  geom_violin()

df_DNA %>% 
  dplyr::select(segment_id, sample, karyotype) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(sample, karyotype) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(f = n / sum(n)) %>% view()



# rna_prot_df = readRDS("results/lfc_prot_and_rna_bind.rds")
# rna_prot_df %>% dplyr::select(lfc, name, karyotype, omic) %>% 
#   dplyr::bind_rows(DNA_res %>% dplyr::mutate(omic = "DNA")) %>% 
#   ggplot(mapping = aes(x = karyotype, y = lfc, fill = omic)) +
#   geom_boxplot(outliers = F) +
#   theme_bw()


