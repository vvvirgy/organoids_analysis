# utils for proteomics analysis

# compute TIC normalization and remove NAs, then rescale everything in log2 
compute_tic = function(raw_data, gene_column, log2 = TRUE, filter_NA = TRUE) {
  # Compute TIC per sample
  tic_factors <- colSums(raw_data, na.rm = TRUE)
  
  # Normalize using TIC
  normalized_data <- sweep(raw_data, 2, tic_factors, FUN = "/") * median(tic_factors)
  
  proteomics_data_norm = normalized_data %>% 
    tibble::rownames_to_column(gene_column) %>% 
    reshape2::melt() %>% 
    dplyr::rename(norm_intensity = value)
  
  if(log2) {
    proteomics_data_norm = proteomics_data_norm %>% 
      dplyr::mutate(norm_intensity = log(norm_intensity, base = 2))
  } 
  
  if(filter_NA) {
    proteomics_data_norm = proteomics_data_norm %>% 
    dplyr::filter(!is.na(norm_intensity)) 
  }
  
  return(proteomics_data_norm)
}

# # ---------------------------------------------------------------------------------------------------------
# # filter to remove genes with small segments or that are highly fragmented --> reduce genes with multiCNA
# filter_cna = function(genes_cnas, gene_column, genes_of_interest, min_length = 10^6) {
#   genes_cnas %>% 
#     dplyr::filter(gene_column %in% genes_of_interest) %>%
#     dplyr::mutate(lenght = to - from) %>% 
#     dplyr::filter(length > min_length)
# }

# ---------------------------------------------------------------------------------------------------------

# PLOTS

# colors = setNames(Polychrome::kelly.colors(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)
# colors = setNames(ggthemes::calc_pal()(length(x$karyotype %>% unique)), nm = x$karyotype %>% unique)


