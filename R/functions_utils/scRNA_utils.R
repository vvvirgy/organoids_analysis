# utils for scRNAs

create_sc_expr = function(seurat) {
  
  # Extract raw counts and metadata to create SingleCellExperiment object
  counts <- seurat@assays$RNA@counts 
  
  metadata <- seurat@meta.data
  
  # Set up metadata as desired for aggregation and DE analysis
  metadata$cluster_id <- factor(seurat@active.ident)
  
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  
  return(sce)
}

aggregate_cells = function(x) {
  groups_expr = colData(x)[,c('sample')]
  
  tt = aggregate.Matrix(t(counts(x)), groupings = groups_expr, fun = 'sum')
  
  return(tt)
}

generate_pseudobulk_count_matrix = function(aggregated) {
  
  pb = lapply(aggregated, function(x) {
    tibble(genes = colnames(x), 
           count = x %>% as.numeric) %>% 
      dplyr::rename_with(., 
                         ~ rownames(x), 
                         dplyr::starts_with('count'))
  }) %>% 
    Reduce(function(x,y) dplyr::full_join(x,y, by = 'genes'), .)
  
  return(pb)
  
}

# batch effect removal 
