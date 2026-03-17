get_genes_by_group = function(full_omics, group) {
  genes_by_function = full_omics %>% 
    dplyr::filter(!!sym(group) != 'None', !is.na(!!sym(group))) %>% 
    group_by(pick(group)) %>% 
    group_split()
  names(genes_by_function) = lapply(genes_by_function, function(x) {x[,group] %>% unique}) %>% unlist
  genes_by_function = lapply(genes_by_function, function(x) {x$hgnc_symbol %>% unique})
  return(genes_by_function)
}

get_groups_coefficients = function(genes_by_function, coefficients_all) {
  res = lapply(genes_by_function, function(x) {
    coefficients_all %>% 
      tibble::rownames_to_column('gene') %>% 
      dplyr::filter(gene %in% x) %>% 
      tibble::column_to_rownames('gene') %>% 
      as.matrix()
  })
  names(res) = names(genes_by_function)
  return(res)
}

plot_heatmap_coefficients = function(coeffs_by_function, 
                                     ann 
                                     # ht_col
                                     ) {
  
  heatmaps = lapply(coeffs_by_function %>% names, function(x) {
    mat = coeffs_by_function[[x]]
    if(nrow(mat) > 0) {
      
      # colors = circlize::colorRamp2(c(min(mat)-0.2, 0, max(mat)+0.2), c('springgreen4', 'snow', 'dodgerblue4'))
      ht_colors = colorRamp2(c(min(mat)-0.1,0,max(mat)+0.1), c('blue', 'snow', 'red'))
      ht = Heatmap(mat, 
                   col = ht_colors,
                   bottom_annotation = ann, 
                   name = 'Predictor coefficients', 
                   show_column_names = F, 
                   column_title = x, 
                   heatmap_legend_param = list(direction = "horizontal", 
                                               at = c(min(mat)-0.1, 0, max(mat)+0.1),
                                               labels = c(round(min(mat), 1), "0", 
                                                          round(max(mat), 1))
                                               )
      )
      draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')}
    
  })
  heatmaps = Filter(Negate(is.null), heatmaps)
  return(heatmaps)
}

