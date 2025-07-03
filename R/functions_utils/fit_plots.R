

plot_fit_statistics = function(x, what) {
  x %>% 
    dplyr::filter(term != '(Intercept)') %>% 
    ggplot(aes(.data[[what]])) + 
    geom_histogram(binwidth = 0.01) + 
    theme_bw()
}

plot_fit = function(data, 
                    x, 
                    y, 
                    facet,
                    color = NULL, 
                    add_fit_line = TRUE, 
                    color_palette = NULL
) {
  
  if (is.null(color)) {
    p = data %>% 
      ggplot(aes(x = .data[[x]], 
                 y = .data[[y]])) + 
               geom_point() +
               facet_wrap(~.data[[facet]], scales = 'free') +
               theme_bw()
  } else {
    p = data %>% 
      ggplot(aes(x = .data[[x]], 
                 y = .data[[y]], 
                 color = .data[[color]])) + 
      geom_point() +
      scale_fill_manual(values = color_palette) +
      facet_wrap(~.data[[facet]], scales = 'free') +
      theme_bw()
    }
  
  if(add_fit_line == TRUE) {
    p = p +
      geom_smooth(method = 'lm') +
      stat_poly_eq(use_label(c('eq', 'R2', 'P'))) 
  } 
  
  return(p)
}


plot_all_fit = function(x, genes, formula) {
  
  x = x %>% 
    dplyr::filter(hgnc_symbol %in% genes)
  
  lm_x = extract_lm_per_gene(x, formula) 
  
  pval = plot_fit_statistics(lm_x, 'p.value') + 
    ggtitle('Pvalues distribution')
  
  estimates = plot_fit_statistics(lm_x, 'estimate') + 
    ggtitle('Linear coefficient distribution')  
  
  r_squared = plot_fit_statistics(lm_x, 'r.squared') + 
    ggtitle('R2 distribution')
  
  plot_fit = plot_fit(x, 
                      x = 'value', 
                      y = 'mean_intensity', 
                      facet = 'hgnc_symbol', 
                      add_fit_line = TRUE) +
    xlab('RNA expression') + 
    ylab('Protein expression (mean of replicates)')
  
  list(
    fit = lm_x, 
    pvalues = pval, 
    estimate = estimates, 
    r_squared = r_squared, 
    plot_res = plot_fit
  )
}
