# plot utils

plot_expression_with_regression_lines = function(data, gene, fit, K_ref = 2) {
  
  alt_fit = fit[[gene]]$fit_alt %>% 
    broom::tidy() %>%
    mutate(gene = gene) %>% 
    select(gene, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    rename(intercept = `(Intercept)`)
  
  null_fit = fit[[gene]]$fit_null %>% 
    broom::tidy() %>%
    mutate(gene = gene) %>% 
    select(gene, estimate) %>%
    rename(intercept_null = estimate)
  
  df = data %>% 
    filter(hgnc_symbol == gene) %>% 
    filter(mutation_status == 'Wild-type') %>%
    group_by(hgnc_symbol) %>%
    full_join(., alt_fit, by = join_by('hgnc_symbol' == 'gene')) %>%
    full_join(., null_fit, by = join_by('hgnc_symbol' == 'gene')) %>%
    mutate(cna_log = log(tot_cna / K_ref))
    
  
  df %>%
    ggplot(aes(x = cna_log, y = expression)) +
    geom_point() +
    geom_abline(aes(intercept = intercept_null,
                    slope = 1),
                linetype = "dashed",
                size = 1.2,
                alpha = 0.8,
                color = "red") +
    geom_abline(aes(intercept = intercept,
                    slope = scaled_cna),
                # linetype = "dashed",
                size = 1.2,
                alpha = 0.8,
                color = "blue") +
    theme_bw() +
    labs(title = paste0('Sublinearity test (', gene, ')'),
         subtitle = "Red dashed: linear expectation (slope = 1)\nBlue: fitted slope",
         x = "log(CNA/2)", y = "log normalised expression")
    
  
}
